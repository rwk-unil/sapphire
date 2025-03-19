# Running SAPPHIRE on DNANexus

The different scripts for the pipeline can be run from a regular Linux/Mac or Windows with WSL, they can also be run from a ttyd instance on DNANexus.

Why is this pipeline so complicated ? Mainly because of scalability, to be able to rephase hundreds of thousands of samples on possibly billions of variants sites by using their WGS data, the work has to be split at many stages in different ways.

## Preparations

- Install jq e.g., `sudo apt install jq`
- Install the DNANexus CLI tool https://documentation.dnanexus.com/getting-started/cli-quickstart (already installed on the ttyd instance)
- Clone the SAPPHIRE git and go in the directory `git clone https://github.com/rwk-unil/sapphire.git && cd sapphire/`

If something else is missing please open an issue and I'll add it to the documentation.

Please note that all scripts below will use the DNANexus CLI tool `dx`.
For the scripts to work the dx CLI tool should have the correct project selected and the user should be logged in. (already done if running from ttyd instance)

Please run `dx login` to login

Please run `dx select` to select the correct project before running the scripts

In this guide we will use a destination folder named `SAPPHIRE` within our project. If you want to change this, update the `--destionation` arguments for all scripts below, as steps depend on the previous steps, don't move files, changes names, destinations etc. Please don't use spaces in the destination path.

Warning: For all steps check that the jobs did finish successfully before running the next step. If some jobs fail relaunch them. If there is a persisting issue, solve it. If you need help, open an issue on the github page.

### Create and run the PP-Toolkit builder

This will create the tools we need on DNANexus

```shell
# Create a directory for SAPPHIRE (we will use it below)
dx mkdir SAPPHIRE
# Create and run the pp-toolkit-builder
dx run $(dx build dnanexus_app/pp-toolkit-builder --destination SAPPHIRE/ | jq -r ".id") --destination SAPPHIRE/ -y
# It will ask for a version, you can input anything e.g., the date YYMMDD
```

### Download the docker image

For some steps a prebuilt docker image is used, this image for example, holds all the reference sequences (human genome assembly) for the CRAM files so that they don't need to be dowloaded during job execution.

```shell
./upload_docker_image.sh --destination SAPPHIRE
```

Will ask to launch a job that downloads the docker image and uploads it in the DNANexus project in the destination folder.

Note: The Docker image has the v1.4 tools but it has the script `/usr/src/pp/Docker/update_pp.sh` that allows it to update itself to the latest version (run from within the Docker). This is used in the scripts below to update the tools automatically. If you use or made a newer docker image, the path is `/usr/src/sapphire/Docker/update_pp.sh`.

### Generate the file with CRAM paths and IDs

Here you can choose which CRAM file release you will be using by providing a UK Biobank Release CRAM path.

(A the path locations should be the subdirectories that starts with the two first numbers of the sample IDs).

```shell
/generate_cram_paths.sh --cram-path "Bulk/DRAGEN WGS/Whole genome CRAM files (DRAGEN) [500k release]" --destination SAPPHIRE
```

Alternatives could be `"Bulk/GATK and GraphTyper WGS/Whole genome GATK CRAM files and indices [500k release]/"` for example.

### Notes

If jobs fail, which could happen when instances don't have enough storage for example, it is possible to change the instance type with `--instance` e.g., `--instance mem2_ssd3_v2_x2`

The instances have been selected to be the cheapest possible to run the phase polishing for the UK Biobank 200k release, for the 500k release the storage of some instances may not be enough.

Check the rate card : https://20779781.fs1.hubspotusercontent-na1.net/hubfs/20779781/Product%20Team%20Folder/Rate%20Cards/BiobankResearchAnalysisPlatform_Rate%20Card_Current.pdf

For price and resource information, and select another instance if needed.

# Phase polishing of a population chromosome

The input VCF/BCF file should be phased (preferably with SHAPEIT5) and indexed with `bcftools index`

Let's assume we will work on the following file :

`Bulk/Previous WGS releases/GATK and GraphTyper WGS/SHAPEIT Phased VCFs/ukb20279_c22_b0_v1.vcf.gz`

Which has the associated index :

`Bulk/Previous WGS releases/GATK and GraphTyper WGS/SHAPEIT Phased VCFs/ukb20279_c22_b0_v1.vcf.gz.tbi`

Note: The scripts below will ask to launch the jobs, you can choose `y` to launch a job, `n` to skip a job, `a` to launch all or `e` to exit.
This makes it possible to double check the job parameters before launch.

## Step0 Split input BCF/VCF file into overlapping chunks

Start by getting the file IDs for both files :

```shell
dx describe "Bulk/Previous WGS releases/GATK and GraphTyper WGS/SHAPEIT Phased VCFs/ukb20279_c22_b0_v1.vcf.gz"
# ID                                file-GX8z0P8JykJpF30X7xXKfxzq
dx describe "Bulk/Previous WGS releases/GATK and GraphTyper WGS/SHAPEIT Phased VCFs/ukb20279_c22_b0_v1.vcf.gz.tbi"
# ID                                file-GX8yvP8JykJf12qPjb18YgqX
```

The split BCF command takes the two file IDs, the chromosome (as in the input VCF/BCF), and stop position (should include the last record of the original VCF/BCF).

```shell
./step0_split_bcf.sh --vcf-id file-GX8z0P8JykJpF30X7xXKfxzq --vcf-idx-id file-GX8yvP8JykJf12qPjb18YgqX \
        --chromosome chr22 --stop-pos 50999999 --destination SAPPHIRE
```

This script will generate files and directories in the `--destination` folder here `SAPPHIRE`

## Step1 Generate the non-overlapping chunks

We need to prepare non-overlapping chunks that will serve as the basis to reconstruct the final BCF file with BCFTools concat.
These files will also be used to construct the BCF with variants only (it could be generated from the input file, but it would be slower).

You'll need to pass the path of step1 for the chromosome, the above jobs generated the files in `SAPPHIRE/SAPPHIRE_step0/chr22`

Also set the destination

```shell
./step1_prepare_bcfs.sh --step0-path SAPPHIRE/SAPPHIRE_step0/chr22 --destination SAPPHIRE
```

## Step2 Prepare the variant list for SAPPHIRE

Preparing the variant file is done with the following command, note the chromosome and destination should be set correctly.

```shell
./step2_prepare_variants.sh --step1-path SAPPHIRE/SAPPHIRE_step1/chr22 --chromosome chr22 --destination SAPPHIRE
```

Note: The actual commands for this step are in `s2_script.sh`, it is parallelized to the number of CPU cores of the instance (by default x36)

## Step3 Extract low-PP / low-MAF variants

For the moment the script provided only allows for low-PP extraction (SHAPEIT5 phased file)

Extraction requires the overlapping BCF files from step0 and the complete variant list from step2
And the pp-extract-split-applet, which is generated by the pp-toolkit-builder (see preparations above)

- You'll need to provide the path to the files of step0
- You'll need to provide the path or ID of the file generated in step2
- You'll need to provide the path or ID of the extraction applet `pp-extract-split-applet` **Don't use another applet**

```shell
./step3_pp_extract.sh --step0-path SAPPHIRE/SAPPHIRE_step0/chr22 \
        --step2-var SAPPHIRE/SAPPHIRE_step2/chr22/ukb20279_c22_b0_v1.vcf.gz_chr22.bcf \
        --applet SAPPHIRE/pp-extract-split-applet --destination SAPPHIRE
```

**TODO:** Need to add option for MAF to pp-extract-split-applet, for non SHAPEIT5 phased files.

## Step4 Merge the extracted regions binary files

This merges the extracted variants from the overlapping chunks to a single binary file at chromosome level, this is done so that the CRAMs will be accessed for the whole chromosome and not in chunks.

Provide the path of step3 and the docker image path (or ID)

```shell
./step4_merge_regions_bin.sh --step3-path SAPPHIRE/SAPPHIRE_step3/chr22 --docker SAPPHIRE/pp_toolkit_v1.4.tar.gz --chromosome chr22 --destination SAPPHIRE
```

## Step5 Split the binary file into batches of samples

Because for each sample a CRAM file is accessed this should be split into smaller batches, this step splits the population level binary files into smaller batches of samples (1,000 samples per batch).

```shell
./step5_split_binary.sh --binary-file SAPPHIRE/SAPPHIRE_step4/chr22/ukb20279_c22_b0_v1.vcf.gz.bin --applet SAPPHIRE/bin-splitter-applet --chromosome chr22 --destination SAPPHIRE
```

## Step6 Run the SAPPHIRE phase caller

The phase caller requires to be launched on each batch of 1000 samples (each splitted binary file), because the phase caller accesses the CRAM files of all samples, they should be mounted in the phase caller instance.

Before we launch the phase caller we need to prepare a file that links the samples in the VCF to their CRAM file.

```shell
./step6a_generate_sample_list.sh --vcf "Bulk/Previous WGS releases/GATK and GraphTyper WGS/SHAPEIT Phased VCFs/ukb20279_c22_b0_v1.vcf.gz" --cram-list SAPPHIRE/cram_paths.csv --chromosome chr22 --destination SAPPHIRE
```

To do the actual phasing provide the path of step5, the variant file ID generated in step2, the sample list generated above, the cram path file generated during the initial preparations, the docker image uploaded in the preparations, and the destination.

```shell
./step6b_phase.sh --step5-path SAPPHIRE/SAPPHIRE_step5/chr22 \
        --step2-var SAPPHIRE/SAPPHIRE_step2/chr22/ukb20279_c22_b0_v1.vcf.gz_chr22.bcf \
        --sample-list SAPPHIRE/SAPPHIRE_data/chr22/sample_list.csv \
        --cram-path-file SAPPHIRE/SAPPHIRE_data/chr22/cram_paths_for_samples.csv \
        --chromosome chr22 --docker-image SAPPHIRE/pp_toolkit_v1.4.tar.gz --destination SAPPHIRE
```

Because for each batch more than 2,000 files need to be mapped (1,000 CRAMs, 1,000 indices, and the other files) launching each job will takes about 30 seconds (for DNANexus to respond that the job was submitted). Therefore the script will launch all jobs in parallel and wait for the `dx run` commands to finish at the end., you can ask the script to launch them all with `a` when asked and wait for all jobs to be submitted.

For some reason the `dx run` command sometimes get stuck after submitting a job, if it is stuck for more than a minute or so, and you didn't launch all jobs with `a` when asked, you can interrupt it with ctrl-c and the job should still be submitted, it is just that the response from `dx run` didn't yet arrive.

The phasing itself should take 1-3h per 1000 sample job.

## Merge the binary files

Check that all rephasing jobs finished successfully, if not solve the issue, and relaunch the ones that failed.

The rephased binary files must now be merged as they were split into sample batches for rephasing.

Provide the path to the previous step, the path to the applet, the chromosome, and the destination folder.

```shell
./step7_merge_binary_files.sh --step6-path SAPPHIRE/SAPPHIRE_step6/chr22 \
        --applet SAPPHIRE/bin-merger-applet --chromosome chr22 --destination SAPPHIRE
```

This step should be fairly quick (less than 10 minutes)

## Update the non-overlapping chunks

Because updating a whole BCF for a whole chromosome for hundreds of thousands of samples takes a long time here we will update the region split BCF files generated in **step 1**. **Do not update the overlapping BCFs of step 0 !**

- This script will need the step1-path for the non-overlapping BCF files
- It will need the variant file from step2
- It will need the merged binary file from step7
- It will need the pp-update applet (generated with pp-toolkit-builder above)
- It will need the chromosome

```shell
./step8_update_bcfs.sh --step1-path SAPPHIRE/SAPPHIRE_step1/chr22 \
        --step2-vars SAPPHIRE/SAPPHIRE_step2/chr22/ukb20279_c22_b0_v1.vcf.gz_chr22.bcf \
        --step7-binary SAPPHIRE/SAPPHIRE_step7/chr22/ukb20279_c22_b0_v1.vcf.gz.rephased.bin_sub_merged.bin \
        --applet SAPPHIRE/pp-update-applet --chromosome chr22 --destination SAPPHIRE
```

## Concatenate to generate the final BCF file

Finally we will concatenate the updated non-overlapping BCF files to generate the final whole-chromosome BCF file.

Note: This step requires an instance with enough storage for the final BCF file, it should be roughly the same size as the original file used as input to the SAPPHIRE pipeline. By default `mem2_ssd2_v2_x2` is used which has 160~GB of storage.

If more is needed chose one from : https://20779781.fs1.hubspotusercontent-na1.net/hubfs/20779781/Product%20Team%20Folder/Rate%20Cards/BiobankResearchAnalysisPlatform_Rate%20Card_Current.pdf

And use `--instance` as argument to the script to specify another instance, e.g., `--instance mem2_ssd2_v2_x16`

```shell
./step9_concat_final_bcf.sh --step8-path SAPPHIRE/SAPPHIRE_step8/chr22 \
        --chromosome chr22 --destination SAPPHIRE
```

Note: If for some reason you prefer a `vcf.gz` file, it is preferable to convert the `.bcf` files from step8 to `.vcf.gz`, index them, and naive concat them (sorted by chromosomal region) to produce the final `vcf.gz` rather than converting the `.bcf` from step9.

# You Are Done !

Congratulations, it's time to celebrate !