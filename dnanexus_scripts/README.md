# DNANexus Scripts

This directory contains scripts to launch jobs on DNANexus.

These scripts require a Docker image (set path in `common.sh`). The Docker image can be built from the `Dockerfile_fat2` Dockerfile, see [README.md](../Docker/README.md) in Docker directory for more info.

## Scripts

It is recommended to use the applets rather than the Docker image for most tools except the phase caller.

1) run_pp_extract.sh
2) run_bin_splitter.sh
3) run_phase_caller_batch.sh **only script really necessary**
4) ~~run_bin_merger.sh~~ TODO
5) run_pp_update.sh
6) split_bcf.sh `./split_bcf.sh -f file-<bcf id> -i file-<index id> -c chr1 -m 247000000 -d destination/folder`