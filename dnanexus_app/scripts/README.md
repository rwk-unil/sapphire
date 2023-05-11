# DNANexus Applications Scripts

These are scripts to launch the DNANexus applets, the applets can be launched from the GUI or from scripts.

* run_pp_extract.sh `./run_pp_extract.sh -i file-XXXYYYZZZ111222333444555 -d phasing_rare/test_scripted/ --batch`
* run_bin_splitter.sh `./run_bin_splitter.sh -b file-AAABBBCCCDDDEEEFFFGGGHHH -d phasing_rare/test_scripted/split`
* For phase calling run the script in dnanexus_scripts `../../dnanexus_scripts/run_phase_caller_batch.sh --vcf-id file-BCFBCFBCFBCFBCFBCFBCFBCF --bin-id file-BBB000BBB000BBB000BBB000 --samples-id file-SSSAAASSSAAASSSAAASSSAAA --project-id 23193 --instance mem2_ssd1_v2_x4 --destination phasing_rare/test_scripted/phase_called`
* run_bin_merger.sh `./run_bin_merger.sh --bin-dir phasing_rare/test_scripted/phase_called -d phasing_rare/test_scripted/merged`
* run_pp_update.sh `./run_pp_update.sh --original-vcf-file file-XXXYYYZZZ111222333444555 --rephased-bin-file file-JJJKKKLLLMMMNNNOOOPPPQQQ -d phasing_rare/test_scripted/rephased/`