#!/bin/bash

# For documentation see doc/DNANexus.md

if ! command -v realpath &> /dev/null
then
    realpath() {
        [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
    }
fi

# Get the path of this script
SCRIPTPATH=$(realpath  $(dirname "$0"))
cd "${SCRIPTPATH}"

# Source common variables and functions
source "${SCRIPTPATH}/common.sh"

VCF_VAR_ID=""
VCF_IDX_ID=""
CHROMOSOME=""
DESTINATION="SAPPHIRE"

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    --vcf)
    VCF_VAR_ID="$2"
    shift # past argument
    shift # past value
    ;;
    --vcf-idx)
    VCF_IDX_ID="$2"
    shift
    shift
    ;;
    -d|--destination)
    DESTINATION="$2"
    shift
    shift
    ;;
    --chromosome)
    CHROMOSOME="$2"
    shift
    shift
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

if [ -z "${VCF_VAR_ID}" ]
then
    echo "Specify an input VCF/BCF file or ID with --vcf <vcf file or id>"
    exit 1
fi

if [ -z "${VCF_IDX_ID}" ]
then
    echo "Specify the index file or ID for the input VCF/BCF file with --vcf-idx <vcf index file or id>"
    exit 1
fi

if [ -z "${CHROMOSOME}" ]
then
    echo "Specify a chromosome with --chromosome <value>"
    exit 1
fi

full_path_filename=$(dx_id_to_dx_path_and_name "${VCF_VAR_ID}")
filename=$(dx_id_to_name "${VCF_VAR_ID}")

echo "Filename: ${filename}"
echo "Filename with path: ${full_path_filename}"

# I use two arrays instead of a dictionary because older versions of bash (<4.x) don't support it
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" \
             "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" \
             "chr20" "chr21" "chr22" "chrX")
stop_positions=("248999999" "242999999" "198999999" "190999999" "181999999" "170999999" \
               "159999999" "145999999" "138999999" "133999999" "135999999" "133999999" \
               "114999999" "107999999" "101999999" "90999999" "83999999" "80999999" \
               "58999999" "64999999" "46999999" "50999999" "156999999")

get_stop_pos() {
    local search_key=$1
    for i in "${!chromosomes[@]}"
    do
        if [[ "${chromosomes[$i]}" == "$search_key" ]]
        then
            echo "${stop_positions[$i]}"
            return
        fi
    done
    echo "$1 chromosome not found in array"
    exit 1
}

command_file="commands_${CHROMOSOME}.md"

echo "# SAPPHIRE commands for chromosome: ${CHROMOSOME}" > "${command_file}"
echo >> "${command_file}"
echo "\`\`\`shell" >> "${command_file}"
echo "# STEP 0" >> "${command_file}"
echo "./step0_split_bcf.sh --vcf-id $(path_to_dx_id "${VCF_VAR_ID}") --vcf-idx-id $(path_to_dx_id "${VCF_IDX_ID}") --chromosome ${CHROMOSOME} --stop-pos $(get_stop_pos "${CHROMOSOME}") --destination ${DESTINATION}" >> "${command_file}"
echo >> "${command_file}"
echo "# STEP 1" >> "${command_file}"
echo "./step1_prepare_bcfs.sh --step0-path ${DESTINATION}/SAPPHIRE_step0/${CHROMOSOME} --destination ${DESTINATION}" >> "${command_file}"
echo >> "${command_file}"
echo "# STEP 2" >> "${command_file}"
echo "./step2_prepare_variants.sh --step1-path ${DESTINATION}/SAPPHIRE_step1/${CHROMOSOME} --chromosome ${CHROMOSOME} --destination ${DESTINATION}" >> "${command_file}"
echo >> "${command_file}"
echo "# STEP 3" >> "${command_file}"
echo "./step3_pp_extract.sh --step0-path ${DESTINATION}/SAPPHIRE_step0/${CHROMOSOME} --step2-var ${DESTINATION}/SAPPHIRE_step2/${CHROMOSOME}/${filename}_${CHROMOSOME}.bcf --applet ${DESTINATION}/pp-extract-split-applet --destination ${DESTINATION}" >> "${command_file}"
echo >> "${command_file}"
echo "# STEP 4" >> "${command_file}"
echo "./step4_merge_regions_bin.sh --step3-path ${DESTINATION}/SAPPHIRE_step3/${CHROMOSOME} --docker ${DESTINATION}/pp_toolkit_v1.4.tar.gz --chromosome ${CHROMOSOME} --destination ${DESTINATION}" >> "${command_file}"
echo >> "${command_file}"
echo "# STEP 5" >> "${command_file}"
echo "./step5_split_binary.sh --binary-file ${DESTINATION}/SAPPHIRE_step4/${CHROMOSOME}/${filename}.bin --applet ${DESTINATION}/bin-splitter-applet --chromosome ${CHROMOSOME} --destination ${DESTINATION}" >> "${command_file}"
echo >> "${command_file}"
echo "# STEP 6a" >> "${command_file}"
echo "./step6a_generate_sample_list.sh --vcf \"${full_path_filename}\" --cram-list ${DESTINATION}/cram_paths.csv --chromosome ${CHROMOSOME} --destination ${DESTINATION}" >> "${command_file}"
echo >> "${command_file}"
echo "# STEP 6b" >> "${command_file}"
echo "./step6b_phase.sh --step5-path ${DESTINATION}/SAPPHIRE_step5/${CHROMOSOME} --step2-var ${DESTINATION}/SAPPHIRE_step2/${CHROMOSOME}/${filename}_${CHROMOSOME}.bcf --sample-list ${DESTINATION}/SAPPHIRE_data/${CHROMOSOME}/sample_list.csv --cram-path-file ${DESTINATION}/SAPPHIRE_data/${CHROMOSOME}/cram_paths_for_samples.csv --chromosome ${CHROMOSOME} --docker-image ${DESTINATION}/pp_toolkit_v1.4.tar.gz --destination ${DESTINATION}" >> "${command_file}"
echo >> "${command_file}"
echo "# STEP 7" >> "${command_file}"
echo "./step7_merge_binary_files.sh --step6-path ${DESTINATION}/SAPPHIRE_step6/${CHROMOSOME} --applet ${DESTINATION}/bin-merger-applet --chromosome ${CHROMOSOME} --destination ${DESTINATION}" >> "${command_file}"
echo >> "${command_file}"
echo "# STEP 8" >> "${command_file}"
echo "./step8_update_bcfs.sh --step1-path ${DESTINATION}/SAPPHIRE_step1/${CHROMOSOME} --step2-vars ${DESTINATION}/SAPPHIRE_step2/${CHROMOSOME}/${filename}_${CHROMOSOME}.bcf --step7-binary ${DESTINATION}/SAPPHIRE_step7/${CHROMOSOME}/${filename}.rephased.bin_sub_merged.bin --applet ${DESTINATION}/pp-update-applet --chromosome ${CHROMOSOME} --destination ${DESTINATION}" >> "${command_file}"
echo >> "${command_file}"
echo "# STEP 9" >> "${command_file}"
echo "./step9_concat_final_bcf.sh --step8-path ${DESTINATION}/SAPPHIRE_step8/${CHROMOSOME} --chromosome ${CHROMOSOME} --destination ${DESTINATION}" >> "${command_file}"
echo "\`\`\`" >> "${command_file}"
echo >> "${command_file}"


echo "Commands are available in file : ${command_file}"