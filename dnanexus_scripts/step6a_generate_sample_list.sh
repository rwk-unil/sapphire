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

INSTANCE="mem2_ssd1_v2_x2"

# Source common variables and functions
source "${SCRIPTPATH}/common.sh"

VCF_FILE=""
CRAM_LIST=""

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# Command line argument parsing from :
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
    --vcf)
    VCF_FILE="$2"
    shift # past argument
    shift # past value
    ;;
    --cram-list)
    CRAM_LIST="$2"
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

if [ -z "${VCF_FILE}" ]
then
    echo "Please provide VCF/BCF filename with --vcf <path to file>"
    exit 1
fi

if [ -z "${CRAM_LIST}" ]
then
    echo "Please specify the cram list file with --cram_list <path to file>"
    exit 1
fi

if [ -z "${CHROMOSOME}" ]
then
    echo "Please specify the chromosome with --chromosome <value>"
    exit 1
fi

inner_dest="${DESTINATION}/SAPPHIRE_data/${CHROMOSOME}"
echo "Output file destination : ${inner_dest}"

filename="$(basename "${VCF_FILE}")"
cram_filename="$(basename "${CRAM_LIST}")"

if [ -z "${TAG}" ]
then
    tag="step6a_sample_list"
else
    tag="${TAG}"
fi


# Generate sample list with bcftool query
# And generate list with information per sample in the list, sorted as in the list, with CRAM info
# This second file will be used to launch the phase caller jobs
# And create a file for the phase caller itself
# Remove the sample list (it was just used to get the correct samples and order from the VCF/BCF)
command="bcftools query --list-samples ${filename} > sample_list.txt && \
         awk -F',' 'NR==FNR {map[\$1] = \$0; next} \$1 in map {print map[\$1]} !(\$1 in map) {print \$0 \",,,,,,,,,,\"}' ${cram_filename} sample_list.txt | awk '{print NR-1 \",\" \$0 }' > cram_paths_for_samples.csv && \
         cut -d, -f1-3 cram_paths_for_samples.csv > sample_list.csv && \
         rm sample_list.txt"

echo "Command : ${command}"
if [ -z "${LAUNCH_ALL}" ]
then
ask_permission_to_launch_all
fi
if [ "${SKIP_JOB}" = "yes" ]
then
unset SKIP_JOB
else
dx run swiss-army-knife -icmd="${command}" \
        -iin="${VCF_FILE}" -iin="${CRAM_LIST}" \
        -imount_inputs=true \
        ${COST_LIMIT_ARG} --name "Step6a: Prepare sample list ${CHROMOSOME}" \
        --tag "${tag}" \
        --destination "${inner_dest}" --priority normal \
        --instance-type ${INSTANCE} -y
fi