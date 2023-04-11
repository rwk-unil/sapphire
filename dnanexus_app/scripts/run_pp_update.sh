#!/bin/bash

if ! command -v realpath &> /dev/null
then
    realpath() {
        [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
    }
fi

# Get the path of this script
SCRIPTPATH=$(realpath $(dirname "$0"))

INSTANCE="mem2_ssd1_v2_x2"
INPUT_ID=""
DESTINATION="/"

# Source common variables and functions
source "${SCRIPTPATH}/common.sh"

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# Command line argument parsing from :
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
    --original-vcf-file)
    ORIGINAL_VCF_FILE="$2"
    shift # past argument
    shift # past value
    ;;
    --rephased-bin-file)
    REPHASED_BIN_FILE="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

if [ -z "${ORIGINAL_VCF_FILE}" ]
then
    echo "Specify an input ID with --original-vcf-file <file_id>"
    exit 1
fi

if [ -z "${REPHASED_BIN_FILE}" ]
then
    echo "Specify an input ID with --rephased-bin-file <file_id>"
    exit 1
fi

echo "Output file destination : ${DESTINATION}"

# Use an array to keep parameters whole (e.g., paths with spaces)
command=(dx run applets/pp-toolkit/pp-update-applet \
    -ioriginal_vcf_file=${ORIGINAL_VCF_FILE} \
    -irephased_binary_file=${REPHASED_BIN_FILE} \
    --destination "${DESTINATION}" \
    ${COST_LIMIT_ARG} --instance-type ${INSTANCE} -y \
    --name "PP-Toolkit : Step 5 - Update" \
    --tag "update")

echo "${command[@]}"

if [ "${BATCH}" = "yes" ]
then
    echo Batch mode - Launching without asking for permission
else
    ask_permission_to_launch
fi

"${command[@]}"
