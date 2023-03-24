#!/bin/bash

# Example : ./run_phase_caller.sh --vcf-id file-GQ1J7VjJ32gq4jVKXzF0gbXg --bin-id file-GQ1J7VjJ32gpJY68Q4K6B1kF --samples-id file-GQ1J7VjJ32gpj30XggFQ671V --samples-list-id file-GQ2xxQjJYyp15ZVQ6q86GPxf

if ! command -v realpath &> /dev/null
then
    realpath() {
        [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
    }
fi

# Get the path of this script
SCRIPTPATH=$(realpath  $(dirname "$0"))

INSTANCE="mem2_ssd1_v2_x4"

# Source common variables and functions
source "${SCRIPTPATH}/common.sh"

VCF_VAR_ID=""
BIN_FILE_ID=""
SAMPLE_FILE_ID=""
SAMPLE_LIST_ID=""
PROJECT_ID=23193
INPUT_ID=""
CRAM_PATH="/mnt/project/Bulk/Whole genome sequences/Whole genome CRAM files"
# 0 threads means "auto" (number of cores)
# Because the phase caller is bottlenecked by network access to CRAM files
# through the mounting point, setting approx 3 times more threads than cores
# results in an acceptable CPU usage (threads wait for data over network)
THREADS_ARG="-t 12"
DESTINATION="phasing_rare"
NEW_BIN_TAG="rephased"

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# Command line argument parsing from :
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
    --vcf-id)
    VCF_VAR_ID="$2"
    shift # past argument
    shift # past value
    ;;
    --bin-id)
    BIN_ID="$2"
    shift # past argument
    shift # past value
    ;;
    --samples-id)
    SAMPLE_FILE_ID="$2"
    shift # past argument
    shift # past value
    ;;
    --samples-list-id)
    SAMPLE_LIST_ID="$2"
    shift # past argument
    shift # past value
    ;;
    --project-id)
    PROJECT_ID="$2"
    shift # past argument
    shift # past value
    ;;
    --cram-path)
    CRAM_PATH="$2"
    shift # past argument
    shift # past value
    ;;
    --new-bin-tag)
    NEW_BIN_TAG="$2"
    shift # past argument
    shift # past value
    ;;
    --threads)
    THREADS_ARG="-t $2"
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

if [ -z "${VCF_VAR_ID}" ]
then
    echo "Specify an input VCF var ID with --vcf-id <actual_input_id>"
    echo "This file must be gnerated with the binary and sample file"
    exit 1
fi
if [ -z "${BIN_ID}" ]
then
    echo "Specify an input binary file ID with --bin-id <actual_input_id>"
    echo "This file must be generated with the VCF var and sample file"
    exit 1
fi
if [ -z "${SAMPLE_FILE_ID}" ]
then
    echo "Specify an input ID with --samples-id <actual_input_id>"
    echo "This file must be generated with the VCF var and binary file"
    exit 1
fi

VCF_FILENAME=$(dx_id_to_path "${VCF_VAR_ID}")
CHROMOSOME=$(basename $(dx describe --json "${VCF_VAR_ID}" | jq -r '.folder'))
echo "FILENAME        = ${VCF_FILENAME}"
echo "CHROMOSOME      = ${CHROMOSOME}"
BIN_FILENAME=$(dx_id_to_path "${BIN_ID}")
echo "BIN FILENAME    = ${BIN_FILENAME}"
SAMPLE_FILENAME=$(dx_id_to_path "${SAMPLE_FILE_ID}")
echo "SAMPLE FILENAME = ${SAMPLE_FILENAME}"

if ! [ -z "${SAMPLE_LIST_ID}" ]
then
    SAMPLE_LIST_FILENAME=$(dx_id_to_path "${SAMPLE_LIST_ID}")
    echo "SAMPLE LIST = ${SAMPLE_LIST_FILENAME}"
    INIIDSAMPLELIST="-iin=${SAMPLE_LIST_ID}"
else
    SAMPLE_LIST_FILENAME=""
fi

tag="phase_caller_v1.2"
echo "dx run with tag : ${tag}"

CRAM_PATH_ARG=""
if ! [ -z "${CRAM_PATH}" ]
then
    CRAM_PATH_ARG="--cram-path \"${CRAM_PATH}\""
fi
SAMPLE_FILENAME_ARG=""
if ! [ -z "${SAMPLE_FILENAME}" ]
then
    SAMPLE_FILENAME_ARG="-S \"${SAMPLE_FILENAME}\""
fi
SAMPLE_LIST_FILENAME_ARG=""
if ! [ -z "${SAMPLE_LIST_FILENAME}" ]
then
    SAMPLE_LIST_FILENAME_ARG="-l \"${SAMPLE_LIST_FILENAME}\""
fi

NEW_BINARY_FILE="$(basename ${BIN_FILENAME})_${NEW_BIN_TAG}.bin"
command="/usr/src/pp/Docker/update_pp.sh; cp ${BIN_FILENAME} ${NEW_BINARY_FILE}; time phase_caller -f ${VCF_FILENAME} -b ${NEW_BINARY_FILE} ${SAMPLE_FILENAME_ARG} -I ${PROJECT_ID} ${THREADS_ARG} ${SAMPLE_LIST_FILENAME_ARG} ${VERBOSE} ${CRAM_PATH_ARG}"

echo "Command : ${command}"
echo "Output file destination : ${DESTINATION}/phase_called/${CHROMOSOME}"
echo "Instance type : ${INSTANCE}"

ask_permission_to_launch

dx run swiss-army-knife -icmd="${command}" \
    ${COST_LIMIT_ARG} --name RephaseCaller \
    -iimage_file=docker/pp_rephase_v1.2.tar.gz --tag "${tag}" \
    --destination "${DESTINATION}/phase_called/${CHROMOSOME}" \
    --instance-type ${INSTANCE} -y