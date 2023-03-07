#!/bin/bash

# Example : ./run_phase_caller.sh --vcf-id file-GQ1J7VjJ32gq4jVKXzF0gbXg --bin-id file-GQ1J7VjJ32gpJY68Q4K6B1kF --samples-id file-GQ1J7VjJ32gpj30XggFQ671V --samples-list-id file-GQ2xxQjJYyp15ZVQ6q86GPxf

if ! command -v realpath &> /dev/null
then
    realpath() {
        [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
    }
fi

if ! command -v jq &> /dev/null
then
    echo "Please install jq"
    exit 1
fi

dx_id_to_path () {
    FFNAME=$(dx describe --json "$1" | jq -r '.name')
    FFPATH=$(dx describe --json "$1" | jq -r '.folder')
    echo "/mnt/project/${FFPATH}/${FFNAME}"
}

# Get the path of this script
SCRIPTPATH=$(realpath  $(dirname "$0"))

VCF_VAR_ID=""
BIN_FILE_ID=""
SAMPLE_FILE_ID=""
SAMPLE_LIST_ID=""
PROJECT_ID=23193
INPUT_ID=""
VERBOSE=""
COST_LIMIT=""
INSTANCE="mem3_ssd2_v2_x8"
CRAM_PATH="/mnt/project/Bulk/Whole genome sequences/Whole genome CRAM files"
# 0 threads means "auto" (number of cores)
THREADS_ARG="-t 0"

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# Command line argument parsing from :
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
    --verbose)
    VERBOSE="-v"
    shift # no value attached
    ;;
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
    --cost-limit)
    COST_LIMIT="$2"
    shift # past argument
    shift # past value
    ;;
    --instance)
    INSTANCE="$2"
    shift # past argument
    shift # past value
    ;;
    --cram-path)
    CRAM_PATH="$2"
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
    SAMPLE_LIST_FILENAME="${SAMPLE_FILENAME}"
fi

tag="phase_caller_v1.1"
echo "dx run with tag : ${tag}"

CRAM_PATH_ARG=""
if ! [ -z "${CRAM_PATH}" ]
then
    CRAM_PATH_ARG="--cram-path \"${CRAM_PATH}\""
fi

NEW_BINARY_FILE="$(basename ${BIN_FILENAME})_rephased.bin"
command="time phase_caller -f ${VCF_FILENAME} -b ${NEW_BINARY_FILE} -S ${SAMPLE_FILENAME} -I ${PROJECT_ID} ${THREADS_ARG} -l ${SAMPLE_LIST_FILENAME} ${VERBOSE} ${CRAM_PATH_ARG}"

echo "Command : ${command}"
echo "Instance type : ${INSTANCE}"

while true; do
    read -p "Do you want to launch on DNANexus? [y/n]" yn
    case $yn in 
        y)
        echo "Launching !";
        break
        ;;
        n)
        echo "exiting...";
        exit
        ;;
        *)
        echo "unexpected input"
        ;;
    esac
done

COST_LIMIT_ARG=""
if ! [ -z "${COST_LIMIT}" ]
then
    COST_LIMIT_ARG="--cost-limit ${COST_LIMIT}"
fi

dx run swiss-army-knife -icmd="ls -al; cp ${BIN_FILENAME} ${NEW_BINARY_FILE}; ${command}" \
    ${COST_LIMIT_ARG} --name RephaseCaller \
    -iimage_file=docker/pp_extract_v1.1.tar.gz --tag "${tag}" \
    --destination phasing_rare/validation_interm_results/ \
    --instance-type ${INSTANCE} -y