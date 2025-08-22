#!/bin/bash

# Example : ./run_phase_caller_batch.sh --vcf-id file-GV4GV78JGb8z5YQxFzjGXZ9x --bin-id file-GV4J09QJ7K6GkKj9Z216p36f --samples-id file-GV4GV6jJGb8V3yzkbXgYGfXz --project-id 23193 --instance mem2_ssd1_v2_x4 -d PhasePolishing/step3_phase_calling/chr3 --tag phase_caller_chr3

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
DESTINATION=""

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
    -d|--destination)
    DESTINATION="$2"
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
if [ -z "${DESTINATION}" ]
then
    echo "Please specify a destination folder with --destination <path>, different than the input file folder"
    exit 1
fi

VCF_FILENAME=$(dx_id_to_path "${VCF_VAR_ID}")
echo "FILENAME        = ${VCF_FILENAME}"
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

if [ -z "${TAG}" ]
then
    tag="phase_caller"
else
    tag="${TAG}"
fi
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

echo "Output file destination : ${DESTINATION}"
echo "Instance type : ${INSTANCE}"

echo "This will launch multiple jobs !"

LAUNCH_ALL=""
NUM_JOBS=$(dx ls $(dx_id_to_dx_path "${BIN_ID}") | wc -l)
echo Expected number of jobs : "${NUM_JOBS}"
for split_bin_file in $(dx ls $(dx_id_to_dx_path "${BIN_ID}"))
do
    INNER_BIN_FILENAME=$(dirname "${BIN_FILENAME}")/"${split_bin_file}"
    INNER_NEW_BINARY_FILE="$(basename ${INNER_BIN_FILENAME})"
    command="cp ${INNER_BIN_FILENAME} ${INNER_NEW_BINARY_FILE}; time phase_caller -f ${VCF_FILENAME} -b ${INNER_NEW_BINARY_FILE} ${SAMPLE_FILENAME_ARG} -I ${PROJECT_ID} ${THREADS_ARG} ${SAMPLE_LIST_FILENAME_ARG} ${VERBOSE} ${CRAM_PATH_ARG}"
    echo "${command}"
    if [ -z "${LAUNCH_ALL}" ]
    then
        ask_permission_to_launch_all
    fi
    if [ "${SKIP_JOB}" = "yes" ]
    then
        unset SKIP_JOB
    else
        dx run swiss-army-knife -icmd="/usr/src/pp/Docker/update_pp.sh; ${command}" \
            ${COST_LIMIT_ARG} --name RephaseCallerBatch \
            -iimage_file="${DOCKER_IMAGE}" --tag "${tag}" \
            --destination "${DESTINATION}" --priority high \
            --instance-type ${INSTANCE} -y
    fi
done

