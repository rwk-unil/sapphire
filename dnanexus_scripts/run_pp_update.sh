#!/bin/bash

if ! command -v realpath &> /dev/null
then
    realpath() {
        [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
    }
fi

# Get the path of this script
SCRIPTPATH=$(realpath  $(dirname "$0"))

VERBOSE=""
COST_LIMIT=""
INSTANCE="mem2_ssd1_v2_x2"

# Source common variables and functions
source "${SCRIPTPATH}/common.sh"

VCF_ID=""
BIN_ID=""
OFNAME=""
DESTINATION="phasing_rare"

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# Command line argument parsing from :
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
    -f|--vcf-id)
    VCF_ID="$2"
    shift # past argument
    shift # past value
    ;;
    -b|--bin-id)
    BIN_ID="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--output)
    OFNAME="$2"
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

if [ -z "${VCF_ID}" ]
then
    echo "Specify an input BCF ID with --vcf-id <actual_input_id>"
    echo "This is the original VCF/BCF file to be rephased"
    exit 1
fi
if [ -z "${BIN_ID}" ]
then
    echo "Specify an input binary file ID with --bin-id <actual_input_id>"
    echo "This file must be generated from the VCF passed"
    exit 1
fi

VCF_FILENAME=$(dx_id_to_path "${VCF_ID}")
CHROMOSOME=$(basename $(dx describe --json "${VCF_ID}" | jq -r '.folder'))
echo "VCF_FILENAME    = ${VCF_FILENAME}"
echo "CHROMOSOME      = ${CHROMOSOME}"
BIN_FILENAME=$(dx_id_to_path "${BIN_ID}")
echo "BIN FILENAME    = ${BIN_FILENAME}"

if [ -z "${VCF_FILENAME}" ]
then
    echo "Cannot get filename..."
    exit 1
fi

tag=pp_update_v1.2
echo "dx run with tag : ${tag}"

NEW_VCF_FILE="$(basename ${VCF_FILENAME})_rephased.bcf"
if [ -z "${OFNAME}" ]
then
    OFNAME="${NEW_VCF_FILE}"
fi

# The update_pp script will pull the git repo, rebuild all tools, and install them
# This allows to update the tools without having to rebuild the fat docker image
command="/usr/src/pp/Docker/update_pp.sh; time pp_update -f ${VCF_FILENAME} -o ${OFNAME} -b ${BIN_FILENAME} ${VERBOSE}"

echo "Command : ${command}"
echo "Output file destination : ${DESTINATION}/rephased/${CHROMOSOME}"
echo "Instance type : ${INSTANCE}"

ask_permission_to_launch

dx run swiss-army-knife -icmd="${command}; bcftools index ${OFNAME}" \
    ${COST_LIMIT_ARG} --name UpdatePhase \
    -iimage_file=docker/pp_rephase_v1.2.tar.gz --tag "${tag}" \
    --destination "${DESTINATION}/rephased/${CHROMOSOME}" \
    --instance-type "${INSTANCE}" -y