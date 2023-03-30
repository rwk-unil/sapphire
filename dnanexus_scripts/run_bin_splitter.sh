#!/bin/bash

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

BIN_ID=""
OFNAME=""
SPLIT_SIZE="1000"
DESTINATION="phasing_rare"

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# Command line argument parsing from :
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
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
    -n|--split-size)
    SPLIT_SIZE="$2"
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

if [ -z "${BIN_ID}" ]
then
    echo "Specify an input binary file ID with --bin-id <actual_input_id>"
    echo "This file must be generated from the VCF passed"
    exit 1
fi


CHROMOSOME=$(basename $(dx describe --json "${BIN_ID}" | jq -r '.folder'))
echo "CHROMOSOME      = ${CHROMOSOME}"
BIN_FILENAME=$(dx_id_to_path "${BIN_ID}")
echo "BIN FILENAME    = ${BIN_FILENAME}"

if [ -z "${OFNAME}" ]
then
    OFNAME=$(basename "${BIN_FILENAME}")
fi

tag=bin_splitter_v1.2
echo "dx run with tag : ${tag}"

# The update_pp script will pull the git repo, rebuild all tools, and install them
# This allows to update the tools without having to rebuild the fat docker image
command="/usr/src/pp/Docker/update_pp.sh; time bin_splitter -b ${BIN_FILENAME} -o ${OFNAME} -n ${SPLIT_SIZE} ${VERBOSE}"

echo "Command : ${command}"
echo "Output file destination : ${DESTINATION}/het_extraction/${CHROMOSOME}/split"
echo "Optput file prefix : ${OFNAME}"
echo "Instance type : ${INSTANCE}"

ask_permission_to_launch

dx run swiss-army-knife -icmd="${command}" \
    ${COST_LIMIT_ARG} --name BinSplitter \
    -iimage_file="${DOCKER_IMAGE}" --tag "${tag}" \
    --destination "${DESTINATION}/het_extraction/${CHROMOSOME}/split" \
    --instance-type "${INSTANCE}" -y