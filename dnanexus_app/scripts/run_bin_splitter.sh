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
BIN_ID=""
OFNAME=""
SPLIT_SIZE="1000"
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

BIN_FILENAME=$(dx_id_to_path "${BIN_ID}")
echo "BIN FILENAME    = ${BIN_FILENAME}"

if [ -z "${OFNAME}" ]
then
    OFNAME=$(basename "${BIN_FILENAME}")
fi

command=(dx run applets/pp-toolkit/bin-splitter-applet \
    -ibinary_file_to_split=${BIN_ID} -isplit_size=${SPLIT_SIZE}\
    --destination "${DESTINATION}" \
    ${COST_LIMIT_ARG} --instance-type ${INSTANCE} -y \
    --name "PP-Toolkit : Step 2 - Bin splitter" \
    --tag "split")


echo "${command[@]}"

if [ "${BATCH}" = "yes" ]
then
    echo Batch mode - Launching without asking for permission
else
    ask_permission_to_launch
fi

"${command[@]}"
