#!/bin/bash

# ./run_bin_splitter_dir.sh --input_dir PhasePolishing/step1_extraction -d PhasePolishing/step2_splitting

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
    -i|--input_dir)
    INPUT_DIR="$2"
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

if [ -z "${INPUT_DIR}" ]
then
    echo "Specify an input directory with --input_dir <actual_input_dir>"
    exit 1
fi

dx select
dirs_to_process=$(dx ls "${INPUT_DIR}")

for dir in ${dirs_to_process}
do
    echo "${INPUT_DIR}/${dir}"
    file=$(dx ls "${INPUT_DIR}/${dir}" | grep -e ".*.bin")
    echo "${file}"
    if ! [ -z "${file}" ]
    then
        BIN_ID=$(dx describe --json "${INPUT_DIR}/${dir}/${file}" | jq -r .id)
        #echo "${BIN_ID}"
        ${SCRIPTPATH}/run_bin_splitter.sh -b "${BIN_ID}" -d "${DESTINATION}/${dir}" --split-size "${SPLIT_SIZE}" ${BATCH_ARG}
    fi
done
