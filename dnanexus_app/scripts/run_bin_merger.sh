#!/bin/bash

if ! command -v realpath &> /dev/null
then
    realpath() {
        [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
    }
fi

# Get the path of this script
SCRIPTPATH=$(realpath  $(dirname "$0"))
DESTINATION="/"
INSTANCE="mem2_ssd1_v2_x2"

# Source common variables and functions
source "${SCRIPTPATH}/common.sh"

BIN_FILE_DIR=""

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# Command line argument parsing from :
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
    --bin-dir)
    BIN_DIR="$2"
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

if [ -z "${BIN_DIR}" ]
then
    echo "Specify an input binary file directory with --bin-dir <actual_input_dir>"
    exit 1
fi

command=(dx run \
    --destination "${DESTINATION}" \
    ${COST_LIMIT_ARG} --instance-type ${INSTANCE} -y \
    --name "PP-Toolkit : Step 4 - Bin merger" \
    --tag "merge" \
    applets/pp-toolkit/bin-merger-applet)

for split_bin_file in $(dx ls ${BIN_DIR})
do
    echo Found file "${split_bin_file}"
    local_id=$(dx describe "${BIN_DIR}/${split_bin_file}" --json | jq -r '.id')
    echo With ID : "${local_id}"
    command+=("-isplitted_binary_files=${local_id}")
done

echo "${command[@]}"

if [ "${BATCH}" = "yes" ]
then
    echo Batch mode - Launching without asking for permission
else
    ask_permission_to_launch
fi

"${command[@]}"
