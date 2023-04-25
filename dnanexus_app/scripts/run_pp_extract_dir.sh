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
INPUT_DIR=""
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
#    -f|--filename)
#    FILENAME="$2"
#    shift # past argument
#    shift # past value
#    ;;
    -i|--input_dir)
    INPUT_DIR="$2"
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
    echo "Specify an input DIR with --input_dir, -i <actual_input_dir>"
    exit 1
fi

if [ -z "${DESTINATION}" ]
then
    echo "Specify an output destination --destination, -d <destination_path>"
    exit 1
fi

dx select
FILES_TO_PROCESS=$(dx ls "${INPUT_DIR}"/*.bcf)

for file in ${FILES_TO_PROCESS}
do
    echo "${file}"
    FILE_ID=$(dx describe --json "${INPUT_DIR}/${file}" | jq -r '.id')
    echo "${FILE_ID}"
    FILE_SIZE=$(dx describe --json ${FILE_ID} | jq -r '.size')
    echo "Input file size : $(numfmt --to=iec-i --suffix=B --format="%9.2f" ${FILE_SIZE})"

    # Choose the cheapest instance based on file size
    if (( FILE_SIZE > 135000000000 ))
    then
        # Handles up to 279 GB
        INSTANCE="mem2_ssd1_v2_x8"
    elif (( FILE_SIZE > 65000000000 ))
    then
        # Handles up to 139 GB
        INSTANCE="mem2_ssd1_v2_x4"
    else
        # Handles up to 69 GB
        INSTANCE="mem2_ssd1_v2_x2"
    fi
    subdir=$(echo "${file}" | sed 's/.*\(chr[0-9]*\).*/\1/')
    if [ -z "${subdir}" ]
    then
        read -p "Please enter subdir :" subdir
    fi
    ${SCRIPTPATH}/run_pp_extract.sh -i "${FILE_ID}" -d "${DESTINATION}/${subdir}" --instance "${INSTANCE}" ${BATCH_ARG}
done

exit 0

