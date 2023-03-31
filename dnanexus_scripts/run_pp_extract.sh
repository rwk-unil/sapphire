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

FILENAME=""
INPUT_ID=""
DESTINATION="phasing_rare/"

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
    -i|--input_id)
    INPUT_ID="$2"
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

if [ -z "${INPUT_ID}" ]
then
    echo "Specify an input ID with --input_id, -i <actual_input_id>"
    exit 1
fi

FILENAME=$(dx describe --json "${INPUT_ID}" | jq -r '.name')
CHROMOSOME=$(basename $(dx describe --json "${INPUT_ID}" | jq -r '.folder'))
echo "FILENAME        = ${FILENAME}"
echo "CHROMOSOME      = ${CHROMOSOME}"

if [ -z "${FILENAME}" ]
then
    echo "Cannot get filename..."
    exit 1
fi

tag=pp_extract_v1
echo "dx run with tag : ${tag}"
command="time run_pp_extract -f ${FILENAME}"

echo "Output file destination : ${DESTINATION}/het_extraction/${CHROMOSOME}"
echo "Command : ${command}"

ask_permission_to_launch

dx run swiss-army-knife -icmd="ls -al; ${command}" \
    -iin="${INPUT_ID}" \
    ${COST_LIMIT_ARG} --name PPExtract \
    -iimage_file="${DOCKER_IMAGE}" -imount_inputs=true --tag "${tag}" \
    --destination "${DESTINATION}/het_extraction/${CHROMOSOME}/" \
    --instance-type ${INSTANCE} -y