#!/bin/bash

# For documentation see doc/DNANexus.md

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

BINARY=""
APPLET=""
SPLIT_SIZE=1000

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# Command line argument parsing from :
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
    -b|--binary-file)
    BINARY="$2"
    shift # past argument
    shift # past value
    ;;
    --applet)
    APPLET_ID="$2"
    shift
    shift
    ;;
    -d|--destination)
    DESTINATION="$2"
    shift
    shift
    ;;
    --chromosome)
    CHROMOSOME="$2"
    shift
    shift
    ;;
    --split-size)
    SPLIT_SIZE="$2"
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

if [ -z "${BINARY}" ]
then
    echo "Specify the binary file path with --binary-file <path to file>"
    exit 1
fi

if [ -z "${APPLET_ID}" ]
then
    echo "Please provide the ID or path of bin-splitter-applet with --applet <value>"
    exit 1
fi

if [ -z "${CHROMOSOME}" ]
then
    echo "Specify the chromosome with --chromosome <value>"
    exit 1
fi

if [ -z "${TAG}" ]
then
    tag="step5_split"
else
    tag="${TAG}"
fi
echo "dx run with tag : ${tag}"

echo "Instance type : ${INSTANCE}"

LAUNCH_ALL=""

position=0


inner_dest="${DESTINATION}/SAPPHIRE_step5/${CHROMOSOME}"

echo "Destination: ${inner_dest}"

if [ -z "${LAUNCH_ALL}" ]
then
    ask_permission_to_launch_all
fi
if [ "${SKIP_JOB}" = "yes" ]
then
    unset SKIP_JOB
else
    dx run "${APPLET_ID}" -ibinary_file_to_split="${BINARY}" -isplit_size="${SPLIT_SIZE}" \
        ${COST_LIMIT_ARG} --name "Step5: Binary split ${CHROMOSOME}" \
        --tag "${tag}" \
        --destination "${inner_dest}" --priority normal \
        --instance-type ${INSTANCE} -y
fi
