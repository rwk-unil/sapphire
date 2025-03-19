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
STEP6_PATH=""

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# Command line argument parsing from :
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
    --step6-path)
    STEP6_PATH="$2"
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
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

if [ -z "${STEP6_PATH}" ]
then
    echo "Please specify step6 path at the chromosome level with --step6-path <path>"
    exit 1
fi

if [ -z "${APPLET_ID}" ]
then
    echo "Please provide the ID or path of bin-merger-applet with --applet <value>"
    exit 1
fi

if [ -z "${CHROMOSOME}" ]
then
    echo "Please specify the chromosome with --chromosome <value>"
    exit 1
fi

if [ -z "${TAG}" ]
then
    tag="step7_merge"
else
    tag="${TAG}"
fi
echo "dx run with tag : ${tag}"

echo "Instance type : ${INSTANCE}"

LAUNCH_ALL=""

array=()
for filename in $(dx ls "${STEP6_PATH}")
do
    # Skip if file is not binary
    if [[ $filename != *.bin_sub_* ]]
    then
        echo "File ${filename} does not end in .bin_sub_*, ignoring..."
        continue
    fi

    file_with_path="${STEP6_PATH}"/"${filename}"
    array+=("-isplitted_binary_files=${file_with_path}")
done

inner_dest="${DESTINATION}/SAPPHIRE_step7/${CHROMOSOME}"

echo "Destination: ${inner_dest}"

if [ -z "${LAUNCH_ALL}" ]
then
    ask_permission_to_launch_all
fi
if [ "${SKIP_JOB}" = "yes" ]
then
    unset SKIP_JOB
else
    dx run "${APPLET_ID}" ${array[@]} \
        ${COST_LIMIT_ARG} --name "Step7: Binary merge ${CHROMOSOME}" \
        --tag "${tag}" \
        --destination "${inner_dest}" --priority normal \
        --instance-type ${INSTANCE} -y
fi
