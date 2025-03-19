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

INSTANCE="mem2_ssd2_v2_x2"

# Source common variables and functions
source "${SCRIPTPATH}/common.sh"

STEP8_PATH=""

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# Command line argument parsing from :
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
    --step8-path)
    STEP8_PATH="$2"
    shift # past argument
    shift # past value
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

if [ -z "${STEP8_PATH}" ]
then
    echo "Provide the path of step8 for the chromosome --step8-path <path>"
    exit 1
fi

if [ -z "${CHROMOSOME}" ]
then
    echo "Specify a chromosome with --chromosome <value>"
    exit 1
fi

if [ -z "${TAG}" ]
then
    tag="step9_final_bcf"
else
    tag="${TAG}"
fi
echo "dx run with tag : ${tag}"

echo "Instance type : ${INSTANCE}"

LAUNCH_ALL=""

# Generate array of files to mount
array=($(dx ls --brief "${STEP8_PATH}" | sed 's/^/-iin=/'))

command="cp ./s9_script.sh s9.sh; chmod +x ./s9.sh; ./s9.sh; rm ./s9.sh"

echo "Command : ${command}"
echo "Will handle the files below :"
echo $(dx ls "${STEP8_PATH}")
if [ -z "${LAUNCH_ALL}" ]
then
    ask_permission_to_launch_all
fi
if [ "${SKIP_JOB}" = "yes" ]
then
    unset SKIP_JOB
else
    if [ ! -z "${DESTINATION}" ]
    then
        dx mkdir -p "${DESTINATION}"
    fi
    script=$(dx upload "${SCRIPTPATH}"/s9_script.sh --path ${DESTINATION}/s9_script.sh --brief)
    DESTINATION="${DESTINATION}/SAPPHIRE_step9/${CHROMOSOME}"
    dx run swiss-army-knife -icmd="${command}" \
        -iin=${script} \
        ${array[@]} \
        -imount_inputs=true \
        ${COST_LIMIT_ARG} --name "Step9: Concat Final BCF" \
        --tag "${tag}" \
        --destination "${DESTINATION}" --priority normal \
        --instance-type ${INSTANCE} -y
fi
