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

INSTANCE="mem1_ssd1_v2_x36"

# Source common variables and functions
source "${SCRIPTPATH}/common.sh"

STEP1_PATH=""

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# Command line argument parsing from :
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
    --step1-path)
    STEP1_PATH="$2"
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

if [ -z "${STEP1_PATH}" ]
then
    echo "Provide the path of step1 for the chromosome --step1-path <path>"
    exit 1
fi

if [ -z "${CHROMOSOME}" ]
then
    echo "Specify a chromosome with --chromosome <value>"
    exit 1
fi

if [ -z "${TAG}" ]
then
    tag="step2_prepare_variants"
else
    tag="${TAG}"
fi
echo "dx run with tag : ${tag}"

echo "Instance type : ${INSTANCE}"

LAUNCH_ALL=""

array=($(dx ls --brief "${STEP1_PATH}" | sed 's/^/-iin=/'))



command="cp ./s2_script.sh s2.sh; chmod +x ./s2.sh; ./s2.sh; rm ./s2.sh"

echo "Command : ${command}"
echo "Will handle the files below :"
echo $(dx ls "${STEP1_PATH}")
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
    script=$(dx upload "${SCRIPTPATH}"/s2_script.sh --path ${DESTINATION}/s2_script.sh --brief)
    DESTINATION="${DESTINATION}/SAPPHIRE_step2/${CHROMOSOME}"
    dx run swiss-army-knife -icmd="${command}" \
        -iin=${script} \
        ${array[@]} \
        -imount_inputs=true \
        ${COST_LIMIT_ARG} --name "Step2: Prepare Variant BCF ${CHROMOSOME}" \
        --tag "${tag}" \
        --destination "${DESTINATION}" --priority normal \
        --instance-type ${INSTANCE} -y
fi

