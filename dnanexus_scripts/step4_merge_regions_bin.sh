#!/bin/bash

# Example : ./run_phase_caller_batch.sh --vcf-id file-GV4GV78JGb8z5YQxFzjGXZ9x --bin-id file-GV4J09QJ7K6GkKj9Z216p36f --samples-id file-GV4GV6jJGb8V3yzkbXgYGfXz --project-id 23193 --instance mem2_ssd1_v2_x4 -d PhasePolishing/step3_phase_calling/chr3 --tag phase_caller_chr3

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

STEP3_PATH=""
DOCKER=""

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# Command line argument parsing from :
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
    --step3-path)
    STEP3_PATH="$2"
    shift # past argument
    shift # past value
    ;;
    -d|--destination)
    DESTINATION="$2"
    shift
    shift
    ;;
    --docker)
    DOCKER="$2"
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

if [ -z "${STEP3_PATH}" ]
then
    echo "Provide the path of step3 for the chromosome --step3-path <path>"
    exit 1
fi

if [ -z "${DOCKER}" ]
then
    echo "Please provide the docker image ID or path with --docker <value>"
    exit 1
fi

if [ -z "${CHROMOSOME}" ]
then
    echo "Specify the chromosome with --chromosome <value>"
    exit 1
fi

if [ -z "${TAG}" ]
then
    tag="step4_merge_regions"
else
    tag="${TAG}"
fi
echo "dx run with tag : ${tag}"

echo "Instance type : ${INSTANCE}"

LAUNCH_ALL=""

echo "Will launch merge for files :"
echo $(dx ls ${STEP3_PATH})
echo "With docker image: ${DOCKER} - $(path_to_dx_id "${DOCKER}")"
DESTINATION="${DESTINATION}/SAPPHIRE_step4/${CHROMOSOME}"
echo "Results will be stored in ${DESTINATION}"

if [ -z "${LAUNCH_ALL}" ]
then
    ask_permission_to_launch_all
fi
if [ "${SKIP_JOB}" = "yes" ]
then
    unset SKIP_JOB
else
    dx run swiss-army-knife -icmd="/usr/src/pp/Docker/update_pp.sh; /usr/src/pp/bin_tools/scripts/region_merge_helper_script.sh --path /mnt/project/${STEP3_PATH}" \
        ${COST_LIMIT_ARG} --name "Step4: Merge regions ${CHROMOSOME}" \
        -iimage_file="${DOCKER}" --tag "${tag}" \
        --destination "${DESTINATION}" --priority normal \
        --instance-type ${INSTANCE} -y
fi
