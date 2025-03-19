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

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# Command line argument parsing from :
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
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

command="wget https://github.com/rwk-unil/sapphire/releases/download/v1.0.0-paper/pp_toolkit_v1.4.tar.gz.part_1 && \
         wget https://github.com/rwk-unil/sapphire/releases/download/v1.0.0-paper/pp_toolkit_v1.4.tar.gz.part_2 && \
         cat pp_toolkit_v1.4.tar.gz.part_1 pp_toolkit_v1.4.tar.gz.part_2 > pp_toolkit_v1.4.tar.gz && \
         rm pp_toolkit_v1.4.tar.gz.part_1 && rm pp_toolkit_v1.4.tar.gz.part_2"

echo "Command : ${command}"
if [ -z "${LAUNCH_ALL}" ]
then
ask_permission_to_launch_all
fi
if [ "${SKIP_JOB}" = "yes" ]
then
unset SKIP_JOB
else
dx run swiss-army-knife -icmd="${command}" \
        ${COST_LIMIT_ARG} --name "Upload SAPPHIRE docker image" \
        --destination "${DESTINATION}" --priority normal \
        --instance-type ${INSTANCE} -y
fi