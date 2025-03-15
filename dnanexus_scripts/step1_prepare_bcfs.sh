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

OVERLAP=2000

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# Command line argument parsing from :
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
    --step0-path)
    STEP0_PATH="$2"
    shift # past argument
    shift # past value
    ;;
    -d|--destination)
    DESTINATION="$2"
    shift
    shift
    ;;
    --stop-pos)
    STOP_POS="$2"
    shift
    shift
    ;;
    --overlap)
    OVERLAP="$2"
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

if [ -z "${STEP0_PATH}" ]
then
    echo "Specify the step0 path with --step0-path <actual_path_at_chromosome_level>"
    exit 1
fi

echo "PATH            = ${STEP0_PATH}"

if [ -z "${TAG}" ]
then
    tag="step1_prepare_bcf"
else
    tag="${TAG}"
fi
echo "dx run with tag : ${tag}"

echo "Instance type : ${INSTANCE}"

LAUNCH_ALL=""

position=0

for filename in $(dx ls ${STEP0_PATH})
do
    # Skip if file ends with .csi
    if [[ $filename == *.csi ]]
    then
        continue
    fi
    echo "File to prepare : ${filename}"
    file_with_path="${STEP0_PATH}/${filename}"
    echo ${file_with_path}
    file_id=$(path_to_dx_id ${file_with_path})
    idx_with_path="${file_with_path}.csi"
    echo ${idx_with_path}
    idx_id=$(path_to_dx_id ${idx_with_path})
    # Extract the prefix (everything before the last underscore)
    prefix="${filename%_*}"
    # Extract the part after the last underscore
    suffix="${filename##*_}"

    # Extract chr, start, and stop using pattern matching
    chr="${suffix%%.*}"   # Everything before the first dot
    rest="${suffix#*.}"   # Everything after the first dot
    start="${rest%-*}"    # Everything before the last '-'
    stop="${rest#*-}"     # Everything after the last '-'
    stop="${stop%.bcf}"   # Remove ".bcf"
    stop=$(((stop - OVERLAP) - 1))

    inner_dest="${DESTINATION}/SAPPHIRE_step1/${chr}"
    echo "Output file destination : ${inner_dest}"

    region="${chr}:$((start))-$((stop))"
    new_file="${prefix}_${chr}.$((start))-$((stop)).bcf"
    echo "${region}"
    command="bcftools view --threads 2 --regions-overlap 0 -r ${region} -Ob -o ${new_file} ${filename} && bcftools index --threads 2 ${new_file}"
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
            -iin=${file_id} -iin=${idx_id} \
            -imount_inputs=true \
            ${COST_LIMIT_ARG} --name "Step1: Prepare BCF ${chr}_${start}" \
            --tag "${tag}" \
            --destination "${inner_dest}" --priority normal \
            --instance-type ${INSTANCE} -y
    fi
done
