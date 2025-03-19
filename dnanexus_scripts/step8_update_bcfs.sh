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
STEP1_PATH=""
STEP2_VARS=""

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# Command line argument parsing from :
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
    --step1-path)
    STEP1_PATH="$2"
    shift
    shift
    ;;
    --step2-vars)
    STEP2_VARS="$2"
    shift # past argument
    shift # past value
    ;;
    --step7-binary)
    BINARY="$2"
    shift
    shift
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

if [ -z "${STEP1_PATH}" ]
then
    echo "Please specify step1 path at the chromosome level with --step1-path <path>"
    exit 1
fi

if [ -z "${STEP2_VARS}" ]
then
    echo "Please specify step2 variant file with --step2-vars <path to file>"
    exit 1
else
    dx describe "${STEP2_VARS}" || { echo "Failed to find file ${STEP2_VARS}"; exit 1; }
fi

if [ -z "${BINARY}" ]
then
    echo "Please specify step7 merged bianry file with --step7-binary <path to file>"
    exit 1
else
    dx describe "${BINARY}" || { echo "Failed to find file ${BINARY}"; exit 1; }
fi

if [ -z "${APPLET_ID}" ]
then
    echo "Please provide the ID or path of pp-update-applet with --applet <value>"
    exit 1
else
    dx describe "${APPLET_ID}" || { echo "Failed to find file ${APPLET_ID}"; exit 1; }
fi

if [ -z "${CHROMOSOME}" ]
then
    echo "Please specify the chromosome with --chromosome <value>"
    exit 1
fi

if [ -z "${TAG}" ]
then
    tag="step8_update"
else
    tag="${TAG}"
fi
echo "dx run with tag : ${tag}"

echo "Instance type : ${INSTANCE}"

LAUNCH_ALL=""

for filename in $(dx ls "${STEP1_PATH}")
do
    # Skip if file is not bcf
    if [[ $filename != *.bcf ]]
    then
        echo "File ${filename} does not end in .bcf, ignoring..."
        continue
    fi

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

    file_with_path="${STEP1_PATH}"/"${filename}"
    array+=("-isplitted_binary_files=${file_with_path}")

    inner_dest="${DESTINATION}/SAPPHIRE_step8/${CHROMOSOME}"

    echo
    echo "File to update : ${file_with_path}"
    echo "Destination: ${inner_dest}"

    if [ -z "${LAUNCH_ALL}" ]
    then
        ask_permission_to_launch_all
    fi
    if [ "${SKIP_JOB}" = "yes" ]
    then
        unset SKIP_JOB
    else
        dx run "${APPLET_ID}" \
            -ioriginal_vcf_file="${file_with_path}" \
            -imain_variant_vcf_file="${STEP2_VARS}" \
            -irephased_binary_file="${BINARY}" \
            -iverbose=true \
            ${COST_LIMIT_ARG} --name "Step8: Update BCF ${CHROMOSOME}_${start}" \
            --tag "${tag}" \
            --destination "${inner_dest}" --priority normal \
            --instance-type "${INSTANCE}" -y
    fi
done


