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

INSTANCE="mem3_ssd1_v2_x2"

# Source common variables and functions
source "${SCRIPTPATH}/common.sh"

STEP0_PATH=""
APPLET_ID=""
BCF_VAR_ID=""
MAF_THRESHOLD=""

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
    --applet)
    APPLET_ID="$2"
    shift
    shift
    ;;
    --step2-var)
    BCF_VAR_ID="$2"
    shift
    shift
    ;;
    -d|--destination)
    DESTINATION="$2"
    shift
    shift
    ;;
    --maf-threshold)
    MAF_THRESHOLD="$2"
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

if [ -z "${BCF_VAR_ID}" ]
then
    echo "Please provide the path or ID of the BCF with the variants from step2 with --step2-var <ID>"
    exit 1
fi

if [ -z "${APPLET_ID}" ]
then
    echo "Please provide the ID of pp-extract-split-applet with --applet <ID>"
    exit 1
fi

MAF_ARG=""
if [ ! -z "$MAF_THRESHOLD" ]
then
    MAF_ARG="-imaf_threshold=${MAF_THRESHOLD}"
fi

echo "PATH            = ${STEP0_PATH}"

if [ -z "${TAG}" ]
then
    tag="step3_pp_extract"
else
    tag="${TAG}"
fi
echo "dx run with tag : ${tag}"

echo "Instance type : ${INSTANCE}"

LAUNCH_ALL=""

position=0

for filename in $(dx ls "${STEP0_PATH}")
do
    # Skip if file is not BCF
    if [[ $filename != *.bcf ]]
    then
        continue
    fi
    echo "File to extract from : ${filename}"
    file_with_path="${STEP0_PATH}/${filename}"
    echo ${file_with_path}
    echo "Complete variant file : $(dx_id_to_dx_path_and_name ${BCF_VAR_ID})"

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

    file_id_to_extract="$(path_to_dx_id ${file_with_path})"
    if [[ -z "${file_id_to_extract}" ]]
    then
        exit 1
    fi

    inner_dest="${DESTINATION}/SAPPHIRE_step3/${chr}"

    echo "Destination: ${inner_dest}"

    if [ -z "${LAUNCH_ALL}" ]
    then
        ask_permission_to_launch_all
    fi
    if [ "${SKIP_JOB}" = "yes" ]
    then
        unset SKIP_JOB
    else
        dx run "${APPLET_ID}" -ivcf_bcf_file="${file_id_to_extract}" -ivars_vcf_bcf_file=${BCF_VAR_ID} ${MAF_ARG} \
            ${COST_LIMIT_ARG} --name "Step3: Extract PP ${chr}_${start}_${stop}" \
            --tag "${tag}" \
            --destination "${inner_dest}" --priority normal \
            --instance-type ${INSTANCE} -y
    fi
done
