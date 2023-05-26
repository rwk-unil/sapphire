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

VCF_FILE_ID=""
MAX_REGION_LENGTH=250000000
CHUNK_SIZE=1000000
REGION_OFFSET=1000
DRY=false

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# Command line argument parsing from :
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
    -f|--vcf-file-id)
    VCF_FILE_ID="$2"
    shift # past argument
    shift # past value
    ;;
    -i|--vcf-index-id)
    VCF_INDEX_ID="$2"
    shift
    shift
    ;;
    -m|--max-region-length)
    MAX_REGION_LENGTH="$2"
    shift
    shift
    ;;
    -c|--chr)
    CHROMOSOME="$2"
    shift
    shift
    ;;
    -d|--destination)
    DESTINATION="$2"
    shift
    shift
    ;;
    --dry-run)
    DRY=true
    shift
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

if [ -z "${VCF_FILE_ID}" ]
then
    echo "Requires VCF/BCF file ID, -f,--vcf-file-id"
    exit 1
fi

if [ -z "${VCF_INDEX_ID}" ]
then
    echo "Requires VCF/BCF index file ID, -i,--vcf-index-id"
    exit 1
fi

if [ -z "${CHROMOSOME}" ]
then
    echo "Requires chromosome -c,--chr"
    exit 1
fi

dx select

REGION=0

vcf_filename=$(dx describe ${VCF_FILE_ID} --json | jq -r .name)

while [ ${REGION} -lt ${MAX_REGION_LENGTH} ]
do
    bcf_region="${CHROMOSOME}:${REGION}-$((${REGION}+${CHUNK_SIZE}))"
    new_name="${vcf_filename}-${bcf_region}.bcf"
    command="date; time bcftools view --threads 2 -Ob -o ${new_name} -r ${bcf_region} --regions-overlap pos ${vcf_filename} && date; time bcftools index ${new_name}"

if [ ${DRY} == true ]
then
    echo "${command}\t Destination=${DESTINATION}"
else
    dx run swiss-army-knife -icmd="${command}" \
        -imount_inputs=true \
        -iin=${VCF_FILE_ID} -iin=${VCF_INDEX_ID} \
        --name "extract - ${bcf_region}" \
        --tag "extract_${CHROMOSOME}" \
        --destination "${DESTINATION}" \
        --instance-type ${INSTANCE} -y
fi

    # Subtract for overlap
    REGION=$((${REGION}+${CHUNK_SIZE}-${REGION_OFFSET}))
done

exit 0