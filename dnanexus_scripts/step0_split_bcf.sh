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

VCF_VAR_ID=""
START_POS=0
STOP_POS=""
CHUNK_SIZE=1000000
OVERLAP=2000
CHROMOSOME=""

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    --vcf-id)
    VCF_VAR_ID="$2"
    shift # past argument
    shift # past value
    ;;
    --vcf-idx-id)
    VCF_IDX_ID="$2"
    shift
    shift
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
    --chunk-size)
    CHUNK_SIZE="$2"
    shift
    shift
    ;;
    --overlap)
    OVERLAP="$2"
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

if [ -z "${VCF_VAR_ID}" ]
then
    echo "Specify an input VCF/BCF file ID with --vcf-id <actual_input_id>"
    exit 1
fi

if [ -z "${VCF_IDX_ID}" ]
then
    echo "Specify the index ID for the input VCF/BCF file with --vcf-idx-id <actual_input_id>"
    exit 1
fi

if [ -z "${CHROMOSOME}" ]
then
    echo "Specify a chromosome with --chromosome <value>"
    exit 1
fi

if [ -z "${STOP_POS}" ]
then
    echo "Specify a chunking stop position with --stop-pos <value>"
    exit 1
fi

VCF_FILENAME=$(dx_id_to_name "${VCF_VAR_ID}")
echo "FILENAME        = ${VCF_FILENAME}"


if [ -z "${TAG}" ]
then
    tag="step0_split_bcf"
else
    tag="${TAG}"
fi
echo "dx run with tag : ${tag}"

DESTINATION="${DESTINATION}/SAPPHIRE_step0/${CHROMOSOME}"

echo "Output file destination : ${DESTINATION}"
echo "Instance type : ${INSTANCE}"

echo "This will launch multiple jobs !"

LAUNCH_ALL=""
NUM_JOBS=$((STOP_POS / CHUNK_SIZE + 1))
echo Expected number of jobs : "${NUM_JOBS}"

position=0

for (( i=0 ; position < STOP_POS; position+=CHUNK_SIZE, i++ ))
do
    region="${CHROMOSOME}:$((position))-$((position+CHUNK_SIZE+OVERLAP))"
    new_file="${VCF_FILENAME}_${CHROMOSOME}.$((position))-$((position+CHUNK_SIZE+OVERLAP)).bcf"
    echo "${region}"
    command="bcftools view --threads 2 --regions-overlap 0 -r ${region} -Ob -o ${new_file} ${VCF_FILENAME} && bcftools index --threads 2 ${new_file}"
    echo "${command}"
    if [ -z "${LAUNCH_ALL}" ]
    then
        ask_permission_to_launch_all
    fi
    if [ "${SKIP_JOB}" = "yes" ]
    then
        unset SKIP_JOB
    else
        dx run swiss-army-knife -icmd="${command}" \
            -iin="${VCF_VAR_ID}" -iin="${VCF_IDX_ID}" \
            -imount_inputs=true \
            ${COST_LIMIT_ARG} --name "Step0: Split BCF ${CHROMOSOME}_${position}" \
            --tag "${tag}" \
            --destination "${DESTINATION}" --priority normal \
            --instance-type ${INSTANCE} -y
    fi
done
