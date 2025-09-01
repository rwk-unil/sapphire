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

INSTANCE="mem2_ssd1_v2_x4"
PRIORITY="high"

# Source common variables and functions
source "${SCRIPTPATH}/common.sh"

BINARY=""
XCF_BCF=""

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# Command line argument parsing from :
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
    --xcf)
    XCF_BCF="$2"
    shift
    shift
    ;;
    --bin)
    BINARY="$2"
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

if [ -z "${XCF_BCF}" ]
then
    echo "Please specify path for the original XCF file (requires .csi, .bin. fam in same directory)"
    exit 1
fi

if [ -z "${BINARY}" ]
then
    echo "Please specify step7 merged bianry file with --step7-binary <path to file>"
    exit 1
else
    dx describe "${BINARY}" || { echo "Failed to find file ${BINARY}"; exit 1; }
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

for filename in $(dx ls "${XCF_BCF}")
do
    # Skip if file is not bcf
    if [[ $filename != *.bcf ]]
    then
        echo "File ${filename} does not end in .bcf, ignoring..."
        continue
    fi

    # Get the path of the XCF BCF file and generate paths for the other files
    xcf_path="$(dx_id_to_dx_path "${XCF_BCF}")"
    base="${xcf_path}/$(dx_id_to_name ${XCF_BCF})"
    base="${base%.*}"
    XCF_BCF_FILENAME="$(basename ${filename})"
    XCF_BIN="${base}.bin"
    XCF_IDX="${base}.bcf.csi"
    XCF_FAM="${base}.fam"
    XCF_BIN_FILENAME="$(basename ${XCF_BIN})"
    XCF_IDX_FILENAME="$(basename ${XCF_IDX})"
    XCF_FAM_FILENAME="$(basename ${XCF_FAM})"
    SAPPHIRE_BIN_FILENAME="$(dx_id_to_name ${BINARY})"
    NEW_XCF_FILENAME="${XCF_BCF_FILENAME%.*}_rephased.bcf"

    inner_dest="${DESTINATION}/SAPPHIRE_step8/${CHROMOSOME}"

    echo
    echo "File to update : ${XCF_BCF_FILENAME}"
    echo "SAPPHIRE binary file : ${SAPPHIRE_BIN_FILENAME}"
    echo "New XCF filename : ${NEW_XCF_FILENAME}"
    echo "Destination: ${inner_dest}"

    command="ls; dx download INTERVENE_PHASING:/csgdata/rick/xcftools_4a87d20 && mv xcftools_4a87d20 xcftools && chmod +x xcftools && mv xcftools /usr/bin/ && \
            xcftools view --threads \$(nproc) -i ${XCF_BCF_FILENAME} --sapphire-bin ${SAPPHIRE_BIN_FILENAME} -Osb -o ${NEW_XCF_FILENAME}"

    echo "Command : ${command}"

    if [ -z "${LAUNCH_ALL}" ]
    then
        ask_permission_to_launch_all
    fi
    if [ "${SKIP_JOB}" = "yes" ]
    then
        unset SKIP_JOB
    else
        dx run "swiss-army-knife" \
            -iin="${XCF_BCF}" -iin="${XCF_IDX}" -iin="${XCF_BIN}" -iin="${XCF_FAM}" \
            -iin="${BINARY}" -imount_inputs=true \
            -icmd="${command}" \
            ${COST_LIMIT_ARG} --name "Step8: Update XCF ${CHROMOSOME}" \
            --tag "${tag}" \
            --destination "${inner_dest}" --priority ${PRIORITY} \
            --instance-type "${INSTANCE}" -y
    fi
done


