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

INSTANCE="mem3_ssd1_v2_x4"
THREADS_ARG="-t 12"

# Source common variables and functions
source "${SCRIPTPATH}/common.sh"

KEEP_TRACK_IN_FILE=""
STEP5_PATH=""
CHROMOSOME=""
BCF_VAR_ID=""
SPLIT_SIZE=1000
SAMPLE_LIST=""
CRAM_PATH_FILE=""

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# Command line argument parsing from :
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
    --step5-path)
    STEP5_PATH="$2"
    shift # past argument
    shift # past value
    ;;
    --step2-var)
    BCF_VAR_ID="$2"
    shift
    shift
    ;;
    --sample-list)
    SAMPLE_LIST="$2"
    shift
    shift
    ;;
    --cram-path-file)
    CRAM_PATH_FILE="$2"
    shift
    shift
    ;;
    --chromosome)
    CHROMOSOME="$2"
    shift
    shift
    ;;
    --split-size)
    SPLIT_SIZE="$2"
    shift
    shift
    ;;
    --keep-track-in-file)
    KEEP_TRACK_IN_FILE="y"
    shift
    ;;
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

if [ -z "${STEP5_PATH}" ]
then
    echo "Specify the step5 path with --step5-path <actual_path_at_chromosome_level>"
    exit 1
fi

if [ -z "${BCF_VAR_ID}" ]
then
    echo "Please provide the path or ID of the BCF with the variants from step2 with --step2-var <ID>"
    exit 1
fi

if [ -z "${SAMPLE_LIST}" ]
then
    echo "Please provide sample list generated in step6a with --sample-list <sample list file>"
    exit 1
fi

if [ -z "${CRAM_PATH_FILE}" ]
then
    echo "Please provide cram path file generated with generate_cram_paths.sh with --cram-path-file <file>"
    exit 1
fi

if [ -z "${CHROMOSOME}" ]
then
    echo "Specify the chromosome with --chromosome <value>"
    exit 1
fi

echo "PATH            = ${STEP5_PATH}"

if [ -z "${TAG}" ]
then
    tag="step6_phase"
else
    tag="${TAG}"
fi
echo "dx run with tag : ${tag}"

echo "Instance type : ${INSTANCE}"

LAUNCH_ALL=""

TMPDIR=temporary_files
mkdir -p "${TMPDIR}"
TACKING_FILE="${TMPDIR}/${CHROMOSOME}_tracking_file.txt"

TEST=$(ls "${TMPDIR}"/cram_paths_for_samples.csv)
if [ -z "${TEST}" ]
then
    cd "${TMPDIR}"
    dx download "${CRAM_PATH_FILE}" || { echo "Failed to download ${CRAM_PATH_FILE}"; exit 1; }
    cd ..
fi

VCF_VAR_ID=$(path_to_dx_id "${BCF_VAR_ID}") || { echo "Failed to get variant file ID"; exit 1; }
SAMPLE_FILE_ID=$(path_to_dx_id "${SAMPLE_LIST}") || { echo "Failed to get sample file ID"; exit 1; }
PROJECT_ID=$(dx describe "${BCF_VAR_ID}" --json | jq -r ".project") || { echo "Failed to get project ID"; exit 1; }
DOCKER_ID=$(path_to_dx_id "${DOCKER_IMAGE}") || { echo "Failed to get Docker ID"; exit 1; }

echo "Variant file ID: ${VCF_VAR_ID}"
echo "Sample file ID: ${SAMPLE_FILE_ID}"
echo "Project ID: ${PROJECT_ID}"
echo "Docker ID: ${DOCKER_ID}"

if [ -z "${VCF_VAR_ID}" ]
then
    echo "Failed to get variant file ID"; exit 1
fi
if [ -z "${SAMPLE_FILE_ID}" ]
then
    echo "Failed to get sample file ID"; exit 1
fi
if [ -z "${PROJECT_ID}" ]
then
    echo "Failed to get project ID"; exit 1
fi
if [ -z "${DOCKER_ID}" ]
then
    echo "Failed to get Docker ID"; exit 1
fi

pids=()
for filename in $(dx ls "${STEP5_PATH}")
do
    # Skip if file is not binary
    if [[ $filename != *.bin_sub_* ]]
    then
        echo "File ${filename} does not end in .bin_sub_*, ignoring..."
        continue
    fi

    if [[ "${KEEP_TRACK_IN_FILE}" == "y" ]]
    then
        if grep -qF "${filename}" "${TRACKING_FILE}"
        then
            echo "${filename} already in tracking file, skipping..."
            continue
        fi
    fi

    echo "File to rephase : ${filename}"
    file_with_path="${STEP5_PATH}/${filename}"
    echo ${file_with_path}
    echo "Complete variant file : $(dx_id_to_dx_path_and_name ${BCF_VAR_ID})"

    # Extract the prefix (everything before _bin_sub_)
    prefix="${filename%.bin_sub_*}"
    # Extract the idx
    idx="${filename##*.bin_sub_}"

    #echo "File index is ${idx}"

    start=$((idx * SPLIT_SIZE + 1))
    stop=$(((idx + 1) * SPLIT_SIZE))

    #echo "Start ${start} Stop ${stop}"

    cram_array=($(sed -n "${start},${stop}p" "${TMPDIR}"/cram_paths_for_samples.csv | cut -d, -f 4))
    cram_index_array=($(sed -n "${start},${stop}p" "${TMPDIR}"/cram_paths_for_samples.csv | cut -d, -f 8))
    # Concat the arrays
    idx_array=( "${cram_array[@]}" "${cram_index_array[@]}" )
    #echo "${cram_array[@]}"
    #echo "${cram_index_array[@]}"
    json_ids=$(printf '%s\n' "${idx_array[@]}" | jq -R -s -c 'split("\n")[:-1] | map({id: .})')
    #echo "${json_ids[@]}"

    LOCAL_BIN_ID=$(path_to_dx_id "${file_with_path}")
    echo "Binary file to rephase ID : ${LOCAL_BIN_ID}"

    inner_dest="${DESTINATION}/SAPPHIRE_step6/${CHROMOSOME}"
    echo "Destination: ${inner_dest}"

    INNER_NEW_BINARY_FILE="${prefix}.rephased.bin_sub_${idx}"

    # Command
    #command="cp ${filename} ${INNER_NEW_BINARY_FILE}; time phase_caller -n -f $(dx_id_to_name ${VCF_VAR_ID}) -b ${INNER_NEW_BINARY_FILE} -S $(dx_id_to_name ${SAMPLE_FILE_ID}) ${THREADS_ARG} ${VERBOSE} --cram-path-from-samples-file"
    command="echo \"\$(date): downloading HTS REF\" && curl -LJO https://github.com/rwk-unil/sapphire/releases/download/v1.0.0-paper/hts-ref.zip && \
	mkdir -p ~/.cache && mv hts-ref.zip ~/.cache/ && cd ~/.cache/ && unzip hts-ref.zip && rm -f hts-ref.zip && cd && \
    git clone https://github.com/rwk-unil/sapphire.git && cd sapphire && sed -i 's/-lzstd//g' common.mk; make && cp phase_caller/phase_caller ~/ && cd ~/out/out/ && \
    cp ${filename} ${INNER_NEW_BINARY_FILE}; time ~/phase_caller -n -f $(dx_id_to_name ${VCF_VAR_ID}) -b ${INNER_NEW_BINARY_FILE} -S $(dx_id_to_name ${SAMPLE_FILE_ID}) ${THREADS_ARG} ${VERBOSE} --cram-path-from-samples-file"
    echo "Command : ${command}"

    if [ -z "${LAUNCH_ALL}" ]
    then
        ask_permission_to_launch_all
    fi
    if [ "${SKIP_JOB}" = "yes" ]
    then
        unset SKIP_JOB
    else

        # This will generate JSON output to mount all the necessary files as not to overload the dx fuse mount point
        jq -n --arg cmd "${command}" \
            --arg projectid "${PROJECT_ID}" \
            --arg fbinid "${LOCAL_BIN_ID}" \
            --arg fbcfid "${VCF_VAR_ID}" \
            --arg fsampleid "${SAMPLE_FILE_ID}" \
            --argjson input_ids "$json_ids" '
            {
                "cmd": $cmd,
                "mount_inputs": true,
                "in": ( $input_ids | map({ "$dnanexus_link": { "project": $projectid, "id": .id } })
                    + [{ "$dnanexus_link": { "project": $projectid, "id": $fbinid } }]
                    + [{ "$dnanexus_link": { "project": $projectid, "id": $fbcfid } }]
                    + [{ "$dnanexus_link": { "project": $projectid, "id": $fsampleid } }])
            }
            ' | dx run swiss-army-knife -f - \
            ${COST_LIMIT_ARG} --name "Step6: Phase ${CHROMOSOME} batch ${idx}" \
            --tag "${tag}" \
            --destination "${inner_dest}" --priority normal \
            --instance-type ${INSTANCE} -y --brief --allow-ssh &
        pids+=($!)

        echo "Job submitted"
        echo

        if [[ "${KEEP_TRACK_IN_FILE}" == "y" ]]
        then
            if grep -qF "${filename}" "${TRACKING_FILE}"
            then
                #echo "${filename} already in tracking file, skipping..."
                continue
            else
                echo "${filename}" >> "${TRACKING_FILE}"
            fi
        fi
    fi
done

fail=0
echo "Waiting for submissions to finish"
for pid in ${pids[*]}
do
    wait $pid || let "fail+=1"
done

if (( fail > 0 ))
then
    echo "${fail} jobs failed to be submitted"
else
    echo "All jobs submitted"
fi