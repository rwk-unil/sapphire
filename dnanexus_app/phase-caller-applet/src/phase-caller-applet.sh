#!/bin/bash
# phase-caller-applet 0.0.1
# Generated by dx-app-wizard.
#
# Basic execution pattern: Your app will run on a single machine from
# beginning to end.
#
# Your job's input variables (if any) will be loaded as environment
# variables before this script runs.  Any array inputs will be loaded
# as bash arrays.
#
# Any code outside of main() (or any entry point you may add) is
# ALWAYS executed, followed by running the entry point itself.
#
# See https://documentation.dnanexus.com/developer for tutorials on how
# to modify this file.

main() {

    # dxfuse
    # Create a manifest file for dxfuse
    echo "{
    \"files\" : [],
    \"directories\" : [
        {
        \"proj_id\" : \"$DX_PROJECT_CONTEXT_ID\",
        \"folder\" : \"/\",
        \"dirname\" : \"/project\"
        }
    ]
    }" > .dxfuse_manifest.json
    echo "DX_PROJECT_CONTEXT_ID : ${DX_PROJECT_CONTEXT_ID}"

    # Create a mount point for the project
    MOUNTDIR=/mnt
    sudo mkdir -p $MOUNTDIR

    # Mount the current project
    dxfuse $MOUNTDIR .dxfuse_manifest.json &

    # Wait for mount to start
    sleep 2

    echo "Value of binary_file: '$binary_file'"
    echo "Value of vcf_file_without_samples: '$vcf_file_without_samples'"
    echo "Value of sample_list: '$sample_list'"
    echo "Value of suffix: '$suffix'"
    echo "Value of cram_path: '$cram_path'"
    echo "Value of verbose: '$verbose'"
    echo "Value of threads: '$threads'"

    # The following line(s) use the dx command-line tool to download your file
    # inputs to the local file system using variable names for the filenames. To
    # recover the original filenames, you can use the output of "dx describe
    # "$variable" --name".

    BINARY_FILENAME="$(dx describe "$binary_file" --name)"
    dx download "$binary_file" -o "${BINARY_FILENAME}"
    VCF_FILENAME="$(dx describe "$vcf_file_without_samples" --name)"
    dx download "$vcf_file_without_samples" -o "${VCF_FILENAME}"
    SAMPLE_FILENAME="$(dx describe "$sample_list" --name)"
    dx download "$sample_list" -o "${SAMPLE_FILENAME}"

    if [ "$verbose" = true ]
    then
        VERBOSE="-v"
    else
        unset VERBOSE
    fi

    if ! [ -z "${suffix}" ]
    then
        NEW_BINARY_FILE="${BINARY_FILENAME}_${suffix}.bin"
    else
        NEW_BINARY_FILE="${BINARY_FILENAME}_rephased.bin"
    fi

    # Hacky way to get project ID number ...
    PROJECT_ID="$(ls "/mnt/project/${cram_path}/20/" 2> /dev/null | grep ".cram" 2> /dev/null | head -n1 | sed -n 's/[0-9]*_\([0-9]*\)_.*$/\1/p' 2> /dev/null)"
    if [ -z "${PROJECT_ID}" ]
    then
        dx-jobutil-report-error "Could not extract project ID from cram path ${cram_path}"
        exit 1
    else
        echo "PROJECT_ID : ${PROJECT_ID}"
    fi

    CRAM_PATH="/mnt/project/${cram_path}"

    # Copy the Reference files
    WORKDIR=$(pwd)
    mkdir -p ~/.cache
    cd ~/.cache
    curl -L https://github.com/rwk-unil/cram_accessor/releases/download/v1.0/hts-ref.zip -o ./hts-ref.zip > /dev/null 2>&1
    unzip hts-ref.zip > /dev/null 2>&1
    rm -f hts-ref.zip
    cd "${WORKDIR}"

    # Copy the binary file as a new file (will be modyfied in-place by the caller)
    cp "${BINARY_FILENAME}" "${NEW_BINARY_FILE}"

    echo "time phase_caller -f \"${VCF_FILENAME}\" -b \"${NEW_BINARY_FILE}\" -S \"${SAMPLE_FILENAME}\" ${VERBOSE} -I ${PROJECT_ID} -t ${threads} --cram-path \"${CRAM_PATH}\""
    time phase_caller -f "${VCF_FILENAME}" -b "${NEW_BINARY_FILE}" -S "${SAMPLE_FILENAME}" ${VERBOSE} -I ${PROJECT_ID} -t ${threads} --cram-path "${CRAM_PATH}"

    # Fill in your application code here.
    #
    # To report any recognized errors in the correct format in
    # $HOME/job_error.json and exit this script, you can use the
    # dx-jobutil-report-error utility as follows:
    #
    #   dx-jobutil-report-error "My error message"
    #
    # Note however that this entire bash script is executed with -e
    # when running in the cloud, so any line which returns a nonzero
    # exit code will prematurely exit the script; if no error was
    # reported in the job_error.json file, then the failure reason
    # will be AppInternalError with a generic error message.

    # The following line(s) use the dx command-line tool to upload your file
    # outputs after you have created them on the local file system.  It assumes
    # that you have used the output field name for the filename for each output,
    # but you can change that behavior to suit your needs.  Run "dx upload -h"
    # to see more options to set metadata.

    rephased_binary_file=$(dx upload "${NEW_BINARY_FILE}" --brief)

    # The following line(s) use the utility dx-jobutil-add-output to format and
    # add output variables to your job's output as appropriate for the output
    # class.  Run "dx-jobutil-add-output -h" for more information on what it
    # does.

    dx-jobutil-add-output "${NEW_BINARY_FILE}" "$rephased_binary_file" --class=file
}
