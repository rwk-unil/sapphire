#!/bin/bash
# bcftools-update-ac-an 0.0.1
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

    echo "Value of bcf_file: '$bcf_file'"

    FILENAME="$(dx describe "$bcf_file" --name)"
    echo "Filename : ${FILENAME}"
    NEW="${FILENAME}"_ac_an_updated.bcf
    echo "New filename : ${NEW}"

    dx download "$bcf_file" -o "${FILENAME}"

    date
    echo "Updating..."
    echo "Please ignore the warning given by bcftools on the empty sample name"
    echo "This is because asking to remove no sample is a fast way to update AC/AN"
    bcftools view -s "^" --force-samples "${FILENAME}" -Ob -o "${NEW}"
    date
    echo "Done updating !"
    echo "Indexing..."
    bcftools index "${NEW}"
    date
    echo "Done indexing !"
    NEW_INDEX=$(ls ${NEW}.*)

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

    updated_bcf_file=$(dx upload "${NEW}" --brief)
    update_bcf_file_index=$(dx upload "${NEW_INDEX}" --brief)

    # The following line(s) use the utility dx-jobutil-add-output to format and
    # add output variables to your job's output as appropriate for the output
    # class.  Run "dx-jobutil-add-output -h" for more information on what it
    # does.

    dx-jobutil-add-output updated_bcf_file "$updated_bcf_file" --class=file
    dx-jobutil-add-output update_bcf_file_index "$update_bcf_file_index" --class=file
}
