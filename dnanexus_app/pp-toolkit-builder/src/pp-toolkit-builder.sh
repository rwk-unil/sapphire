#!/bin/bash
# pp-toolkit-builder 0.0.1
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
    # Clone and build the tools
    git clone https://rwk-unil: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@github.com/rwk-unil/pp.git
    cd pp
    ./prepare_build.sh

    # Copy the tools
    mkdir -p dnanexus_app/pp-extract-applet/resources/usr/bin
    cp bin_tools/bin_splitter dnanexus_app/pp-extract-applet/resources/usr/bin/
    cp bin_tools/bin_merger dnanexus_app/pp-extract-applet/resources/usr/bin/
    cp pp_extractor/pp_extract dnanexus_app/pp-extract-applet/resources/usr/bin/
    cp pp_update/pp_update dnanexus_app/pp-extract-applet/resources/usr/bin/
    
    # Package the applet
    dx build dnanexus_app/pp-extract-applet
    #dx mv pp-extract-applet pp-toolkit

    # Since we cannot do dx upload --brief to get the ID (because dx build uploads)
    # Note : The dx build --brief doesn't output in the correct format !
    # We need to retrieve the ID ourselves
    pp_toolkit="$(dx describe pp-extract-applet --json | jq -r .id)"

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

    # The following line(s) use the utility dx-jobutil-add-output to format and
    # add output variables to your job's output as appropriate for the output
    # class.  Run "dx-jobutil-add-output -h" for more information on what it
    # does.

    dx-jobutil-add-output pp_toolkit "$pp_toolkit" --class=applet
}