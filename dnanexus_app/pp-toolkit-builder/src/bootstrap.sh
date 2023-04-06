#!/bin/bash

# Copy the tools
mkdir -p dnanexus_app/pp-extract-applet/resources/usr/bin
cp pp_extractor/pp_extract dnanexus_app/pp-extract-applet/resources/usr/bin/
mkdir -p dnanexus_app/bin-splitter-applet/resources/usr/bin
cp bin_tools/bin_splitter dnanexus_app/bin-splitter-applet/resources/usr/bin/
mkdir -p dnanexus_app/phase-caller-applet/resources/usr/bin
cp phase_caller/phase_caller dnanexus_app/phase-caller-applet/resources/usr/bin/
# TODO cp dxfuse !
# TODO cp htsref
# curl -L https://github.com/rwk-unil/cram_accessor/releases/download/v1.0/hts-ref.zip -o ./hts-ref.zip
mkdir -p dnanexus_app/bin-merger-applet/resources/usr/bin
cp bin_tools/bin_merger dnanexus_app/bin-merger-applet/resources/usr/bin/
mkdir -p dnanexus_app/pp-update-applet/resources/usr/bin
cp pp_update/pp_update dnanexus_app/pp-update-applet/resources/usr/bin/

# Package the applets
dx build dnanexus_app/pp-extract-applet
dx build dnanexus_app/bin-splitter-applet
dx build dnanexus_app/phase-caller-applet
dx build dnanexus_app/bin-merger-applet
dx build dnanexus_app/pp-update-applet

# Since we cannot do dx upload --brief to get the ID (because dx build uploads)
# Note : The dx build --brief doesn't output in the correct format !
# We need to retrieve the ID ourselves
pp_extract_applet="$(dx describe pp-extract-applet --json | jq -r .id)"
bin_splitter_applet="$(dx describe bin-splitter-applet --json | jq -r .id)"
phase_caller_applet="$(dx describe phase-caller-applet --json | jq -r .id)"
bin_merger_applet="$(dx describe bin-merger-applet --json | jq -r .id)"
pp_update_applet="$(dx describe pp-update-applet --json | jq -r .id)"

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

dx-jobutil-add-output pp_extract "$pp_extract_applet" --class=applet
dx-jobutil-add-output bin_splitter "$bin_splitter_applet" --class=applet
dx-jobutil-add-output phase_caller "$phase_caller_applet" --class=applet
dx-jobutil-add-output bin_merger "$bin_merger_applet" --class=applet
dx-jobutil-add-output pp_update "$pp_update_applet" --class=applet
