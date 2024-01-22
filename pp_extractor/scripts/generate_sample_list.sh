#!/bin/bash

FILENAME=""

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# Command line argument parsing from :
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
    -f|--filename)
    FILENAME="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

if [ -z "${FILENAME}" ]
then
    echo "Specify a filename with --filename, -f <filename>"
    exit 1
fi

echo "FILENAME        = ${FILENAME}"

# Extract sample list (fast so not launched in background)
#bcftools query --list-samples "${FILENAME}" | nl -s ',' | awk '{$1=$1};1' > "$(basename ${FILENAME}).samples.csv"
bcftools query --list-samples "${FILENAME}" | nl -s ',' | awk '{$1=$1};1' > "${FILENAME}.samples.csv"
