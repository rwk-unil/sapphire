#!/bin/bash

if ! command -v realpath &> /dev/null
then
    realpath() {
        [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
    }
fi

# Get the path of this script
SCRIPTPATH=$(realpath  $(dirname "$0"))

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
bcftools query --list-samples "${FILENAME}" > "${FILENAME}.samples.txt"
# Extract variants (can take a long time)
bcftools view -G -Ob "${FILENAME}" -o "${FILENAME}_vars.bcf" &
# Extract the low PPs
pp_extract -f "${FILENAME}" -o "${FILENAME}_hets.bin" &
# Wait for the subprocesses to finish
wait