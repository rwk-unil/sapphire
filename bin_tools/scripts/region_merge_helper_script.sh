#!/bin/bash

if ! command -v realpath &> /dev/null
then
    realpath() {
        [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
    }
fi

# Get the path of this script
SCRIPTPATH=$(realpath  $(dirname "$0"))

TMPDIR=/tmp/temporary_dir

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# Command line argument parsing from :
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
    -p|--path)
    BIN_PATH="$2"
    shift
    shift
    ;;
    -d|--directory)
    TMPDIR="$2"
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

if [ -z "${BIN_PATH}" ]
then
    echo "Please provide path to bin files with --path <path>"
    exit 1
fi

mkdir -p ${TMPDIR}
echo "Temporary directory : ${TMPDIR}"

# Sort the array based on the chromosomal start position, sep "_.-" with NF-4 because it ends with _chr.start-stop.bcf_hets.bin
array=($(ls "${BIN_PATH}"/*.bin | awk -F'[_.-]' '{print $(NF-4), $0}' | sort -n | cut -d' ' -f2-))

out_dir=$(pwd)

# Create symlinks (in sorted order) in temporary directory with naming convention for the merger
cd ${TMPDIR}
i=0
for file in ${array[@]}
do
    filename=$(basename "${file}")
    # Remove the hets.bin suffix
    filename="${filename%_hets.bin}"
    # Extract the prefix (everything before the last underscore)
    prefix="${filename%_*}"

    ln -s "${file}" "${prefix}".bin_${i}
    i=$((i+1))
done

ls -l

# Do the vertical bin merge
vertical_bin_merger -b "${prefix}".bin -o "${out_dir}"/"${prefix}".bin