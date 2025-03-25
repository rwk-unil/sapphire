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

CRAM_PATH=""
DESTINATION="/"

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# Command line argument parsing from :
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
    --cram-path)
    CRAM_PATH="$2"
    shift # past argument
    shift # past value
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

if [ -z "${CRAM_PATH}" ]
then
    echo "Provide the path of where the CRAM files are on DNANexus --cram-path <path>"
    echo "This path should have the sub folders starting with two digits, and withing the CRAMs"
    echo "E.g., \"Bulk/DRAGEN WGS/Whole genome CRAM files (DRAGEN) [500k release]\""
    exit 1
fi

TMPDIR=$(mktemp -d -t pp_XXXXXX) || { echo "Failed to create temporary directory"; exit 1; }

echo "Temporary directory : ${TMPDIR}"

function exit_fail_rm_tmp {
    echo "Removing directory : ${TMPDIR}"
    rm -r ${TMPDIR}
    exit 1
}

echo "This will take a few minutes"

for directory in $(dx ls "${CRAM_PATH}")
do
    echo "Querying ${CRAM_PATH}/${directory}"
    dx ls -l "${CRAM_PATH}"/${directory} | tail -n+4 >> "${TMPDIR}"/cram_paths.txt
done

cat "${TMPDIR}"/cram_paths.txt | \
 awk '{size=$(NF-3) " " $(NF-2); split($(NF-1), a, "_"); print a[1] "," $(NF-1) "," $NF "," size }' | \
 sed 's/[()]//g' | \
 awk 'NR%2==1 {first=$0; next} {print first "," $0}' > \
 "${TMPDIR}"/cram_paths.csv

mkdir -p files
cp "${TMPDIR}"/cram_paths.csv files/cram_paths.csv

dx upload "${TMPDIR}"/cram_paths.csv --destination "${DESTINATION}/cram_paths.csv"

# awk -F',' 'NR==FNR {map[$1] = $0; next} $2 in map {print $0 "," map[$2]} !($2 in map) {print $0 ",,,,,,,,,,"}' cram_paths.csv sample_list.csv > test.txt
# awk -F',' '{OFS=","; $3=$5; print $0}' test.txt | less


# Merge sample list file (with just samples)
# awk -F',' 'NR==FNR {map[$1] = $0; next} $1 in map {print map[$1]} !($1 in map) {print $0 ",,,,,,,,,,"}' cram_paths.csv sample_list.txt | awk '{print NR-1 "," $0 }' > cram_paths_for_samples.csv
# cut -d, -f1-3 cram_paths_for_samples.csv > sample_list.csv

# Print lines between N and M (1 based)
# sed -n 'N,Mp' sample_list.csv
# sed -n '200000,201000p' sample_list.csv
# sed -n '1,1000p' sample_list.csv

# Generate the array of CRAM files
# sed -n '1,1000p' cram_paths_for_samples.csv | cut -d, -f 4 | sed 's/^/-iin=/'

# Generate the array of CRAM index files
# sed -n '1,1000p' cram_paths_for_samples.csv | cut -d, -f 8 | sed 's/^/-iin=/'

rm -r ${TMPDIR}