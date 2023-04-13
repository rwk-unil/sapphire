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
REFERENCE=""
unset -v FIFO_SIZE

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
    -r|--reference)
    REFERENCE="$2"
    shift # past argument
    shift # past value
    ;;
    --fifo-size)
    FIFO_SIZE="$2"
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

if [ -z "${FILENAME}" ]
then
    echo "Specify a filename with --filename, -f <filename>"
    exit 1
fi

if [ -z "${REFERENCE}" ]
then
    echo "Specify a filename with --reference, -r <filename>"
    exit 1
fi

if ! [ -z "${FIFO_SIZE}" ]
then
    FIFO_ARG=--fifo-size
fi

echo "FILENAME        = ${FILENAME}"
echo "REFERENCE       = ${REFERENCE}"

TMPDIR=$(mktemp -d -t pp_XXXXXX) || { echo "Failed to create temporary directory"; exit 1; }

echo "Temporary directory : ${TMPDIR}"

function exit_fail_rm_tmp {
    echo "Removing directory : ${TMPDIR}"
    rm -r ${TMPDIR}
    exit 1
}

OUTPUTNAME="$(basename "${FILENAME}")".bin

"${SCRIPTPATH}"/../../pp_extractor/pp_extract ${FIFO_ARG} ${FIFO_SIZE} -f "${FILENAME}" -o ${TMPDIR}/"${OUTPUTNAME}" || { echo "Failed to extract ${FILENAME}"; exit_fail_rm_tmp; }
cmp "${REFERENCE}" ${TMPDIR}/"${OUTPUTNAME}" || { echo "[KO] Output file and reference are different"; exit_fail_rm_tmp; }

echo "[OK] The extracted file and reference are the same"

rm -r $TMPDIR
exit 0
