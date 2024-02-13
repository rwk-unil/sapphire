#!/bin/bash

# Depends on bcftools
# Depends on XCFTools
# Depends on SAPPHIRE

threads=1
overlap=2000
maf=0.00001
directory="."

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# Command line argument parsing from :
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
    -f|--vcf-filename)
    VCF_FILENAME="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--output)
    OFNAME="$2"
    shift # past argument
    shift # past value
    ;;
    -d|--dir)
    directory="$2"
    shift
    shift
    ;;
    --threads)
    threads="$2"
    shift
    shift
    ;;
    --overlap)
    overlap="$2"
    shift
    shift
    ;;
    -m|--maf)
    maf="$2"
    shift
    shift
    ;;
    --verbose)
    VERBOSE="-v"
    shift # no value attached
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

set -e
set -o pipefail

start_position=0
position=0
end_position=$(bcftools view --threads ${threads} -H -G ${VCF_FILENAME} | tail -n 1 | cut -f 2) || \
    { echo failed to get end position; exit 1; }
chromosome=$(bcftools index -s ${VCF_FILENAME} | tail -n 1 | cut -f 1) || \
    { echo failed to get chromosome; exit 1; }
vcf_basename=$(basename "${VCF_FILENAME}")

if [ -z "$OFNAME" ]
then
    OFNAME="${vcf_basename}_sapphire.bin"
fi

if ((overlap >= 10000))
then
    echo "Overlap too big (should be <10000)"
    exit 1
fi

mkdir -p "${directory}" || \
    { echo "failed to create directory : ${directory}"; exit 1; }

start=$(date)
echo Starting processing : ${start}

if (( end_position < (10000 * threads) || threads == 1 ))
then
    # Don't bother with multi-threading
    echo "single thread"
    xcftools view -i "${VCF_FILENAME}" -m "${maf}" -o "${directory}/${OFNAME}" --log "${directory}/${OFNAME}.log" -Obs
else
    echo "Running with ${threads} threads"
    fail=0
    increment=$(((end_position+overlap)/threads))
    #echo $threads
    #echo $increment
    #echo $end_position
    for (( i=0 ; position < end_position; position+=increment, i++))
    do
        echo "${chromosome}:$((position))-$((position+increment+overlap))"
        xcftools view --line-from-vcf -r "${chromosome}:$((position))-$((position+increment+overlap))" -i "${VCF_FILENAME}" -m "${maf}" -o "${directory}/${OFNAME}_$i" --log "${directory}/${OFNAME}_$i.log" -Obs > /dev/null 2>&1 &
        pids[$i]=$!
    done

    for pid in ${pids[*]}
    do
        wait $pid || let "fail+=1"
    done

    if (( fail > 0 ))
    then
        echo ${fail} jobs failed
    else
        vertical_bin_merger -b "${directory}/${OFNAME}" -o "${directory}/${OFNAME}"
    fi
fi

stop=$(date)
echo Done processing : ${stop}