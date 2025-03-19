#!/bin/bash

echo "$(date): Removing GT..."

mkdir -p work

# The number of parallel jobs, allocate as many as cpu cores
N=$(nproc)
RUNNING_JOBS=0
pids=()

for file in $(ls *.bcf)
do
    if (( RUNNING_JOBS >= N ))
    then
        wait -n
        ((RUNNING_JOBS--))
    fi
    echo "$(date): File to remove GT : ${file}"
    filename=$(basename ${file})

    # Extract the prefix (everything before the last underscore)
    prefix="${filename%_*}"
    # Extract the part after the last underscore
    suffix="${filename##*_}"

    # Extract chr, start, and stop using pattern matching
    chr="${suffix%%.*}"   # Everything before the first dot

    { bcftools view -G -Ob -o "work/${filename}" ${filename} && bcftools index "work/${filename}"; } &
    pids+=($!)
    ((RUNNING_JOBS++))
done

# Wait for all jobs to finish
echo "Waiting for threads to finish"
for pid in ${pids[*]}
do
    wait $pid || let "fail+=1"
done

if (( fail > 0 ))
then
    echo "${fail} threads failed"
else
    echo "All threads successfully completed"
fi

# Sort (on genomic region in filename) for the naive concat
cd work
array=($(ls *.bcf | awk -F'[_.-]' '{print $(NF-2), $0}' | sort -n | cut -d' ' -f2-))
echo "${array[@]}"

echo "$(date): Concatenating..."

new_file="${prefix}_${chr}.bcf"
bcftools concat --naive -o "../${new_file}" -Ob ${array[@]}
cd ..
rm -rf work
