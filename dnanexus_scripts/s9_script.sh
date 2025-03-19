#!/bin/bash

echo $1
echo "$(date): Concatenating BCFs..."

# Sort (on genomic region in filename) for the naive concat
# This sorting here is very hacky, it uses the "_.-" delimitors to extract the chromosome locations
# From the filename by discarding the 3 last elements (.bcf_rephased.bcf)
# Before that will be chrN.start-stop and here we sort numerically by stop position
# This will break very easily if the filenames are different, so be aware...
array=($(ls *.bcf | awk -F'[_.-]' '{print $(NF-3), $0}' | sort -n | cut -d' ' -f2-))
echo "${array[@]}"

if [ -z "$1" ]
then
    # Dump the genomics region and rephased extension
    prefix="${array[1]%_*.*-*.bcf_rephased.bcf}"
else
    prefix="$1"
fi
# Use the prefix for the new file name
new_file="${prefix}.bcf"

echo "$(date): Concatenating into ${new_file} ..."

bcftools concat --naive -o "${new_file}" -Ob ${array[@]}
bcftools index --threads $(nproc) "${new_file}"
