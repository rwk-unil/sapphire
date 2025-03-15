#!/bin/bash

echo "$(date): Removing GT..."

mkdir -p work

for file in $(ls *.bcf)
do
    echo "$(date): File to remove GT : ${file}"
    filename=$(basename ${file})

    # Extract the prefix (everything before the last underscore)
    prefix="${filename%_*}"
    # Extract the part after the last underscore
    suffix="${filename##*_}"

    # Extract chr, start, and stop using pattern matching
    chr="${suffix%%.*}"   # Everything before the first dot

    bcftools view -G -Ob -o "work/${filename}" ${filename} && bcftools index "work/${filename}"
done

# Sort for the naive concat
cd work
array=($(ls *.bcf | awk -F'[_.-]' '{print $(NF-2), $0}' | sort -n | cut -d' ' -f2-))
echo "${array[@]}"

echo "$(date): Concatenating..."

new_file="${prefix}_${chr}.bcf"
bcftools concat --naive -o "../${new_file}" -Ob ${array[@]}
cd ..
rm -rf work