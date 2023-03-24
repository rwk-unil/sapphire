# PP Extractor tool

This is a tool meant to be run on BCF files phased with SHAPEIT5 https://github.com/odelaneau/shapeit5 it will extract heterozygous variants with their `PP` (phasing probability) field below 0.99 alongside heterozygous variants that come before and after. The extracted heterozygous variants are placed in a sparse binary file (see doc/Binary_Format.md).

This file will allow extremely fast access to the variants and is used as input to rephase them using sequencing data (BAM/CRAM).