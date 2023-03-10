#!/bin/bash

if ! command -v realpath &> /dev/null
then
    realpath() {
        [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
    }
fi

# Get the path of this script
SCRIPTPATH=$(realpath  $(dirname "$0"))

cd "${SCRIPTPATH}"
cd ..

git submodule update --init --recursive xSqueezeIt
cd xSqueezeIt
cd htslib
autoreconf -i
./configure
make
make install
ldconfig
cd ..


git clone https://github.com/facebook/zstd.git
cd zstd
make
cd ..