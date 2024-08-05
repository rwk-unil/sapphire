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

git submodule update --init --recursive xSqueezeIt
cd xSqueezeIt
cd htslib
autoreconf -i
./configure
make -j$(nproc)
cd ..

git clone https://github.com/facebook/zstd.git
cd zstd
make -j$(nproc)
cd ..
cd ..

touch .dependencies_ready