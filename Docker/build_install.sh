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

git pull
make

chmod +x pp_extractor/pp_extract
cp pp_extractor/pp_extract /usr/local/bin/
chmod +x phase_caller/phase_caller
cp phase_caller/phase_caller /usr/local/bin
chmod +x pp_update/pp_update
cp pp_update/pp_update /usr/local/bin