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

# To get the latest version
git pull
make STATIC_BINS=y

chmod +x pp_extractor/pp_extract
cp pp_extractor/pp_extract /usr/local/bin/
chmod +x bin_tools/bin_splitter
cp bin_tools/bin_splitter /usr/local/bin/
chmod +x bin_tools/bin_merger
cp bin_tools/bin_merger /usr/local/bin/
chmod +x phase_caller/phase_caller
cp phase_caller/phase_caller /usr/local/bin
chmod +x pp_update/pp_update
cp pp_update/pp_update /usr/local/bin
chmod +x bin_tools/vertical_bin_merger
cp bin_tools/vertical_bin_merger /usr/local/bin