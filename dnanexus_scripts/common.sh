#/bin/bash

DOCKER_IMAGE=docker/pp_toolkit_v1.4.tar.gz

if ! command -v jq &> /dev/null
then
    echo "Please install jq"
    echo "E.g., sudo apt install jq"
    exit 1
fi

dx_id_to_name () {
    echo "$(dx describe --json "$1" | jq -r '.name')"
}

dx_id_to_path () {
    FFNAME=$(dx describe --json "$1" | jq -r '.name')
    FFPATH=$(dx describe --json "$1" | jq -r '.folder')
    echo "/mnt/project/${FFPATH}/${FFNAME}"
}

dx_id_to_dx_path () {
    echo $(dx describe --json "$1" | jq -r '.folder')
}

dx_id_to_dx_path_and_name () {
    FFNAME=$(dx describe --json "$1" | jq -r '.name')
    FFPATH=$(dx describe --json "$1" | jq -r '.folder')
    echo "${FFPATH}/${FFNAME}"
}

path_to_dx_id () {
    echo "$(dx describe --json "$1" | jq -r '.id')"
}

ask_permission_to_launch() {
    while true; do
        read -p "Do you want to launch on DNANexus? [y/n]" yn
        case $yn in
            y)
            echo "Launching !";
            break
            ;;
            n)
            echo "exiting...";
            exit
            ;;
            *)
            echo "unexpected input"
            ;;
        esac
    done
}

ask_permission_to_launch_all() {
    while true; do
        read -p "Do you want to launch on DNANexus? [y/n/a/e]" yn
        case $yn in
            y)
            echo "Launching !";
            break
            ;;
            n)
            echo "Skippping !";
            SKIP_JOB="yes"
            break
            ;;
            e)
            echo "exiting...";
            exit
            ;;
            a)
            echo "Launching all !";
            LAUNCH_ALL="yes";
            break
            ;;
            *)
            echo "unexpected input"
            ;;
        esac
    done
}

ask_permission_to_launch_message() {
    while true; do
        read -p "Do you want to $1? [y/n/a/e]" yn
        case $yn in
            y)
            echo "Launching !";
            break
            ;;
            n)
            echo "Skippping !";
            SKIP_JOB="yes"
            break
            ;;
            e)
            echo "exiting...";
            exit
            ;;
            a)
            echo "Launching all !";
            LAUNCH_ALL="yes";
            break
            ;;
            *)
            echo "unexpected input"
            ;;
        esac
    done
}

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

# Command line argument parsing from :
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
    --cost-limit)
    COST_LIMIT="$2"
    shift # past argument
    shift # past value
    ;;
    --instance)
    INSTANCE="$2"
    shift # past argument
    shift # past value
    ;;
    --docker-image)
    DOCKER_IMAGE="$2"
    shift # past argument
    shift # past value
    ;;
    --verbose)
    VERBOSE="-v"
    shift # no value attached
    ;;
    --tag)
    TAG="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

COST_LIMIT_ARG=""
if ! [ -z "${COST_LIMIT}" ]
then
    COST_LIMIT_ARG="--cost-limit ${COST_LIMIT}"
fi
