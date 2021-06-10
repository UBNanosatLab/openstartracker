#!/bin/bash

# requires docker to be installed

# build the Dockerfile to include the current director and have internet access
docker build --network=host -t startracker1 .

# directory of this setup.sh script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Make an alias to start the docker environment
#   in interactive mode
#   in the proper directory
alias dstart='docker run -it -v '$DIR'/:/home startracker1'
