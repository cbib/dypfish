#!/usr/bin/env bash
# setup file for local venv
# setup venv
CURRENT_DIR="$(cd "$(dirname "$0")" && pwd)"

ROOT_DIR=$CURRENT_DIR/../
VENV_DIR=$CURRENT_DIR/../venv/

# some pb with the default virtualenv install not working, had to force upgrade
# create the venv
cd ${ROOT_DIR} && virtualenv -p python3.7 venv --prompt="(DYPFISH)"

# setup in the venv

source ${VENV_DIR}/bin/activate && pip3.7 install -r ${ROOT_DIR}/requirements.txt
