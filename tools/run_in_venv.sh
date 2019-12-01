#!/usr/bin/env bash
# setup venv
CURRENT_DIR="$(cd "$(dirname "$0")" && pwd)"

ROOT_DIR=$CURRENT_DIR/../
VENV_DIR=$CURRENT_DIR/../venv/

source ${VENV_DIR}/bin/activate
export PYTHONPATH=${BACK_DIR}:$PYTHONPATH
python $1 "${@:2}"