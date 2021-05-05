#!/bin/bash
# Find base dir from this file
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
base_dir="$(dirname "$SCRIPT_DIR")"

software_dir=${base_dir}/SOFTWARE

cd ${software_dir}
git clone https://github.com/gkichaev/PAINTOR_V3.0.git

cd PAINTOR_V3.0
bash install.sh
