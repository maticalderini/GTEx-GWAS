# Find base dir from this file
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
base_dir="$(dirname "$SCRIPT_DIR")"

base_url=https://alkesgroup.broadinstitute.org/LDSCORE
plink_3=1000G_Phase3_plinkfiles.tgz
data_raw_dir=${base_dir}/data/raw

cd ${data_raw_dir}
echo ${data_raw_dir}
wget ${base_url}/${plink_3} && tar -xvzf ${plink_3} && mv 1000G_EUR_Phase3_plink plink_files && rm ${plink_3}
