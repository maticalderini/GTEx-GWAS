SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
base_dir="$(dirname "$SCRIPT_DIR")"
src_dir=${base_dir}/src

##########################
##### Download files #####
##########################
echo "Downloading files"
echo "Downloading Plink files"
bash ${src_dir}/get_plinkfiles.sh

echo "Downloading PAINTOR"
bash ${src_dir}/get_paintor.sh

echo "Downloading gtex tissue data"
bash ${src_dir}/get_GTEx_data.sh

######################
##### Virtualenv #####
######################
echo "Installing python virtual environement"
module load python/3.7
virtualenv --no-download ${base_dir}/env
source ${base_dir}/env/bin/activate
pip install --no-index --upgrade pip

pip install --no-index -r ${base_dir}/requirements.txt


echo "Setup done. Don't forget to add the raw GWAS and risk loci files"
