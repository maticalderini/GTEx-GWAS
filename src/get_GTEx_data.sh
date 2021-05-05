# Find base dir from this file
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
base_dir="$(dirname "$SCRIPT_DIR")"


raw_dir=${base_dir}/data/raw
echo ${raw_dir}

tissues_file=${raw_dir}/GTEx_tissues_list.txt

gtex_dir=${raw_dir}/GTEx_v7
tissues_dir=${gtex_dir}/tissues
ref_dir=${gtex_dir}

mkdir -p ${tissues_dir}
cd ${tissues_dir}
wget -N --input ${tissues_file}

mkdir -p ${ref_dir}
cd ${ref_dir}
wget -N https://storage.googleapis.com/gtex_analysis_v7/reference/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz
wget -N https://storage.googleapis.com/gtex_analysis_v7/reference/gencode.v19.genes.v7.patched_contigs.gtf

echo 'Unzipping tissue files (needed for pipeline)'
for pathname in ${tissues_dir}/*.txt.gz
do

  [[ -e ${pathname%.gz} ]] || (echo "unzipping ${pathname}" && gunzip -kd "${pathname}")
done
echo 'Done'
