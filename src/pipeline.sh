#!/bin/bash
#SBATCH -A def-audginny
#SBATCH --job-name=gtex_pipeline
#SBATCH --mem-per-cpu=40G
#SBATCH --time=1:00:00
#SBATCH -e base_dir_path/slurm/outputs/pipeline.error
#SBATCH -o base_dir_path/slurm/outputs/pipeline.out

# Change the last two lines for your own jobs!

######################
##### Parameters #####
######################

base_dir=base_dir_path

src_dir=${base_dir}/src
data_dir=${base_dir}/data
data_raw_dir=${data_dir}/raw
proc_data_dir=${data_dir}/proc
data_temp_dir=${data_dir}/temp

jobs_dir=${base_dir}/slurm/jobs
error_dir=${base_dir}/slurm/outputs
out_dir=${base_dir}/slurm/outputs

# Gwas Prep
gwas_prep_py=${src_dir}/prepare_gwas.py
raw_gwas_path=${data_raw_dir}/raw_GWAS.txt #MCP_Eur_A21.tsv
breaks_path=${data_raw_dir}/risk_loci.txt
gwas_save_path=${proc_data_dir}/gwas.txt

# Gwas fmap Prep
gwas_fmap_temp=${data_dir}/temp/gwas
fmap_prep_path=${src_dir}/cmd_fmap.py
G_dir=${data_dir}/raw/plink_files
G_pattern='1000G.EUR.QC.*.bed'
gwas_results_dir=${gwas_fmap_temp}/results

# Gwas finemap jobs
fmap_method='paintor'
gwas_job='gwas_fmap'
paintor_path=${base_dir}/SOFTWARE/PAINTOR_V3.0/PAINTOR

# GTEx pipeline


gtex_ref_dir=${proc_data_dir}/GTEx_v7_ref

raw_gene_ref=${gtex_data_dir}/gencode.v19.genes.v7.patched_contigs.gtf
raw_snps_ref=${gtex_data_dir}/GTEx_ref.txt


gene_ref_path=${gtex_ref_dir}/ref_genes.txt.gz
snps_ref_path=${gtex_ref_dir}/ref_snps.txt.gz

tissue_preproc_save_dir=${tissue_dir}/gene_groups
tissue_fmap_dir=${tissue_dir}/fmap_files
tissues_temp_dir=${data_temp_dir}/gtex_tissues




# # Drg preproc
# drg_preproc_py=${src_dir}/drg_preproc.py
# raw_drg_path=${data_raw_dir}/DRG_raw.txt
# n_groups=1000
# gene_groups_dir=${proc_data_dir}/gene_groups

# job creator
# job_creator_sh=${src_dir}/job_creator.sh

####################################
##### Activate job environment #####
####################################

source ${base_dir}/env/bin/activate

################
##### GWAS #####
################

####### Commenting ######
##### Clean GWAS #####
python ${gwas_prep_py} \
  --ss_path ${raw_gwas_path} \
  --breaks_path ${breaks_path} \
  --save_path ${gwas_save_path}

##### Prepare GWAS #####
python ${fmap_prep_path} \
  --method ${fmap_method} \
  --save_dir ${gwas_fmap_temp} \
  --ss_path ${gwas_save_path} \
  --G_dir ${G_dir} \
  --G_pattern ${G_pattern}

##### Send GWAS fine-mapping job #####
mkdir -p ${gwas_results_dir}
cat > ${jobs_dir}/${gwas_job}.job << EOT
#!/bin/bash
#SBATCH --job-name=${gwas_job}.job
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=10G
#SBATCH -t 10:00:00
#SBATCH -A def-audginny
#SBATCH --error=${error_dir}/${gwas_job}.error
#SBATCH --output=${out_dir}/${gwas_job}.out

echo 'Fine-mapping'
if [[ ${fmap_method} == 'finemap' ]]; then
  echo 'finemap not defined yet'

elif [[ '${fmap_method}' == 'paintor' ]]; then
  echo 'Fine-mapping with Paintor'
  ${paintor_path} \
  -input ${gwas_fmap_temp}/input.files \
  -in ${gwas_fmap_temp} \
  -out ${gwas_fmap_temp}/results \
  -Zhead z \
  -LDname LD \
  -annotations base \
  -mcmc
fi
echo 'done'
EOT
sbatch ${jobs_dir}/${gwas_job}.job

################
##### GTEx #####
################
# Ref files
gtex_prep_ref_py=${src_dir}/prep_ref.py
gtex_raw_data_dir=${data_raw_dir}/GTEx_v7
gtex_tissues_dir=${gtex_raw_data_dir}/tissues

raw_gene_ref=${gtex_raw_data_dir}/gencode.v19.genes.v7.patched_contigs.gtf
raw_snps_ref=${gtex_raw_data_dir}/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz
ref_save_dir=${base_dir}/data/proc

# pre-process tissues
gtex_preproc_py=${src_dir}/gtex_preproc.py

# Process tissue gene group
gtex_proc_py=${src_dir}/gtex_proc.py

# 1. Create ref-files
python ${gtex_prep_ref_py} \
  --gene_ref ${raw_gene_ref} \
  --snps_ref ${raw_snps_ref} \
  --save_dir ${ref_save_dir}

# 2. Create Tissue job
for tissue_path in ${gtex_tissues_dir}/*.txt
do

echo ${tissue_path}
tissue_name=${tissue_path##*/}
tissue_name="${tissue_name%.allpairs.txt}"

echo "pre-processing ${tissue_name}"

cat > ${jobs_dir}/${tissue_name}.job << EOT
#!/bin/bash
#SBATCH --job-name=${tissue_name}.job
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=40G
#SBATCH -t 1:00:00
#SBATCH -A def-audginny
#SBATCH --error=${error_dir}/${tissue_name}.error
#SBATCH --output=${out_dir}/${tissue_name}.out

source ${base_dir}/env/bin/activate

python ${gtex_preproc_py} \
--gene_ref_path ${ref_save_dir}/ref_genes.txt.gz \
--snps_ref_path ${ref_save_dir}/ref_snps.txt.gz \
--gtex_path ${tissue_path} \
--gwas_path ${gwas_save_path} \
--ngroups 100 \
--save_dir ${tissues_temp_dir}/${tissue_name}/gene_groups


for gene_group_path in ${tissues_temp_dir}/${tissue_name}/gene_groups/*.txt
do

group_name=\${gene_group_path##*/}
group_name="\${group_name%.txt}"

echo "\${group_name}"
job_name=${tissue_name}_\${group_name}
cat > ${jobs_dir}/\${job_name}.job << geneEOT
#!/bin/bash
#SBATCH --job-name=\${job_name}.job
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=40G
#SBATCH -t 10:00:00
#SBATCH -A def-audginny
#SBATCH --error=${error_dir}/\${job_name}.error
#SBATCH --output=${out_dir}/\${job_name}.out

source ${base_dir}/env/bin/activate

# gene-group processing
python ${gtex_proc_py} \
--gene_group_path \${gene_group_path} \
--gtex_path ${tissue_path} \
--snps_ref_path ${ref_save_dir}/ref_snps.txt.gz \
--save_path ${tissues_temp_dir}/${tissue_name}/\${group_name}/\${group_name}.ss.txt

# gene-group fmap-prep
python ${fmap_prep_path} \
--method '${fmap_method}' \
--save_dir ${tissues_temp_dir}/${tissue_name}/\${group_name}/ \
--ss_path ${tissues_temp_dir}/${tissue_name}/\${group_name}/\${group_name}.ss.txt \
--G_dir ${G_dir} \
--G_pattern ${G_pattern}

mkdir -p ${tissues_temp_dir}/${tissue_name}/\${group_name}/results

echo 'Fine-mapping'
if [[ ${fmap_method} == 'finemap' ]]; then
echo 'Fine-mapping with FINEMAP'
${finemap_path} --sss --in-files ${tissues_temp_dir}/${tissue_name}/${group_name}/master

elif [[ '${fmap_method}' == 'paintor' ]]; then
echo 'Fine-mapping with Paintor'

${paintor_path} \
-input ${tissues_temp_dir}/${tissue_name}/\${group_name}/input.files \
-in ${tissues_temp_dir}/${tissue_name}/\${group_name} \
-out ${tissues_temp_dir}/${tissue_name}/\${group_name}/results \
-Zhead z \
-LDname LD \
-annotations base \
-mcmc
fi
geneEOT

sbatch ${jobs_dir}/\${job_name}.job
done
EOT
sbatch ${jobs_dir}/${tissue_name}.job
done
