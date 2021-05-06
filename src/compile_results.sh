#!/bin/bash
#SBATCH -A def-audginny
#SBATCH --job-name=gtex_compile_results
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1:00:00
#SBATCH -e base_dir_path/slurm/outputs/gtex_compile_results.error
#SBATCH -o base_dir_path/slurm/outputs/gtex_compile_results.out

base_dir=base_dir_path
compile_results_py=${base_dir}/src/compile_results.py

data_dir=${base_dir}/data
tissues_dir=${data_dir}/temp/gtex_tissues
save_dir=${data_dir}/proc

gwas_results_dir=${base_dir}/data/temp/gwas/results
gwas_path=${base_dir}/data/proc/gwas.txt
gwas_save_path=${base_dir}/data/proc/gwas_results.txt
tissues_dir=${base_dir}/data/temp/gtex_tissues
tissues_save_dir=${base_dir}/data/proc
scores_save_path=${base_dir}/data/proc/scores.txt


source ${base_dir}/env/bin/activate

python ${compile_results_py} \
--gwas_results_dir ${gwas_results_dir} \
--gwas_path ${gwas_path} \
--gwas_save_path ${gwas_save_path} \
--tissues_dir ${tissues_dir} \
--tissues_save_dir ${tissues_save_dir} \
--scores_save_path ${scores_save_path}
