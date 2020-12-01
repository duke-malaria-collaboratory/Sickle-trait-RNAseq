#!/bin/bash
#SBATCH -e /hpc/group/haaselab/rcm49/sickle-trait-rnaseq/TAKT_DCC/results/takt_ring_stage_output/FUP/takt_FUP_ring.err
#SBATCH -o /hpc/group/haaselab/rcm49/sickle-trait-rnaseq/TAKT_DCC/results/takt_ring_stage_output/FUP/takt_FUP_ring.out
#SBATCH -a 1-4657
#SBATCH -J FUP_ring
#SBATCH --mem=2G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=robert.moseley@duke.edu

module load Python/3.8.1

ring_file_dir='/hpc/group/haaselab/rcm49/sickle-trait-rnaseq/TAKT_DCC/data/ring_stage_data'
ring_out_dir='/hpc/group/haaselab/rcm49/sickle-trait-rnaseq/TAKT_DCC/results/takt_ring_stage_output'

python -u /hpc/group/haaselab/rcm49/sickle-trait-rnaseq/TAKT/python/rna_rna_diss_slurm.py -fd ${ring_file_dir} -c 720 -s 5 -o ${ring_out_dir}
