#!/bin/tcsh
#
##SBATCH -p common                # Partition to submit to (comma separated)
#SBATCH -J salmon         # Job name
#SBATCH -n 16                     # Number of cores
#SBATCH -N 1                     # Ensure that all cores are on one machine
#SBATCH -t 2-00:00                # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 50000               # Memory in MB
#SBATCH -o _salmonQuant_%j.out # File for STDOUT (with jobid = %j) 
#SBATCH -e salmonQuant_%j.err       # File for STDERR (with jobid = %j) 
#SBATCH --mail-type=ALL          # Type of email notification: BEGIN,END,FAIL,ALL
#SBATCH --mail-user=jws48@duke.edu  # Email where notifications will be sent
#Your actual work goes after this line


./salmonQuant.pl /gpfs/fs1/data/taylorlab/Genomes/3D7_ASM276v2/cDNA/pfal_ASM276v2_index /gpfs/fs1/data/taylorlab/HbAS/salmonQuant/out /gpfs/fs1/data/taylorlab/HbAS/HsCounts/3D7/reads/1 /gpfs/fs1/data/taylorlab/HbAS/HsCounts/3D7/reads/2 /gpfs/fs1/data/taylorlab/Genomes/Homo_sapiens/cDNA/Homo_sapiens_GRCh38_index /gpfs/fs1/data/taylorlab/Genomes/3D7_ASM276v2/GTF/Plasmodium_falciparum.ASM276v2.45.gtf
