#!/bin/tcsh
#
##SBATCH -p common                # Partition to submit to (comma separated)
#SBATCH -J trim_align         # Job name
#SBATCH -n 16                     # Number of cores
#SBATCH -N 1                     # Ensure that all cores are on one machine
#SBATCH -t 2-00:00                # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 50000               # Memory in MB
#SBATCH -o _FUPtp1_align_%j.out # File for STDOUT (with jobid = %j) 
#SBATCH -e FUP-tp1_%j.err       # File for STDERR (with jobid = %j) 
#SBATCH --mail-type=ALL          # Type of email notification: BEGIN,END,FAIL,ALL
#SBATCH --mail-user=jws48@duke.edu  # Email where notifications will be sent
#Your actual work goes after this line


module load STAR/2.5.0c-fasrc01
module load samtools/1.3.1-gcb01
module load cufflinks/2.2.1-fasrc01 

./align_trim.pl /gpfs/fs1/data/taylorlab/Genomes/3D7_2019/FASTA /gpfs/fs1/data/taylorlab/HbAS/FUP/tp1/out /gpfs/fs1/data/taylorlab/HbAS/FUP/reads/tp1/1 /gpfs/fs1/data/taylorlab/HbAS/FUP/reads/tp1/2 /gpfs/fs1/data/taylorlab/Genomes/Homo_sapiens/FASTA/ /gpfs/fs1/data/taylorlab/Genomes/3D7_2019/GTF/Plasmodium_falciparum.EPr1.43.gtf
