#!/bin/bash
#
##SBATCH -p common                # Partition to submit to (comma separated)
#SBATCH -J takt         # Job name
#SBATCH -n 1                     # Number of cores
#SBATCH -N 1                     # Ensure that all cores are on one machine
#SBATCH -t 2-00:00                # Runtime in D-HH:MM (or use minutes)
#SBATCH --mem 2000               # Memory in MB
#SBATCH -o _TAKT_%j.out # File for STDOUT (with jobid = %j) 
#SBATCH -e TAKT_%j.err       # File for STDERR (with jobid = %j) 
#SBATCH --mail-type=ALL          # Type of email notification: BEGIN,END,FAIL,ALL
#SBATCH --mail-user=jws48@duke.edu  # Email where notifications will be sent
#Your actual work goes after this line
module load Anaconda3/4.3.0-gcb01
source /data/taylorlab/jws48-conda/envs/takt/bin/activate

echo 3D7
echo AA13 vs AA17
python  /gpfs/fs1/home/jws48/TAKT/python/rna_rna_diss.py -t1 /gpfs/fs1/home/jws48/TAKT/DTW_genes/Data/AA13_thresh_3D7_DTW.txt -t2 /gpfs/fs1/home/jws48/TAKT/DTW_genes/Data/AA17_thresh_3D7_DTW.txt -c 10000 -s 15 -o /gpfs/fs1/home/jws48/TAKT/DTW_genes/Results
echo ------------------------------------

echo 3D7
echo AA13 vs AS15
python  /gpfs/fs1/home/jws48/TAKT/python/rna_rna_diss.py -t1 /gpfs/fs1/home/jws48/TAKT/DTW_genes/Data/AA13_thresh_3D7_DTW.txt -t2 /gpfs/fs1/home/jws48/TAKT/DTW_genes/Data/AS15_thresh_3D7_DTW.txt -c 10000 -s 15 -o /gpfs/fs1/home/jws48/TAKT/DTW_genes/Results
echo ------------------------------------

echo 3D7
echo AA13 vs AS16
python  /gpfs/fs1/home/jws48/TAKT/python/rna_rna_diss.py -t1 /gpfs/fs1/home/jws48/TAKT/DTW_genes/Data/AA13_thresh_3D7_DTW.txt -t2 /gpfs/fs1/home/jws48/TAKT/DTW_genes/Data/AS16_thresh_3D7_DTW.txt -c 10000 -s 15 -o /gpfs/fs1/home/jws48/TAKT/DTW_genes/Results
echo ------------------------------------

echo 3D7
echo AA17 vs AS16
python  /gpfs/fs1/home/jws48/TAKT/python/rna_rna_diss.py -t1 /gpfs/fs1/home/jws48/TAKT/DTW_genes/Data/AA17_thresh_3D7_DTW.txt -t2 /gpfs/fs1/home/jws48/TAKT/DTW_genes/Data/AS15_thresh_3D7_DTW.txt -c 10000 -s 15 -o /gpfs/fs1/home/jws48/TAKT/DTW_genes/Results
echo ------------------------------------

echo 3D7
echo AA17 vs AS16
python  /gpfs/fs1/home/jws48/TAKT/python/rna_rna_diss.py -t1 /gpfs/fs1/home/jws48/TAKT/DTW_genes/Data/AA17_thresh_3D7_DTW.txt -t2 /gpfs/fs1/home/jws48/TAKT/DTW_genes/Data/AS16_thresh_3D7_DTW.txt -c 10000 -s 15 -o /gpfs/fs1/home/jws48/TAKT/DTW_genes/Results
echo ------------------------------------

echo FUP
echo AA13 vs AA17
python  /gpfs/fs1/home/jws48/TAKT/python/rna_rna_diss.py -t1 /gpfs/fs1/home/jws48/TAKT/DTW_genes/Data/AA13_thresh_FUP_DTW.txt -t2 /gpfs/fs1/home/jws48/TAKT/DTW_genes/Data/AA17_thresh_FUP_DTW.txt -c 10000 -s 15 -o /gpfs/fs1/home/jws48/TAKT/DTW_genes/Results
echo ------------------------------------

echo FUP
echo AA13 vs AS18
python  /gpfs/fs1/home/jws48/TAKT/python/rna_rna_diss.py -t1 /gpfs/fs1/home/jws48/TAKT/DTW_genes/Data/AA13_thresh_FUP_DTW.txt -t2 /gpfs/fs1/home/jws48/TAKT/DTW_genes/Data/AS18_thresh_FUP_DTW.txt -c 10000 -s 15 -o /gpfs/fs1/home/jws48/TAKT/DTW_genes/Results
echo ------------------------------------

echo FUP
echo AA13 vs AS19
python  /gpfs/fs1/home/jws48/TAKT/python/rna_rna_diss.py -t1 /gpfs/fs1/home/jws48/TAKT/DTW_genes/Data/AA13_thresh_FUP_DTW.txt -t2 /gpfs/fs1/home/jws48/TAKT/DTW_genes/Data/AS19_thresh_FUP_DTW.txt -c 10000 -s 15 -o /gpfs/fs1/home/jws48/TAKT/DTW_genes/Results
echo ------------------------------------

echo FUP
echo AA17 vs AS18
python  /gpfs/fs1/home/jws48/TAKT/python/rna_rna_diss.py -t1 /gpfs/fs1/home/jws48/TAKT/DTW_genes/Data/AA17_thresh_FUP_DTW.txt -t2 /gpfs/fs1/home/jws48/TAKT/DTW_genes/Data/AS18_thresh_FUP_DTW.txt -c 10000 -s 15 -o /gpfs/fs1/home/jws48/TAKT/DTW_genes/Results
echo ------------------------------------

echo FUP
echo AA17 vs AS19
python  /gpfs/fs1/home/jws48/TAKT/python/rna_rna_diss.py -t1 /gpfs/fs1/home/jws48/TAKT/DTW_genes/Data/AA17_thresh_FUP_DTW.txt -t2 /gpfs/fs1/home/jws48/TAKT/DTW_genes/Data/AS19_thresh_FUP_DTW.txt -c 10000 -s 15 -o /gpfs/fs1/home/jws48/TAKT/DTW_genes/Results
echo ------------------------------------


