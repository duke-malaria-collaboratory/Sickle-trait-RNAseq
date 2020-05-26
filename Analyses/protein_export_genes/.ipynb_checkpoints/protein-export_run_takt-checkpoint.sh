#!/bin/sh

# Sickle cell AS
# Normal AA

# python rna_rna_diss_para.py -t1  -t2  -c 10000 -s 15 -p 5 -o 'Results'

# Questions
# are RNA delayed between traits of each parasite: AS 3D7 vs AA 3D7, AS FUP vs AA FUP.

ECHO 3D7
ECHO AA13 vs AA17
python ../Sickle-trait-RNAseq/TAKT/python/rna_rna_diss.py -t1 /Users/robertmoseley/Desktop/Malaria-Sickle/Sickle-trait-RNAseq/Analyses/protein_export_genes/Data/AA13_3D7_PE.txt -t2 /Users/robertmoseley/Desktop/Malaria-Sickle/Sickle-trait-RNAseq/Analyses/protein_export_genes/Data/AA17_3D7_PE.txt -c 10000 -s 15 -p 5 -o 'Results'
ECHO ------------------------------------

ECHO AA13 vs AS15
python ../Sickle-trait-RNAseq/TAKT/python/rna_rna_diss.py -t1 /Users/robertmoseley/Desktop/Malaria-Sickle/Sickle-trait-RNAseq/Analyses/protein_export_genes/Data/AA13_3D7_PE.txt -t2 /Users/robertmoseley/Desktop/Malaria-Sickle/Sickle-trait-RNAseq/Analyses/protein_export_genes/Data/AS15_3D7_PE.txt -c 10000 -s 15 -p 5 -o 'Results'
ECHO ------------------------------------

ECHO AA13 vs AS16
python ../Sickle-trait-RNAseq/TAKT/python/rna_rna_diss.py -t1 /Users/robertmoseley/Desktop/Malaria-Sickle/Sickle-trait-RNAseq/Analyses/protein_export_genes/Data/AA13_3D7_PE.txt -t2 /Users/robertmoseley/Desktop/Malaria-Sickle/Sickle-trait-RNAseq/Analyses/protein_export_genes/Data/AS16_3D7_PE.txt -c 10000 -s 15 -p 5 -o 'Results'
ECHO ------------------------------------

ECHO AA17 vs AS15
python ../Sickle-trait-RNAseq/TAKT/python/rna_rna_diss.py -t1 /Users/robertmoseley/Desktop/Malaria-Sickle/Sickle-trait-RNAseq/Analyses/protein_export_genes/Data/AA17_3D7_PE.txt -t2 /Users/robertmoseley/Desktop/Malaria-Sickle/Sickle-trait-RNAseq/Analyses/protein_export_genes/Data/AS15_3D7_PE.txt -c 10000 -s 15 -p 5 -o 'Results'
ECHO ------------------------------------

ECHO AA17 vs AS16
python ../Sickle-trait-RNAseq/TAKT/python/rna_rna_diss.py -t1 /Users/robertmoseley/Desktop/Malaria-Sickle/Sickle-trait-RNAseq/Analyses/protein_export_genes/Data/AA17_3D7_PE.txt -t2 /Users/robertmoseley/Desktop/Malaria-Sickle/Sickle-trait-RNAseq/Analyses/protein_export_genes/Data/AS16_3D7_PE.txt -c 10000 -s 15 -p 5 -o 'Results'
ECHO ------------------------------------
ECHO ------------------------------------

ECHO FUB
ECHO AA13 vs AA17
python ../Sickle-trait-RNAseq/TAKT/python/rna_rna_diss.py -t1 /Users/robertmoseley/Desktop/Malaria-Sickle/Sickle-trait-RNAseq/Analyses/protein_export_genes/Data/AA13_FUP_PE.txt -t2 /Users/robertmoseley/Desktop/Malaria-Sickle/Sickle-trait-RNAseq/Analyses/protein_export_genes/Data/AA17_FUP_PE.txt -c 10000 -s 15 -p 5 -o 'Results'
ECHO ------------------------------------

ECHO AA13 vs AS18
python ../Sickle-trait-RNAseq/TAKT/python/rna_rna_diss.py -t1 /Users/robertmoseley/Desktop/Malaria-Sickle/Sickle-trait-RNAseq/Analyses/protein_export_genes/Data/AA13_FUP_PE.txt -t2 /Users/robertmoseley/Desktop/Malaria-Sickle/Sickle-trait-RNAseq/Analyses/protein_export_genes/Data/AS18_FUP_PE.txt -c 10000 -s 15 -p 5 -o 'Results'
ECHO ------------------------------------

ECHO AA13 vs AS19
python ../Sickle-trait-RNAseq/TAKT/python/rna_rna_diss.py -t1 /Users/robertmoseley/Desktop/Malaria-Sickle/Sickle-trait-RNAseq/Analyses/protein_export_genes/Data/AA13_FUP_PE.txt -t2 /Users/robertmoseley/Desktop/Malaria-Sickle/Sickle-trait-RNAseq/Analyses/protein_export_genes/Data/AS19_FUP_PE.txt -c 10000 -s 15 -p 5 -o 'Results'
ECHO ------------------------------------

ECHO AA17 vs AS18
python ../Sickle-trait-RNAseq/TAKT/python/rna_rna_diss.py -t1 /Users/robertmoseley/Desktop/Malaria-Sickle/Sickle-trait-RNAseq/Analyses/protein_export_genes/Data/AA17_FUP_PE.txt -t2 /Users/robertmoseley/Desktop/Malaria-Sickle/Sickle-trait-RNAseq/Analyses/protein_export_genes/Data/AS18_FUP_PE.txt -c 10000 -s 15 -p 5 -o 'Results'
ECHO ------------------------------------

ECHO AA17 vs AS19
python ../Sickle-trait-RNAseq/TAKT/python/rna_rna_diss.py -t1 /Users/robertmoseley/Desktop/Malaria-Sickle/Sickle-trait-RNAseq/Analyses/protein_export_genes/Data/AA17_FUP_PE.txt -t2 /Users/robertmoseley/Desktop/Malaria-Sickle/Sickle-trait-RNAseq/Analyses/protein_export_genes/Data/AS19_FUP_PE.txt -c 10000 -s 15 -p 5 -o 'Results'
