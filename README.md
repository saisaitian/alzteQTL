
teQTL data analysis pipeline

1) 01_TE_Expression_Quantification

  This pipeline is used to identify differentially expressed TE between AD cases and cognitive healthy controls. It contains four steps:
  
  (1) step01_squire_fetch.sh: Downloads input files from RefGene and generates STAR index
  
  (2) step02_squire_clean.sh: Filters Repeatmasker file for Repeats of interest, collapses overlapping repeats, and returns as BED file
  
  (3) step03_squire_Map_Count_Single.sh: Aligns single-end RNAseq data to reference genome index, you can run this file by 'qsub -cwd step03_squire_Map_Count_Single.sh raw_dir sample_name'. raw_dir refer to the directory of your raw fastq files
  
  (4) step03_squire_Map_Count_Paired.sh:Aligns paired-end RNAseq data to reference genome index, you can run this file by 'qsub -cwd step03_squire_Map_Count_Paired.sh raw_dir sample_name'. raw_dir refer to the directory of your raw fastq files
  
  (5) step04_squire_Call_sf.sh: Quantifies RNAseq reads aligning to subfamily TEs and genes. Use 'qsub -cwd step04_squire_Call_sf.sh AD_subjects HC_subjects' to run this script
  
  (6) step04_squire_Call_locus.sh: Quantifies RNAseq reads aligning to locus TEs and genes. Use 'qsub -cwd step04_squire_Call_locus.sh AD_subjects HC_subjects' to run this script


2) Identification of cis-teQTLs using MatrixEQTL
  
  run this script for each chromosome, use chr1 as an example:
  nohup Rscript Matrix_eQTL_pepline.R chr1_genotyping_dir chr1_transcript_dir output_dir chr1 > chr1.log 2>&1 &


3) colocalization analysis between GWAS loci and xQTLs
  Prepare GWAS data and xQTL data first:
  
  The column names of GWAS data: "snp","chr","position","ref","alt","maf","gene","beta","varbeta","pvalue".
  The column names of xQTL data: "snp", "chr", "position", "ref", "alt", "maf", "probes", "beta", "varbeta", "pvalue". 

  run this script using: nohup Rscript coloc.R GWAS xQTL > log 2>&1 &

4) Mayo_teQTL_199samples.txt.gz: teQTLs identified from Mayo brain biobank. This file contains teQTLs with P < 1E-5. Coordinates shown in this file are annotated with GRCh38.

5) ROSMAP_teQTL_57samples.txt.gz: teQTLs identified from Mayo brain biobank. This file contains teQTLs with P < 1E-5. Coordinates shown in this file are annotated with GRCh38.



