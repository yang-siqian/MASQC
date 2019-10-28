# MASQC
Next Generation Sequencing assists Third Generation Sequencing for quality control in N6-Methyladenine DNA identification
## Description
DNA N6-methyladenine (6mA) modification is the most prevalent DNA modification in prokaryotes and eukaryotes, involving gene expression, DNA replication and repair, and host-pathogen interactions. Single-molecule real-time sequencing (SMRT-seq) can detect 6mA events at single-nucleotide resolution in prokaryotic genomes, but its application to detecting 6mA in eukaryotic genomes still faces great challenges. We analyzed the distribution of 6mA methylomes in both eukaryotic and prokaryotes, and proposed MASQC method (MeDIP-seq assists SMRT-seq for quality control in 6mA identification), to identify 6mA events without doing the whole genome amplification (WGA) sequencing. The proposed MASQC method was evaluated on two eukaryotic genomes and six bacterial genomes, our results demonstrate that MASQC performs well in quality control of false positive 6mA identification for both eukaryotic and prokaryotic genomes.
## Dependencies
- PacBio SMRT-analysis-2.3.0
- BWA
- MACS2
- R-3.5 
- Python-3.6
## RUN
Please read the scripts in detail and modify the path and other configures ccording to your demand.
### 1.preprocess
```
# generate modification.gff
bash smrt-analysis.sh

# generate peak file
bash macs2_peak.sh
```
### 2.threshold_scripts
```
# modify config.ini and species_table.txt
bash run.sh
```
### 3.extract_info
```
# modify species.txt and path in generate_csv.py
bash csv.sh
```
## Desciption of generated file in step 2
#### $species_dataset.txt
```
chr | m6A_loc | score | coverage | context | frac | IPDratio | enrichment
```
#### $species_peakfilter.txt
- filter out enrichment<1,-log10(pvalue)<2,-log10(qvalue)<2 peaks
####  $species_modfilter.txt
- filter out score<30,coverage<50 6mA_modifications
#### mod_thres_filter.txt
total m6a_modifications after threshold filtering
#### $species_$motif_motifsall.txt
- all motifs
```
chr | m6A | m6A_loc | m6A_loc | score | coverage | context | frac | IPDratio
```
#### $species_$motif_motifsinpeak.txt
- motifs in peaks
```
chr | m6A | m6A_loc | m6A_loc | score | coverage | context | frac | IPDratio
```
#### $species_info.txt
- the information file (used in step 3)

