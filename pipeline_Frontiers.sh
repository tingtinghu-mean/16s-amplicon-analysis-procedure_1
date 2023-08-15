[TOC]

# Reference
# Name Name: Easy Amplicon (EasyAmplicon)
# Author Authors: Liu Yongxin (Yong-Xin Liu), Chen Tong (Tong Chen), Zhou Xin (Xin Zhou), etc

# 1. Set up process EasyAmplicon (ea) and working directory (work directory, wd), add environment variables and enter wd
# ** The following 4 lines must run each time Rstudio * *
ea=/d/gmu/EasyAmplicon
wd=/d/gmu/lung part1/
PATH=$PATH:${ea}/win/
cd ${wd}

# Amplicon analysis procedure Analysis pipeline for 16S amplicon

# Preparation before the run
head -n2 metadata.txt
# Raw sequencing data are deposited in the seq catalog, usually in `.fq. End of gz `, one pair of files per sample
ls -sh seq

# Create temporary and result directories, which can be deleted after the temporary directory analysis
mkdir -p temp result

### 1.1. metadata. Experimental design files of txt

# cat View the first 3 rows, -A shows symbols
cat -A metadata.txt | head -n3
# windows If the user has ^ M at the end, run the sed command to remove, and cat-A checks the results
sed -i 's/\r/\n/' metadata.txt
cat -A metadata.txt | head -n3

### 1.2. seq/*.fq. And gz the raw sequencing data

# Sequencing results returned by companies, usually a pair of fq / fastq. Compressed file in gz format

### 1.3. pipeline. The sh process depends on the database

# The first use of the database requires decompression, and the decompression can skip this segment

## 2. merged double-end sequences and renamed Merge paired reads and label samples after the sample

# Batch-process and merge according to the experimental design.
# Table-n + 2 to the header, cut-f 1 take the first column, to obtain the sample list. The merging of 18 samples of 50,000 pairs of sequences took about 2 minutes.
# Windows Copy command Ctrl + C terminate under Linux, to prevent abnormal interruption, end add & turn the command to the background.

for i in `tail -n+2 metadata.txt | cut -f 1`;do
  Vsearch --fastq_mergepairs seq/${i}_1.fastq.gz --reverse seq/${i}_2.fastq.gz --fastqout temp/${i}.merged.fq --relabel ${i}.
done &

# View the file sequence name
head -n1 seq/L_1_1.fastq.gz

# Sequence is named by sample and exported to a new folder
for i in `tail -n+2 metadata2.txt | cut -f 1`;do
  vsearch --fastq_convert seq/${i}.fastq.gz --fastqout temp/${i}.merged.fq --relabel ${i}.
done

# View the converted 33 encoding format, the quality value is mostly capital letters
head -n1 FAQ/relabel/L_1_1.fq

### 3. primer excision and quality control Cut primers and quality filter

# The primer is excised in the same length. If the label is available, the length of the label needs to be counted.
# As shown in this case, the left end tag 10 bp + forward primer V5 19 bp has a total of 29 bp and the reverse primer V7 18 bp.
# Cut barcode 10bp + V5 19bp in left and V7 18bp in right
# Note: Be sure to know the primer and tag length in the experimental design. If the primer has been removed, fill in 0 at the parameters below to indicate that no removal is required.760,000 sequences 37s- -fastq _ maxlen 400- -fastq _ minlen 150

for i in `tail -n+2 metadata.txt | cut -f 1`;do
  vsearch --fastq_filter temp/${i}.merged.fq --fastq_stripleft 23 --fastq_stripright 27 -fastq_maxee_rate 0.01 --fastaout temp/${i}.filtered.fa --fastq_qmax 42
done

# Combine all the samples to the same file
cat temp/*.filtered.fa > temp/filtered.fa
cat temp/NC/*.filtered.fa > temp/Wall.filtered.fa

# View file size 634M, results from different versions of the software
ls -lsh temp/filtered.fa

# View the sequence name. Whether it is the sample name before, and the sample name is never allowed.
head -n 6 temp/filtered.fa|cut -c1-60

# View the file for the fa file format
head temp/filtered.fa
less temp/filtered. Fa | grep '>' -c 
# View number of sequences # usearch-fastx _ relabel seq / $ {i}.fq -fastqout temp/${i}.merged.fq -prefix ${i}.
# less temp/Wall.filtered.fa|grep '>' -c
# Another approach to refer to "FAQ 2"

### 4. Sequence redundancy and select representative sequence (OTU / ASV) Dereplicate and cluster / denoise

### 4.1 sequence deredundancy Dereplication
vsearch --derep_fulllength temp/filtered.fa \
  --output temp/uniques.fa \
  --relabel Uni --minuniquesize 10 --sizeout
# High-abundance non-redundant sequences are very small (<2 Mb), with size and frequency after the name
ls -lsh temp/uniques.fa
head -n 2 temp/uniques.fa
less temp/uniques.fa|grep '>' -c

### 4.2 cluster OTU / denoising ASV Cluster OTUs / denoise ASV

# Redundancy data decreases by at least 1 order of magnitude, reduces the workload of downstream analysis, and is more suitable for identifying true feature sequences based on abundance, such as operational taxonomic units (OTUs) or amplified sequence changes (ASVs).
# Parameter- -miniuniqusize setting uses the minimum number of occurrences of the sequence, the default is 8, here is set to 10, the recommended minimum is one million of the total data volume, can remove low abundance noise and increase the calculation speed.
# -sizeout Output abundance, - -relabel with an added sequence prefix.

# OTU and ASV are collectively known as features (Feature), which:
# OTU usually selects the representative sequences with the highest abundance or center after 97% clustering;
# ASV is denoised based on sequences (excluding or correcting incorrect sequences and selecting credible sequences with higher abundance) as representative sequences

# There are two methods: recommended unoise3 denoising to obtain single base accuracy ASV, and the traditional 97% cluster OTU (horizontal accuracy) for alternative
# usearch Both feature selection methods have their own de novo desmographs

# Method ASV denoising Denoise: predict biological sequences and filter chimeras
#59s, 2920 good, 227 chimeras
usearch -unoise3 temp/uniques.fa \
  -zotus temp/zotus.fa
# Modify the sequence name: Zotu is changed to ASV for easy identification
sed 's/Zotu/ASV_/g' temp/zotus.fa > temp/otus.fa
head -n 2 temp/otus.fa
# To remove the singleton sequence
# vsearch --sortbysize otus1.fa --output otus.fa --minsize 2

### 4.3 reference-based deschimera Reference-based chimera detect

# Not recommended, which can easily cause false negatives because the reference database has no abundance information,
# However, de novo dechimerism requires more than 16 times the parental abundance of the chimeras to prevent false negatives
# Because known sequences will not be removed, larger database selection is more reasonable with the lowest false negative rate
mkdir -p result/raw

# Method 1. vsearch + rdp (fast but easily false negatives), or
# Silva removal (silva_16s_v123.fa), recommended (slow, 15m ~ 3h, but better)
vsearch --uchime_ref temp/otus.fa \
  --db ${ea}/usearch/rdp_16s_v16_sp.fa \
  --nonchimeras result/raw/otus.fa
# RDP: 51s, 250 (8.6%) chimeras; SILVA：10m, 255 (8.7%) chimeras

# Win vsearch Results added windows line changer ^ M to delete, mac do not execute this command
sed -i 's/\r//g' result/raw/otus.fa

# Method 2., no chimerism
cp -f temp/otus.fa result/raw/otus.fa
less result/raw/otus.fa|grep '>' -c

### 5. Features table generation and screening of Feature table create & filter
# Sequences with contamination rate higher than 0.99, set on demand
vsearch --usearch_global result/raw/otus.fa \
  --db temp/Wall.filtered.fa \
  --id 0.99 --query_cov 0.99 \
  --strand both --biomout result/raw/otus.txt --alnout result/raw/otus.aln \
  --blast6out result/raw/otus_blast.xls --fastapairs result/raw/otus_pairs.fa \
  --notmatched result/raw/otus_notmatched.fa --userfields query+target+id+qcov+tcov \
  --userout result/raw/otus_stat.xls
less result/raw/otus.fa|grep '>' -c
less result/raw/otus_notmatched.fa|grep '>' -c
less temp/filtered.fa|grep '>' -c
# 5.1 generated the feature Table Creat Feature table

# Method vsearch generated the feature table
vsearch --usearch_global temp/filtered.fa --db result/raw/otus_notmatched.fa \
  --otutabout result/raw/otutab.txt --id 0.97 --threads 12

# 660953 of 761432 (86.80%) comparable at 8m
# windows User removes the line character ^ M
sed -i 's/\r//' result/raw/otutab.txt
head -n3 result/raw/otutab.txt |cat -A
wc -l result/raw/otutab.txt
less result/raw/otus.fa|grep '>' -c
less result/raw/otus_notmatched.fa|grep '>' -c

### 5.2 Species annotation-Remove plastid and non-bacterial / archaea and statistical proportion (optional) Remove plastid and non-Bacteria

# RDP species annotation (rdp _ 16s _ v 16 _ sp) database is small, extremely fast, but lacks complete eukaryotic source data, taking 15s;
# SILVA database (silva_16s_v123.fa) for better annotated eukaryotic and plastid sequences, 3h
# --sintax_cutoff sets the confidence threshold for classification, which is usually 0.6 / 0.8. The larger, the lower the proportion of annotation.

vsearch --sintax result/raw/otus_notmatched.fa --db ${ea}/usearch/silva_16s_v123.fa.gz \ 	--tabbedout result/raw/otus.sintax --sintax_cutoff 0.6
vsearch --sintax result/raw/otus_notmatched.fa --db ${ea}/usearch/ltp_16s_v123.fa.gz \
  --tabbedout result/raw/ltp/otus.sintax --sintax_cutoff 0.6
vsearch --sintax result/raw/otus_notmatched.fa --db ${ea}/usearch/rdp_16s_v18_sp.fa.gz \
  --tabbedout result/raw/rdp/otus.sintax --sintax_cutoff 0.6
less result/raw/otus_notmatched.fa|grep '>' -c

# Number of original feature table rows
# wc -l result/raw/otutab.txt

# To remove non-specific amplification and assign contamination from 16S rDNA sequencing,
# Bacteria and archaea (prokaryotes), remove chloroplasts and mitochondria and statistical proportion, output OTU tables of filtered and sorted by abundance.
# Enter for the feature table (result / raw / otutab. The txt) and species annotation (result / raw / otus.sintax)，
# Output the filtered and sorted feature table (result / otutab. And txt), statistical pollution proportion file (result / raw / otutab_nonBac. The txt) and filter details (otus.sintax.discard)。
# Note: For fungal ITS data, use otutab _ filter _ nonFungi instead. R script, screening only for sequences annotated as fungi.
# For viewing the script help, please run Rscript $ {ea} / script / otutab_filter_nonBac.R -h
Rscript ${ea}/script/otutab_filter_nonBac.R \
  --input result/raw/otutab.txt \
  --taxonomy result/raw/otus.sintax \
  --output result/otutab.txt\
  --stat result/raw/otutab_nonBac.stat \
  --discard result/raw/otus.sintax.discard
# Number of feature table rows after screening, 2921-2904, deleted 17 ASV
wc -l result/otutab.txt

# Filter the sequences corresponding to the feature table
cut -f 1 result/otutab.txt | tail -n+2 > result/otutab.id
usearch -fastx_getseqs result/raw/otus.fa -labels result/otutab.id -fastaout result/otus.fa

# Filter feature table corresponding sequence annotation
awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}'\
result/raw/otus.sintax result/otutab.id \
  > result/otus.sintax

# Fill in the end column
sed -i 's/\t$/\td:Unassigned/' result/otus.sintax
head -n2 result/otus.sintax
less result/otus.fa|grep '>' -c
wc -l result/otutab.txt
wc -l result/otutab.id
# wc -l result/otutab_rare.txt
# wc -l result/taxonomy_rare.txt
# cp result/raw/otu* result/

# Optional statistical method: OTU table simple statistics Summary OTUs table
usearch -otutab_stats result/otutab.txt \
  -output result/otutab.stat
cat result/otutab.stat
# Note the minimum, quantile, or view result / raw / otutab_nonBac. Sample detailed data volume in stat, used for resampling

### 5.3 Equal sampling is standardized for the normlize by subsample

# Equal resampling using the vegan package, enter the reads count format Feature Table result / otutab.txt
# You can specify input file, sampling quantity and random number, and output draw table result / otutab _ are. The txt and diversity alpha / vegantxt
mkdir -p result/alpha/
Rscript ${ea}/script/otutab_rare.R --input result/otutab.txt \
  --depth 45310 --seed 123 \
  --normalize result/otutab_rare.txt \
  --output result/alpha/vegan.txt
usearch -otutab_stats result/otutab_rare.txt \
  -output result/otutab_rare.stat
cat result/otutab_rare.stat

cut -f 1-2 metadata.txt > temp/group.txt

### 6. Alpha diversity Alpha diversity

### 6.1. To calculate the diversity index, Calculate alpha diversity index
# Calculate all alpha diversity index (Chao1 without errors)
# details in http://www.drive5.com/usearch/manual/alpha_metrics.html
usearch -alpha_div result/otutab.txt \
  -output result/alpha/alpha.txt

### 6.2. To calculate the richness change of the dilution process, Rarefaction
# Dilution curve: Take the number of OTUs in 1% -100% of the sequences, 20s
# Rarefaction from 1%, 2% ..100% in richness (observed OTUs)-method fast / with_replacement / without_replacement https://drive5.com/usearch/manual/cmd_otutab_subsample.html
usearch -alpha_div_rare result/otutab_rare.txt \
  -output result/alpha/alpha_rare.txt \
  -method without_replacement

### 6.3. screen each group for comparison
# Calculate the mean of each feature, and then calculate the group mean, metadata according to the experimental design. The txt modified the group column names
# Input file is feautre table result / otutab. And txt, the experimental design of metadata.txt
# Output as feature table by group-an experiment may be grouped in multiple ways
Rscript ${ea}/script/otu_mean.R --input result/otutab_rare.txt \
  --design metadata.txt \
  --group Group --thre 1 \
  --output result/otutab_mean.txt
head -n3 result/otutab_mean.txt

# If the average abundance frequency is higher than one thousand (0.1%) is used as the screening criterion, five thousand or five thousand can be selected to obtain the OTU combinations of each group
awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i;} \
else {for(i=2;i<=NF;i++) if($i>0.1) print $1, a[i];}}' \
result/otutab_mean.txt > result/alpha/otu_group_exist.txt
head result/alpha/otu_group_exist.txt
# Results can be obtained directly at the http: / / www.ehbio. And com / ImageGP plot Venn, upSetView and Sanky
# Available at http: / / ehbio. In com / test / venn / plot and display common and unique elements to each group

### 7. Beta diversity Beta diversity

# Results with multiple files, required directory
mkdir -p result/beta/
# Build the evolutionary tree Make OTU tree based on the OTU, 30s
usearch -cluster_agg result/otus.fa -treeout result/otus.tree
# Generate 5 distance matrices: bray _ curtis, euclidean, jaccard, manhatten, unifrac, 3s
usearch -beta_div result/otutab_rare.txt -tree result/otus.tree \
  -filename_prefix result/beta/ # 1s
usearch -beta_div result/otutab_rare.txt -filename_prefix result/beta/

### 8. Species annotation and classification summary Taxonomy summary

# OTU corresponding species annotation 2 column format: remove the confidence value in sintax, and only keep the species annotation, replace: _, delete the quotation marks
cut -f 1,4 result/otus.sintax \
  |sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' \
  > result/taxonomy2.txt
head -n3 result/taxonomy2.txt

# OTU corresponds to the species 8-column format: Note that the annotation is non-neat
# The blank filling in OTU / ASV is Unassigned
awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
  split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
  print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
  result/taxonomy2.txt > temp/otus.tax
sed 's/;/\t/g;s/.__//g;' temp/otus.tax|cut -f 1-8 | \
  sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
  > result/taxonomy.txt
head -n3 result/taxonomy.txt

# Statistical taxonomic family genera, using the rank parameter p c o f g, for phylum, class, order, family, genome abbreviation
mkdir -p result/tax
for i in p c o f g;do
usearch -sintax_summary result/otus.sintax \
  -otutabin result/otutab_rare.txt -rank ${i} \
  -output result/tax/sum_${i}.txt
done
sed -i 's/(//g;s/)//g;s/\"//g;s/\#//g;s/\/Chloroplast//g' result/tax/sum_*.txt
# List all of the files
ls -sh result/tax/sum_*.txt
head -n3 result/tax/sum_g.txt

### 9. Space cleaning and data submission

rm -rf temp/*.fq

# Short-term without counting library compression to save space
gzip ${ea}/usearch/*.fa
gzip ${ea}/gg/*.fasta
gzip seq/*.fq
gzip seq/*.fastq
# Raw sequence statistics for md 5 values, used for data submission
cd seq
md5sum *_1.fq.gz> md5sum1.txt
md5sum *_2.fq.gz > md5sum2.txt
paste md5sum1.txt md5sum2.txt | awk '{print $2"\t"$1"\t"$4"\t"$3}' | \
sed 's/*//g' > ../result/md5sum.txt
rm md5sum*
cd ..
cat result/md5sum.txt

### R language diversity and species analysis

### 1.Alpha diversity

### 1.1 Alpha diversity
# View help (ANOVA + LSD calculation for normal distribution and homogeneity of variance)
Rscript ${ea}/script/alpha_boxplot.R -h
# Full parameters, diversity index optional richness chao1 ACE shannon simpson invsimpson
Rscript ${ea}/script/alpha_boxplot.R --alpha_index richness \
  --input result/alpha/vegan.txt --design metadata.txt \
  --group Group --output result/alpha/ \
  --width 89 --height 59
  
# Draw using a loop the 6 common indices
mkdir -p result/alpha/
for i in `head -n1 result/alpha/vegan.txt|cut -f 2-`;do
  Rscript ${ea}/script/alpha_boxplot.R --alpha_index ${i} \
    --input result/alpha/vegan.txt --design metadata.txt \
    --group Group --output result/alpha/ \
    --width 139 --height 109;done
for i in `head -n1 result/alpha/alpha.txt|cut -f 2-`;do
  Rscript ${ea}/script/alpha_boxplot.R --alpha_index ${i} \
    --input result/alpha/alpha.txt --design metadata.txt \
    --group Group --output result/alpha/ \
    --width 139 --height 109;done

# Alpha diversity bar graph + standard deviation
for i in `head -n1 result/alpha/alpha.txt|cut -f 2-`;do
  Rscript ${ea}/script/alpha_barplot.R --alpha_index ${i} \
    --input result/alpha/alpha.txt --design metadata.txt \
    --group Group --output result/alpha/ \
    --width 89 --height 59;done
for i in `head -n1 result/alpha/vegan.txt|cut -f 2-`;do
  Rscript ${ea}/script/alpha_barplot.R --alpha_index ${i} \
    --input result/alpha/vegan.txt --design metadata.txt \
    --group Group --output result/alpha/ \
    --width 89 --height 59;done

### 1.2 dilution curve
Rscript ${ea}/script/alpha_rare_curve.R \
  --input result/alpha/alpha_rare.txt --design metadata.txt \
  --group Group --output result/alpha/ \
  --width 139 --height 109

### 1.3 Diversity Venn diagram
# Comparison of three groups: -f input file, -a / b / c / d / g group name, -w / u is width and height inch, -p output file name suffix
bash ${ea}/script/sp_vennDiagram.sh \
  -f result/alpha/otu_group_exist.txt \
  -a WT -b KO -c OE \
  -w 3 -u 3 \
  -p WT_KO_OE
# For the four comparisons, see the input file directory for figures and code, and the running directory is the current project root directory
bash ${ea}/script/sp_vennDiagram.sh \
  -f result/alpha/otu_group_exist.txt \
  -a WT -b KO -c OE -d All \
  -w 3 -u 3 \
  -p WT_KO_OE_All

### 2.Beta diversity

### 2.1 Distance matrix heatmap pheatmap
# Add grouped annotations such as genotypes and sites in columns 2,4
cut -f 1-2 metadata.txt > temp/group.txt
# Take bray _ curtis as an example, -f input file, whether does-h cluster TRUE / FALSE, -u / v is the width-height inch
for i in weighted_unifrac unifrac bray_curtis jaccard;do
  bash ${ea}/script/sp_pheatmap.sh \
    -f result/beta/${i}.txt \
    -H 'TRUE' -u 5 -v 5;done

### 2.2 Principal Coordinate Analysis for PCoA
# Input file, select group, output file, picture size mm, statistics see beta _ pcoa _ stat.txt
for i in weighted_unifrac unifrac bray_curtis;do
  Rscript ${ea}/script/beta_pcoa.R \
    --input result/beta/${i}.txt --design metadata.txt \
    --group Group --output result/beta/${i}.pcoa.pdf \
    --width 139 --height 109;done

### 2.3 Restricted principal coordinate analysis of CPCoA
for i in wweighted_unifrac unifrac bray_curtis;do
  Rscript ${ea}/script/beta_cpcoa.R \
    --input result/beta/${i}.txt --design metadata.txt \
    --group Group --output result/beta/${i}.cpcoa.pdf \
    --width 139 --height 109;done

### 3. Species composition of the Taxonomy

### 3.1 Stacked bar graph Stackplot
# Taking the phylum (p) level as an example, the results include the output.sample/group. Two files for pdf
Rscript ${ea}/script/tax_stackplot.R \
  --input result/tax/sum_p.txt --design metadata.txt \
  --group Group --output result/tax/sum_p.stackplot \
  --legend 5 --width 89 --height 59
  
# The batch drawing input includes 5 levels of p / c / o / f / g
for i in k p c o f g s; do
  Rscript ${ea}/script/tax_stackplot.R \
    --input result/tax/sum_${i}.txt --design metadata.txt \
    --group Group --output result/tax/sum_${i}.stackplot \
    --legend 838 --width 909 --height 209; done

### 4. STAMP input file preparation

### 4.1 The command line generates the input file
Rscript ${ea}/script/format2stamp.R -h
mkdir -p result/stamp
Rscript ${ea}/script/format2stamp.R --input result/otutab_rare.txt \
  --taxonomy result/taxonomy.txt --threshold 0 \
  --output result/stamp/tax

### 5. LEfSe input file preparation

### 5.1. Command-line generation file

# Optional command line generates the input file
Rscript ${ea}/script/format2lefse.R -h
mkdir -p result/lefse
Rscript ${ea}/script/format2lefse.R --input result/otutab_rare.txt \
  --taxonomy result/taxonomy.txt --design metadata.txt \
  --group Group --threshold 0 \
  --output result/lefse/LEfSe

### 5.2 LEfSe analysis
# Method LEfSe official website online use

### 6. the QIIME 2 analysis process
# Code is detailed in 25QIIME2/pipeline_qiime2.sh
# External imported feature tables and representative sequences (commonly used)

# Statistical characteristics table.
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file metadata.txt

# Statistics represent the sequence.
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

### Alpha and beta diversity analysis

# Construct evolutionary trees for diversity analysis for 53s
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

# Calculate the core diversity
# 13s, the sampling depth usually selects the minimum value, coming from the table.qzv##
qiime diversity core-metrics-phylogenetic \
	--i-phylogeny rooted-tree.qza \
	--i-table feature-frequency-filtered-table.qza \
	--p-sampling-depth 45310 \
	--m-metadata-file metadata.txt \
	--output-dir core-metrics-results

# Significance analysis and visualization between groups of Alpha diversity
# 7s, the optional alpha indices are false _ pd, shannon, observed_features, evenness
index=shannon
index=faith_pd
index=observed_features
index=evenness
qiime diversity alpha-group-significance \
	--i-alpha-diversity core-metrics-results/${index}_vector.qza \
	--m-metadata-file metadata.txt \
	--o-visualization core-metrics-results/${index}-group-significance.qzv

# Significance analysis and visualization between groups of Beta diversity
# The optional beta indices are unweighted_unifrac, bray _ curtis, weighted_unifrac, and jaccard
# 7s, the assignment of grouping is to reduce the computational amount, the permutation test is time consuming

distance=weighted_unifrac
distance=unweighted_unifrac
distance=bray_curtis
column=Group
qiime diversity beta-group-significance \
	--i-distance-matrix core-metrics-results/${distance}_distance_matrix.qza \
	--m-metadata-file metadata.txt \
	--m-metadata-column ${column} \
	--o-visualization core-metrics-results/${distance}-${column}-significance.qzv \
	--p-pairwise

# leading-out table
qiime tools export \
  --input-path core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --output-path weighted_unifrac_distance_matrix

