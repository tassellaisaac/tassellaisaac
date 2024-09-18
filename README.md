This practical aims to process and analyze ChIP-seq data from mouse hemangioblast cells, focusing on the Tal1 and Lmo2 transcription factors, which are essential in hematopoiesis. Through steps such as duplicate removal, peak calling, visualization, and annotation of ChIP-seq data to identify where Tal1 and Lmo2 proteins bind the chromatin.

Remove duplicates from aligned ChIP-seq data using Picard.
Call peaks using MACS2 to identify protein-DNA interactions.
Visualize peaks using IGV.
Filter peaks by removing blacklisted regions using BEDTools.
Annotate peaks and identify associated genes using HOMER.
Compare DNA binding patterns between Tal1 and Lmo2 proteins.

The following software modules need to be loaded on the BlueBear HPC environment:

Picard: For duplicate removal.
SAMtools: For BAM file inspection.
MACS2: For peak calling.
BEDTools: For filtering peaks.
IGV: For visualizing peaks.
HOMER: For peak annotation.

The data provided includes:

ChIP-seq aligned BAM files for Tal1 and Lmo2 in mouse hemangioblast cells.
The mouse mm10 genome for alignment.
Annotation files such as blacklist and simple repeats for filtering.

#Duplicate Removal in BASH
java -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT=HB_Lmo2_aligned_sorted.bam OUTPUT=HB_Lmo2_aligned_sorted_noDups.bam METRICS_FILE=HB_Lmo2_aligned_sorted_noDups.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true

#Peak Calling in BASH
macs2 callpeak -t HB_Lmo2_aligned_sorted_noDups.bam -c HB_input_aligned_sorted_noDups.bam -f BAM -g mm -n lmo2 -q 0.05 --keep-dup auto -B --trackline

#Visualising peaks
igvtools toTDF tal1_treat_pileup.bdg tal1_treat_pileup.tdf mm10
igvtools toTDF lmo2_treat_pileup.bdg lmo2_treat_pileup.tdf mm10

#Remove Blacklisted Regions
bedtools intersect -a lmo2_summits.bed -b Annotation_files/mm10_simpleRepeat.bed -v | bedtools intersect -a - -b Annotation_files/mm10-blacklist.bed -v > lmo2_noBL_summits.bed


#Peak Filtering
cut -f4 lmo2_noBL_summits.bed | grep -Fwf - lmo2_peaks.narrowPeak | grep -v chrUn | grep -v chrM | grep -v random > lmo2_filtered_peaks.bed

#Compare Binding Patterns
bedtools intersect -a lmo2_filtered_peaks.bed -b tal1_filtered_peaks.bed -v > lmo2_specific.bed
bedtools intersect -a lmo2_filtered_peaks.bed -b tal1_filtered_peaks.bed -u > shared_peaks.bed
bedtools intersect -a tal1_filtered_peaks.bed -b lmo2_filtered_peaks.bed -v > tal1_specific.bed

#Peak Annotation
annotatePeaks.pl lmo2_specific.bed mm10 -annStats lmo2_peaks_annstats.tsv > lmo2_peaks_annotation.tsv 2> lmo2_peaks_homer.out

#Extract Gene Lists
awk -F "\t" '(NR>1){print $16}' lmo2_peaks_annotation.tsv | sort | uniq > lmo2_specific_genes.txt

Output Files
Duplicate-free BAM files: HB_Lmo2_aligned_sorted_noDups.bam, HB_Tal1_aligned_sorted_noDups.bam.
Peak calling results: lmo2_peaks.narrowPeak, tal1_peaks.narrowPeak.
Filtered peaks: lmo2_filtered_peaks.bed, tal1_filtered_peaks.bed, shared_peaks.bed.
Gene lists: lmo2_specific_genes.txt, tal1_specific_genes.txt, shared_genes.txt.
