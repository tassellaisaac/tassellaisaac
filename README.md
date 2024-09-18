This project implements a Python-based tool to analyze VCF (Variant Call Format) files, focusing on trio samples (parents and child) for identifying de novo mutations and annotating variants using external tools like ANNOVAR. The script is designed to parse the VCF file, summarize the variants, identify de novo mutations in trios, and annotate the variants with relevant biological information.

Before running the script, make sure you have the following packages installed:

pysam: A Python module for reading, manipulating, and writing SAM, BAM, CRAM, and VCF files.
pandas: A Python library for data manipulation and analysis.
subprocess: A built-in Python library to run external commands (used to call ANNOVAR).

pip install pysam pandas
You will also need ANNOVAR to annotate the variants.

Parameters:
vcf_file: The path to the VCF file (in this case, trios.vcf).
trio_samples: The sample names of the trio (parents and child) as listed in the VCF file.
annovar_path: The path to the ANNOVAR tool directory.
genome_build: The genome build (e.g., hg19, hg38) that corresponds to the reference data for annotation.


After running the script, the following output files will be generated:

variants_summary.csv: A CSV file summarizing the variants found in the VCF file.
de_novo_mutations.csv: A CSV file containing variants identified as de novo mutations.
annotated_variants.csv: A CSV file with the variants annotated using ANNOVAR, including gene annotations and functional effects.

The script uses pysam to parse the VCF file. It extracts useful information such as:

Chromosome (CHROM)
Position (POS)
Reference allele (REF)
Alternate allele(s) (ALT)
Quality score (QUAL)
Variant filters (FILTER)
Info field (INFO)
The parsed data is stored in a pandas DataFrame for easier manipulation and further analysis.


The script identifies de novo mutations in trios (parent1, parent2, and child) by comparing their genotypes. A de novo mutation is defined when:

Both parents are homozygous for the reference allele (0/0).
The child has a variant allele (0/1 or 1/1).
This helps identify new mutations that occurred in the child but were not inherited from the parents.

The script integrates with ANNOVAR to annotate the variants. ANNOVAR provides useful biological information such as:

The gene affected by the variant.
The functional impact (e.g., missense, nonsense, synonymous mutation).
The region (e.g., exonic, intronic) where the variant is located.
Running ANNOVAR:
The script generates a VCF-like file with the variants and calls ANNOVAR using the table_annovar.pl script. It annotates the variants using reference databases like refGene.

You can adjust the command and databases used in the annotate_variants() function of the script.

import pysam
import pandas as pd
import subprocess
import os

# Function to parse the VCF file and summarize variants
def summarize_variants(vcf_file):
    # Open the VCF file using pysam
    vcf = pysam.VariantFile(vcf_file)

    # Create lists to store data
    variants = []

    for record in vcf:
        variant_info = {
            "CHROM": record.chrom,
            "POS": record.pos,
            "REF": record.ref,
            "ALT": str(record.alts),
            "QUAL": record.qual,
            "FILTER": record.filter.keys(),
            "INFO": record.info.keys()
        }
        variants.append(variant_info)

    # Convert the list of variants to a DataFrame for easier analysis
    variants_df = pd.DataFrame(variants)
    
    # Summarize data
    summary = {
        "Total Variants": len(variants_df),
        "Chromosome Counts": variants_df["CHROM"].value_counts(),
        "Variant Types": variants_df["ALT"].value_counts(),
    }

    return summary, variants_df

# Function to identify de novo mutations in trios
def identify_de_novo_mutations(variants_df, trio_samples):
    # Assuming trio_samples is a list of the trio sample names in the VCF (parent1, parent2, child)
    parent1, parent2, child = trio_samples

    de_novo_mutations = []

    for idx, row in variants_df.iterrows():
        # Check the genotype of the parents and the child
        parent1_gt = row["FORMAT"][parent1]["GT"]
        parent2_gt = row["FORMAT"][parent2]["GT"]
        child_gt = row["FORMAT"][child]["GT"]

        # De novo mutation: Both parents are homozygous reference, but the child has the variant
        if parent1_gt == "0/0" and parent2_gt == "0/0" and child_gt != "0/0":
            de_novo_mutations.append(row)

    # Convert de novo mutations to DataFrame for further analysis
    de_novo_df = pd.DataFrame(de_novo_mutations)

    return de_novo_df

# Function to annotate variants using ANNOVAR
def annotate_variants(variants_df, annovar_path, genome_build):
    # Convert the variants DataFrame to a VCF-like format that ANNOVAR understands
    variants_df.to_csv("variants_for_annovar.vcf", sep="\t", index=False, header=False)

    # Run ANNOVAR to annotate the variants
    command = f"{annovar_path}/table_annovar.pl variants_for_annovar.vcf {annovar_path}/humandb/ -buildver {genome_build} -out annotated_variants -remove -protocol refGene -operation g -nastring ."
    subprocess.run(command, shell=True)

    # Load the annotated variants
    annotated_variants = pd.read_csv("annotated_variants.hg19_multianno.txt", sep="\t")

    return annotated_variants

# Main function to run the analysis
def main(vcf_file, trio_samples, annovar_path, genome_build):
    # Step 1: Summarize the variants
    summary, variants_df = summarize_variants(vcf_file)
    print("Variant Summary:", summary)

    # Step 2: Identify de novo mutations
    de_novo_df = identify_de_novo_mutations(variants_df, trio_samples)
    print("De Novo Mutations Identified:", len(de_novo_df))

    # Step 3: Annotate variants using ANNOVAR
    annotated_variants = annotate_variants(variants_df, annovar_path, genome_build)
    print("Annotated Variants:", annotated_variants.head())

    # Save results to CSV
    variants_df.to_csv("variants_summary.csv", index=False)
    de_novo_df.to_csv("de_novo_mutations.csv", index=False)
    annotated_variants.to_csv("annotated_variants.csv", index=False)


    # Path to the VCF file (be sure to copy it to your home directory)
    vcf_file = "/home/user/trios.vcf"
    
    # Sample names of the trio (parent1, parent2, child) as listed in the VCF
    trio_samples = ["parent1", "parent2", "child"]

    # Path to ANNOVAR and genome build version (e.g., hg19, hg38)
    annovar_path = "/path/to/annovar"
    genome_build = "hg19"

    # Run the main analysis
    main(vcf_file, trio_samples, annovar_path, genome_build)


