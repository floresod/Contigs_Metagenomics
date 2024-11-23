########################################################
#### Functional and Taxonomic Annotation of Contigs ####
########################################################

################################
#### Load Python Libraries ####
################################
import glob
import os

##########################
#### Global Variables ####
##########################
# Correctly capture all sample names from the Data directory
SAMPLE,=glob_wildcards("../resources/Data/{sample}.fasta")

########################
#### Global Outputs ####
########################
rule all:
    input:
        # Kraken2 outputs
        expand("../results/kraken2/outputs/{sample}.output", sample=SAMPLE),
        expand("../results/kraken2/reports/{sample}.report", sample=SAMPLE),
        # CARD alignment outputs
        expand("../results/card/{sample}.tsv", sample=SAMPLE)

##########################################
#### Taxonomic Classification Kraken2 ####
##########################################
rule kraken2:
    input: 
        contigs = "../resources/Data/{sample}.fasta"
    output:
        report = "../results/kraken2/reports/{sample}.report",
        kraken_output = "../results/kraken2/outputs/{sample}.output"
    params:
        database = "../../../Databases/k2_standard_08gb_20240605/"
    threads:
        6
    log:
        "../resources/Logs/kraken2/{sample}.log"
    conda:
        "../envs/kraken2_env.yaml"
    shell:
        """
        mkdir -p $(dirname {output.report}) $(dirname {output.kraken_output})
        kraken2 --db {params.database} \
                {input.contigs} \
                --report {output.report} \
                --output {output.kraken_output} \
                --threads {threads} \
                --confidence 0.005 \
                --use-names \
                > {log} 2>&1
        """

########################
#### CARD Alignment ####
########################
rule diamond_card:
    input:
        contigs = "../resources/Data/{sample}.fasta"
    output:
        card_output = "../results/card/{sample}.tsv"
    params:
        database = "../../../Databases/card/card_v3.3.0.dmnd"
    resources:
        threads = 6
    log:
        "../resources/Logs/card/{sample}.log"
    conda:
        "../envs/diamond_env.yaml"
    shell:
        """
        mkdir -p $(dirname {output.card_output})
        diamond blastx -d {params.database} \
                      -q {input.contigs} \
                      -o {output.card_output} \
                      --id 95 \
                      --subject-cover 90 \
                      > {log} 2>&1
        """

