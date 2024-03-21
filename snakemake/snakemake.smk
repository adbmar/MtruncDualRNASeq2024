#for running on a cluster with SLURM
"""
snakemake -s snakemake.smk\
 --use-conda\
 --latency-wait 30\
 --conda-frontend conda\
 --rerun-incomplete\
 --keep-going\
 --jobs 12\
 --default-resource mem_mb=12000\
 --cluster 'sbatch --cpus-per-task {threads} --mem {resources.mem_mb} --output ~/SLURM_Logs/snakemake/{rule}.{wildcards.sample}.o --error ~/SLURM_Logs/snakemake/{rule}.{wildcards.sample}.e --time 1-0:0:0'
"""

import itertools
import os
import glob


samples = list(range(1, 46, 4)) + list(range(2, 47, 4))
samples.sort()
samples = [str(sample) for sample in samples]

organisms = ["Medicago", "Nematode", "Rhizob21", "Rhizob22"]

count_targets = [str(target[0]) + "/" + str(target[1]) + "_count" for target in itertools.product([organisms[0]], samples)]
count_targets = count_targets + [str(target[0]) + "/" + str(target[1]) + "_count" for target in itertools.product([organisms[1]], [sample for sample in samples if sample in range(13,23)])]
count_targets.append([str(target[0]) + "/" + str(target[1]) + "_count" for target in itertools.product([organisms[2]], [sample for sample in samples if sample in range(25,46,4)])])
count_targets.append([str(target[0]) + "/" + str(target[1]) + "_count" for target in itertools.product([organisms[3]], [sample for sample in samples if sample in range(26,47,4)])])

count_targets = ["counts/" + "Medicago" + "/" + sample + "_count" for sample in samples]
count_targets.append(["counts/" + "Nematode" + "/" + sample + "_count" for sample in samples if int(sample) in range(13,23)])
count_targets.append(["counts/" + "Rhizob21" + "/" + sample + "_count" for sample in samples if int(sample) in range(25,46,4)])
count_targets.append(["counts/" + "Rhizob22" + "/" + sample + "_count" for sample in samples if int(sample) in range(26,47,4)])


print(count_targets)

rawinputs1 = ["rawReads/" + str(sample) + "_R1_001.fastq.gz" for sample in samples]
rawinputs2 = ["rawReads/" + str(sample) + "_R2_001.fastq.gz" for sample in samples]
rawinputs = rawinputs1 + rawinputs2

trimmed_targets1 = ["trimmedReads/" + str(sample) + "_1.fastq.gz" for sample in samples]
trimmed_targets2 = ["trimmedReads/" + str(sample) + "_2.fastq.gz" for sample in samples]
trimmed_targets = trimmed_targets1 + trimmed_targets2

sorted_targets1 = ["trimmedAndFilteredReads/" + str(sample) + "_fwd.fq.gz" for sample in samples]
sorted_targets2 = ["trimmedAndFilteredReads/" + str(sample) + "_rev.fq.gz" for sample in samples]
sorted_targets = sorted_targets1 + sorted_targets2

rule all:
    input:
        #targets = count_targets + ["rawReads/multiqc_report.html", "trimmedReads/multiqc_report.html"] + trimmed_targets + sorted_targets
        targets = count_targets + trimmed_targets + sorted_targets


#rule rawQC:
#Gets quality for raw input files to use as comparison after trimming
#    input:
#        rawinputs
#    output:
#        "rawReads/multiqc_report.html"
#    conda:
#        "env_qc.yml"
#    threads:
#        1
#    shell:
#        """
#        fastqc -o rawReads/ {input}
#        multiqc -o rawReads/ rawReads/
#        """

rule trimmomatic:
#Basic quality control on reads, removing adapters, and removing low quality reads
    input:
        read1 = "rawReads/{sample}_R1_001.fastq.gz",
        read2 = "rawReads/{sample}_R2_001.fastq.gz"
    output:
        paired1 = temp("trimmedReads/{sample}_1.fastq.gz"),
        paired2 = temp("trimmedReads/{sample}_2.fastq.gz"),
        sterrLog = "trimmedReads/sterr.{sample}.log",
        stoutLog = "trimmedReads/stout.{sample}.log"
    threads:
        8
    resources:
        mem_mb = 24000
    conda:
        "env_trimmomatic.yml"
    params:
        unpaired1 = "discardReads/{sample}_trimmomatic_unpaired_1.fastq.gz",
        unpaired2 = "discardReads/{sample}_trimmomatic_unpaired_2.fastq.gz",
        adapterdirectory = "adapter/"

    shell:
        """
        trimmomatic PE -threads 8 -phred33 {input.read1} {input.read2} {output.paired1} {params.unpaired1} {output.paired2} {params.unpaired2} ILLUMINACLIP:{params.adapterdirectory}TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:15 TRAILING:15 SLIDINGWINDOW:5:15 MINLEN:50 1> {output.sterrLog} 2> {output.stoutLog}
        rm {params.unpaired1}
        rm {params.unpaired2}
        """

#rule trimmedQC:
#Gets quality for trimmed reads to use as comparison to before trimming
#Quality control check
#    input:
#        trimmed_targets
#    output:
#        "trimmedReads/multiqc_report.html"
#    threads:
#        1
#    conda:
#        "env_qc.yml"
#    shell:
#        """
#        fastqc -o trimmedReads/ {input} 
#        multiqc -o trimmedReads/ trimmedReads/
#        """


rule sortmeRNA:
#Removes rRNA reads from the dataset
#rRNA reads put into discardReads; non-rRNA reads used downstream
    input:
        trimmed1 = "trimmedReads/{sample}_1.fastq.gz",
        trimmed2 = "trimmedReads/{sample}_2.fastq.gz",
        sterrLog = "trimmedReads/sterr.{sample}.log",
        stoutLog = "trimmedReads/stout.{sample}.log"
    output:
        sterrLog = "trimmedAndFilteredReads/sortmerna.sterr.{sample}.log",
        stoutLog = "trimmedAndFilteredReads/sortmerna.stout.{sample}.log",
        rRNAFiltered_fwd = temp("trimmedAndFilteredReads/{sample}_fwd.fq.gz"),
        rRNAFiltered_rev = temp("trimmedAndFilteredReads/{sample}_rev.fq.gz")
    params:
        work_dir = "~/sortmerna/{sample}/run",
        rRNAFiltered_dir = "trimmedAndFilteredReads/{sample}",
        rRNAReads_dir = "rRNAReads/{sample}",
        rRNAReferenceSequence = "AllrRNAReferences.fasta",
        rRNAReads_fwd = "rRNAReads/{sample}_fwd.fq.gz",
        rRNAReads_rev = "rRNAReads/{sample}_rev.fq.gz"
    threads:
        8
    resources:
        mem_mb = 48000
    conda:
        "env_sortmerna.yml"
    shell:
        """
        rm -rf {params.work_dir}
        sortmerna --ref {params.rRNAReferenceSequence} --reads {input.trimmed1} --reads {input.trimmed2} --print_all_reads True --fastx True --aligned {params.rRNAReads_dir} --other {params.rRNAFiltered_dir} --out2 True --workdir {params.work_dir} --paired_in --zip-out true -a 8 1> {output.sterrLog} 2> {output.stoutLog}
        rm {params.rRNAReads_fwd}
        rm {params.rRNAReads_rev}
        """


rule sortedQC:
#Gets quality for trimmed and sorted reads to use as comparison to other QC
#Quality control check
    input:
        sorted_targets
    output:
        "trimmedAndFilteredReads/multiqc_report.html"
    conda:
        "env_qc.yml"
    shell:
        """
        fastqc -o trimmedAndFilteredReads/ {input} 
        multiqc -o trimmedAndFilteredReads/ trimmedAndFilteredReads/
        """

rule alignMedicago:
#Alignment to Medicago genome for downstream gene expression quantification
    input:
        fwd_reads = "trimmedAndFilteredReads/{sample}_fwd.fq.gz",
        rev_reads = "trimmedAndFilteredReads/{sample}_rev.fq.gz"
    output:
        bam = "MedicagoAlignment/{sample}.bam",
        summary = "MedicagoAlignment/{sample}.log",
        sterrLog = "MedicagoAlignment/hisat2.sterr.{sample}.log",
        stoutLog = "MedicagoAlignment/hisat2.stout.{sample}.log"
    params:
        sam = "MedicagoAlignment/{sample}.sam",
        hisatindex = "MedicagoIndex"
    threads:
        8
    resources:
        mem_mb = 48000
    conda:
        "env_hisat.yml"
    shell:
        """
        hisat2 --min-intronlen 40 --max-intronlen 8000 --very-sensitive --threads 8 --summary-file {output.summary} -x {params.hisatindex} -S {params.sam} -1 {input.fwd_reads} -2 {input.rev_reads} 1> {output.sterrLog} 2> {output.stoutLog}
        samtools view {params.sam} -@ 8 -S -b -o {output.bam}
        rm {params.sam}
        """


rule filteringAlignmentsMedicago:
#Filtering:
# 1. low quality alignments (-q 30 in samtools view)
# 2. reads whose mate aligned to a different chromosme (awk command)
#Then sorting in preparation for HTSeq quantification
    input:
        bam = "MedicagoAlignment/{sample}.bam"
    output:
        bam = "MedicagoAlignment/{sample}.filtered.sorted.bam"
    threads:
        1
    conda:
        "env_samtools.yml"
    shell:
        """
        samtools view -h -@ 1 -q 30 {input.bam} | awk '{{ if ((substr($1,1,1) == /@/ ) || (substr($1,1,1) !~ /@/ && $7 !~ /=/ )) {{printf ""}} else {{print $0}} }}' | samtools sort - -O bam -o {output.bam} -@ 1
        """

rule countMedicago:
#Quantification of gene expression using HTSeq
    input:
        "MedicagoAlignment/{sample}.filtered.sorted.bam"
    output:
        "counts/Medicago/{sample}_count"
    params:
        gtf = "MedicagoAnnotation.gtf"
    threads:
        1
    conda:
        "env_htseq.yml"
    shell:
        """
        htseq-count -r pos -f bam -m union --stranded=no {input} {params.gtf} > {output}
        """

rule alignNematode:
#Alignment to nematode genome for downstream gene expression quantification
    input:
        fwd_reads = "trimmedAndFilteredReads/{sample}_fwd.fq.gz",
        rev_reads = "trimmedAndFilteredReads/{sample}_rev.fq.gz"
    output:
        bam = "NematodeAlignment/{sample}.bam",
        summary = "NematodeAlignment/{sample}.log",
        sterrLog = "NematodeAlignment/hisat2.sterr.{sample}.log",
        stoutLog = "NematodeAlignment/hisat2.stout.{sample}.log"
    params:
        sam = "NematodeAlignment/{sample}.sam",
        hisatindex = "NematodeIndex"
    threads:
        8
    resources:
        mem_mb = 48000
    conda:
        "env_hisat.yml"
    shell:
        """
        hisat2 --min-intronlen 40 --max-intronlen 8000 --very-sensitive --threads 8 --summary-file {output.summary} -x {params.hisatindex} -S {params.sam} -1 {input.fwd_reads} -2 {input.rev_reads} 1> {output.sterrLog} 2> {output.stoutLog}
        samtools view {params.sam} -@ 8 -S -b -o {output.bam}
        rm {params.sam}
        """

rule filteringAlignmentsNematode:
#Filtering:
# 1. low quality alignment (-q 30 in samtools view)
# 2. reads whose mate aligned to a different chromosme (awk command)
#Then sorting in preparation for HTSeq quantification
    input:
        bam = "NematodeAlignment/{sample}.bam"
    output:
        bam = "NematodeAlignment/{sample}.filtered.sorted.bam"
    threads:
        1
    conda:
        "env_samtools.yml"
    shell:
        """
        samtools view -h -@ 1 -q 30 {input.bam} | awk '{{ if ((substr($1,1,1) == /@/ ) || (substr($1,1,1) !~ /@/ && $7 !~ /=/ )) {{printf ""}} else {{print $0}} }}' | samtools sort - -O bam -o {output.bam} -@ 1
        """

rule countNematode:
#Quantification of gene expression using HTSeq
    input:
        "NematodeAlignment/{sample}.filtered.sorted.bam"
    output:
        "counts/Nematode/{sample}_count"
    params:
        gtf = "NematodeAnnotation.gtf"
    conda:
        "env_htseq.yml"
    threads:
        1
    shell:
        """
        htseq-count -r pos -f bam -m union --stranded=no {input} {params.gtf} > {output}
        """

rule alignRhizob21:
#Alignment to Em1021 genome for downstream gene expression quantification
    input:
        fwd_reads = "trimmedAndFilteredReads/{sample}_fwd.fq.gz",
        rev_reads = "trimmedAndFilteredReads/{sample}_rev.fq.gz"
    output:
        bam = "Rhizob21Alignment/{sample}.bam",
        summary = "Rhizob21Alignment/{sample}.log",
        sterrLog = "Rhizob21Alignment/bowtie2.sterr.{sample}.log",
        stoutLog = "Rhizob21Alignment/bowtie2.stout.{sample}.log"
    threads:
        8
    resources:
        mem_mb = 48000
    params:
        sam = "Rhizob21Alignment/{sample}.sam",
        bowtie2index = "Rhizob21Index"
    conda:
        "env_bowtie.yml"
    shell:
        """
        bowtie2 -q -p 8 --very-sensitive -x {params.bowtie2index} -1 {input.fwd_reads} -2 {input.rev_reads} --met-file {output.summary} -S {params.sam} 1> {output.sterrLog} 2> {output.stoutLog}
        samtools view {params.sam} -@ 8 -S -b -o {output.bam}
        rm {params.sam}
        """


rule filteringAlignmentsRhizob21:
#Filtering:
# 1. low quality alignment (-q 30 in samtools view)
# 2. reads whose mate aligned to a different chromosme (awk command)
#Then sorting in preparation for HTSeq quantification
    input:
        bam = "Rhizob21Alignment/{sample}.bam"
    output:
        bam = "Rhizob21Alignment/{sample}.filtered.sorted.bam"
    threads:
        1
    conda:
        "env_samtools.yml"
    shell:
        """
        samtools view -h -@ 1 -q 30 {input.bam} | awk '{{ if ((substr($1,1,1) == /@/ ) || (substr($1,1,1) !~ /@/ && $7 !~ /=/ )) {{printf ""}} else {{print $0}} }}' | samtools sort - -O bam -o {output.bam} -@ 1
        """

rule countRhizob21:
#Quantification of gene expression using HTSeq
    input:
        "Rhizob21Alignment/{sample}.filtered.sorted.bam"
    output:
        "counts/Rhizob21/{sample}_count"
    params:
        gtf = "Rhizob21Annotation.gtf"
    conda:
        "env_htseq.yml"
    shell:
        """
        htseq-count -r pos -f bam -m union --stranded=no --type=gene {input} {params.gtf} > {output}
        """

rule alignRhizob22:
#Alignment to Em1022 genome for downstream gene expression quantification
    input:
        fwd_reads = "trimmedAndFilteredReads/{sample}_fwd.fq.gz",
        rev_reads = "trimmedAndFilteredReads/{sample}_rev.fq.gz"
    output:
        bam = "Rhizob22Alignment/{sample}.bam",
        summary = "Rhizob22Alignment/{sample}.log",
        sterrLog = "Rhizob22Alignment/bowtie2.sterr.{sample}.log",
        stoutLog = "Rhizob22Alignment/bowtie2.stout.{sample}.log"
    params:
        sam = "Rhizob22Alignment/{sample}.sam",
        bowtie2index = "Rhizob22Index"
    threads:
        8
    resources:
        mem_mb = 48000
    conda:
        "env_bowtie.yml"
    shell:
        """
        bowtie2 -q -p 8 --very-sensitive -x {params.bowtie2index} -1 {input.fwd_reads} -2 {input.rev_reads} --met-file {output.summary} -S {params.sam} 1> {output.sterrLog} 2> {output.stoutLog}
        samtools view {params.sam} -@ 8 -S -b -o {output.bam}
        rm {params.sam}
        """


rule filteringAlignmentsRhizob22:
#Filtering:
# 1. low quality alignments (-q 30 in samtools view)
# 2. reads whose mate aligned to a different chromosme (awk command)
#Then sorting in preparation for HTSeq quantification
    input:
        bam = "Rhizob22Alignment/{sample}.bam"
    output:
        bam = "Rhizob22Alignment/{sample}.filtered.sorted.bam"
    threads:
        1
    conda:
        "env_samtools.yml"
    shell:
        """
        samtools view -h -@ 1 -q 30 {input.bam} | awk '{{ if ((substr($1,1,1) == /@/ ) || (substr($1,1,1) !~ /@/ && $7 !~ /=/ )) {{printf ""}} else {{print $0}} }}' | samtools sort - -O bam -o {output.bam} -@ 1
        """

rule countRhizob22:
#Quantification of gene expression using HTSeq
    input:
        "Rhizob22Alignment/{sample}.filtered.sorted.bam"
    output:
        "counts/Rhizob22/{sample}_count"
    params:
        gtf = "Rhizob22Annotation.gtf"
    conda:
        "env_htseq.yml"
    shell:
        """
        htseq-count -r pos -f bam -m union --type=gene --stranded=no {input} {params.gtf} > {output}
        """
