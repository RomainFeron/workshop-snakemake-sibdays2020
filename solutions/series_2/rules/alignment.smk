def get_sample_path(wildcards):
    return config['samples'][wildcards.sample]


rule bwa_map_gen:
    input:
        config['genome_in'],
        get_sample_path
    output:
        temp('results/{sample}.bam')
    threads:
        config['thread_bwa']
    log:
        'logs/{sample}_bwa_map.log'
    benchmark:
        'benchmarks/{sample}_bwa_map.txt'
    conda:
        '../envs/alignment.yaml'
    shell:
        'bwa mem -t {threads} {input} | samtools view -b > {output} 2>{log}'

rule sort_bam:
    input:
        rules.bwa_map_gen.output
    output:
        'results/{sample}.sorted.bam'
    log:
        'logs/{sample}_sort.log'
    benchmark:
        'benchmarks/{sample}_sort.txt'
    conda:
        '../envs/alignment.yaml'
    shell:
        'samtools sort -O bam {input} > {output} 2>{log}'

rule samtools_idx:
    input:
        rules.sort_bam.output
    output:
        'results/{sample}.sorted.bam.bai'
    log:
        'logs/{sample}_samtools_idx.log'
    benchmark:
        'benchmarks/{sample}_samtools_idx.txt'
    conda:
        '../envs/alignment.yaml'
    shell:
        'samtools index {input} 2>{log}'
