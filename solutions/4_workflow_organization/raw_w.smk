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
    wrapper:
        'v0.41.0/bio/samtools/view',
        'v0.41.0/bio/bwa/mem'
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
    wrapper:
        'v0.41.0/bio/samtools/sort'
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
    wrapper:
        'v0.41.0/bio/samtools/index'
    shell:
        'samtools index {input} 2>{log}'
