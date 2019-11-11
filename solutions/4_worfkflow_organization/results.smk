rule bcftools:
    input:
        genome = config['genome_in'],
        aln_idx = expand('results/{sample}.sorted.bam.bai', sample=config["samples"]),
        aln = expand('results/{sample}.sorted.bam', sample=config["samples"])
    output:
        'results/variants.vcf'
    params:
        subrate = config['subrate_bcftools']
    log:
        'logs/bcftools.log'
    benchmark:
        'benchmarks/bcftools.txt'
    conda:
        'raw.yaml'
    shell:
        'bcftools mpileup -f {input.genome} {input.aln} | bcftools call -P {params.subrate} -mv - > {output} 2>{log}'
