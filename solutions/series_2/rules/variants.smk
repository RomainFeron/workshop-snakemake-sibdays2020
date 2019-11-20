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
        '../envs/variants.yaml'
    shell:
        'bcftools mpileup -f {input.genome} {input.aln} | bcftools call -P {params.subrate} -mv - > {output} 2>{log}'

rule parse_bcftools:
    input:
        rules.bcftools.output
    output:
        'results/variants.tsv'
    script:
        '../scripts/create_substitution_table.py'
