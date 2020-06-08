rule bcftools:
    input:
        ref = config['genome_in'],
        indexes = expand('results/{sample}.sorted.bam.bai', sample=config["samples"]),
        samples = expand('results/{sample}.sorted.bam', sample=config["samples"])
    output:
        'results/variants.vcf'
    params:
        call = f'-P {config["subrate_bcftools"]}'
    log:
        'logs/bcftools.log'
    benchmark:
        'benchmarks/bcftools.txt'
    wrapper:
        'v0.41.0/bio/bcftools/call'

rule parse_bcftools:
    input:
        rules.bcftools.output
    output:
        'results/variants.tsv'
    script:
        '../scripts/create_substitution_table.py'
