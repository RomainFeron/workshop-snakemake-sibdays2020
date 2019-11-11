# snakemake header
input_file = open(snakemake.input[0], 'r')
output_file = open(snakemake.output[0], 'w')

d={}

for line in input_file:
    if not line.startswith("#"):
        li=line.strip().split('\t')
        ref = li[3]
        alt = li[4]
        if ref != alt and len(ref) == 1 and len(alt) == 1:
            if ref not in d:
                d[ref] = {}
                d[ref][alt] = 0
            elif alt not in d[ref]:
                d[ref][alt] = 0
            d[ref][alt] += 1


with output_file as f:
    for k, v in d.items():
        f.write(','.join([str(x) for x in v.values()]) + '\n')
