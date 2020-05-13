aed_arrays = {}
for filepath in snakemake.input:
    strain = filepath.split('/')[-1].split('.')[0]
    strain = strain.replace('-', ' ').title()
    aed_arrays[strain] = [0 for i in range(0, 101)]
    with open(filepath) as f:
        for line in f.readlines():
            try:
                if line.split()[2] == 'mRNA':
                    aed_pos = line.find('AED=')
                    if aed_pos == -1:
                        continue
                    aed = float(line[aed_pos + 4:].split(';')[0])
                    aed_arrays[strain][int(aed*100)] += 1
            except IndexError:
                continue

with open(snakemake.output[0], 'w') as f:
    f.write('Strain\tAED\tfreq\tcumulative_freq\tcumulative_frac\n')
    for strain, aed_array in aed_arrays.items():
        total = sum(aed_array)
        cumulative = 0
        for i, freq in enumerate(aed_array):
            cumulative += freq
            f.write('\t'.join([str(n) for n in [strain, round(i / 100, 2), freq, cumulative, round(cumulative / total, 5)]]) + '\n')

