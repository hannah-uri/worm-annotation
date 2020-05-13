import os

output_filename = snakemake.output[0]
with open(output_filename, 'w') as f:
    f.write('Strain, Replicate, Reads, Genes\n')

summary_files = []
read_files = []
for i, inp in enumerate(snakemake.input):
    if i < len(snakemake.input) / 2:
        read_files.append(inp)
    else:
        summary_files.append(inp)
for filepath in summary_files:
    filename = filepath.split('/')[-1]
    strain = filename.split('_')[0]
    try:
        replicate = filename.split('_')[1][1:].split('.')[0]
    except IndexError:
        print('could not identify replicate of ' + filename + ', skipping')
        continue
        
    reads = 0
    for n in replicate.split('+'):
        read_count_filename = 'output/cutadapt/' + strain + '/' + strain + '_r' + n + '.txt'
        try:
            with open(read_count_filename, 'r') as f:
                reads += int(f.readline().strip())
        except FileNotFoundError:
            print('could not find read record at ' + read_count_filename + ', skipping ' + filename) 
            continue
    reads = str(reads)

    with open(filepath, 'r') as f:
        print(filepath)
        for line in f.readlines()[1:]:
            if 'gene' in line.split()[1]:
                genes = line.split()[3].strip('\"').strip()
        
    strain = strain.capitalize()
    print(', '.join([strain, replicate, reads, genes]))
    with open(output_filename, 'a') as f:
        f.write(', '.join([strain, replicate, reads, genes]) + '\n')

