import collections
import matplotlib.pyplot as plt
import numpy as np
import HTSeq


def parse_kegg(amino_acid = 'lys'):
    with open('../../output/ref/kegg/' + amino_acid + '-tRNA.txt') as f:
        return [field.strip().split('_')[1] for field in f.read().split('+')]


def merge_genes(genes):
    merged_genes = collections.Counter()
    for row in genes:
        if row[1] != '_no_feature' and row[1] != '_unmapped':
            print(row, row[1])
            merged_genes[row[0]] += int(row[2])
    return merged_genes


def assign_tRNA_classes(count_file):
    genes = []
    with open(count_file, 'r') as f:
        for line in f.readlines():
            fields = line.strip().split(',')
            count = fields[1]
            tRNA_class = 'unk'
            try:
                gene_id = fields[0].split(':')[1]
                for amino_acid, test_gene_ids in tRNA_class_dict.items():
                    if gene_id in test_gene_ids:
                        tRNA_class = amino_acid
            except IndexError:
                gene_id = fields[0]
            genes.append([tRNA_class, gene_id, count])

    print(sorted(genes, reverse=True))
    with open(count_file.split('.')[0] + '_assigned.csv', 'w') as f:
        for row in sorted(genes):
            f.write(','.join(row) + '\n')

    merged_genes = merge_genes(genes)
    with open(count_file.split('.')[0] + '_merged-tRNA-counts.csv', 'w') as f:
        for key, val in merged_genes.items():
            f.write(f'{key},{val}\n')
    return merged_genes


def plot_counts(plot_data):
    print('plotting...')
    labels = amino_acids + ['unk']
    counts = []
    for i, (count_file, merged_genes) in enumerate(plot_data):
        counts.append([])
        for amino_acid in labels:
            counts[i].append(merged_genes[amino_acid])

    x = np.arange(len(labels))  # the label locations
    margin = 0.15
    width = (0.5 - margin) * 2 / len(plot_data)

    fig, ax = plt.subplots()
    rects = []
    for i, (count_file, merged_genes) in enumerate(plot_data):
        rects.append(ax.bar(x - 0.5 + margin + width*i, counts[i], width, label=count_file))

    ax.set_ylabel('% of tRNA Counts')
    ax.set_title('tRNA counts by class')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    ax.set_yscale('log')

    def autolabel(rects):
        for rect in rects:
            height = rect.get_height()
            ax.annotate('{}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')
    for rect in rects:
        autolabel(rect)

    fig.tight_layout()

    plt.show()


def normalise_counts(plot_data):
    print('plot_data:')
    print(plot_data)
    normalised_plot_data = []
    for i, (count_file, merged_genes) in enumerate(plot_data):
        normalised_plot_data.append([count_file, collections.Counter()])
        print(count_file, merged_genes)
        tot = sum(merged_genes.values())
        for aa in merged_genes.keys():
            normalised_plot_data[i][1][aa] = round(merged_genes[aa] / tot * 100, 3)
    print('normalised_plot_data: ')
    print(normalised_plot_data)
    return normalised_plot_data

#########
annotations = HTSeq.GFF_Reader('../../output/ref/transcriptome/PRJNA13758/c_elegans.PRJNA13758.WS270.annotations.gff3')
tRNA_features = HTSeq.GenomicArrayOfSets('auto', stranded=True)
for feature in annotations:
    if feature.type == 'tRNA':
        tRNA_features[feature.iv] += feature.name
        print(feature, f'\niv: {feature.iv}\nname: {feature.name}\nscore: {feature.score}\nsource: {feature.source}\n'
                   f'type: {feature.type}\nattr: {feature.attr}\n\n\n')
print('tRNA list constructed')
#####
sam_paths = ['../../output/bwa/auckland/auckland_small-RNA_5-indep_1.sam',
             '../../output/bwa/bristol/bristol_small-rna_5-indep_1.sam']
#####
for sam_path in sam_paths:
    try:
        output_filename = sam_path.split('/')[len(sam_path.split('/')) - 1].split('.')[0] + '.csv'
    except IndexError:
        output_filename = sam_path + '.csv'

    counts = collections.Counter()
    almnt_file = HTSeq.SAM_Reader(sam_path)
    for almnt in almnt_file:
        if not almnt.aligned:
            counts['_unmapped'] += 1
            continue
        aligned_tRNA_IDs = set()
        for iv, val in tRNA_features[almnt.iv].steps():
            aligned_tRNA_IDs |= val  # constructs a set of all tRNAs IDs which the alignment could map to
        if len(aligned_tRNA_IDs) == 1:
            tRNA_ID = list(aligned_tRNA_IDs)[0]
            counts[tRNA_ID] += 1
        elif len(aligned_tRNA_IDs) == 0:
            counts['_no_feature'] += 1
        else:
            counts['_ambiguous'] += 1

    print(counts)

    with open(output_filename, 'w') as f:
        for key, val in counts.items():
            f.write(f'{key},{val}\n')
######

amino_acids = ['ala', 'asn', 'cys', 'glu', 'his', 'leu', 'met', 'pro', 'thr', 'tyr', 'arg', 'asp', 'gln', 'gly', 'ile',
               'lys', 'phe', 'ser', 'trp', 'val']
tRNA_class_dict = dict()
total_tRNAs = []
for amino_acid in amino_acids:
    for gene_id in parse_kegg(amino_acid):
        if gene_id in total_tRNAs:
            print(f'warning: duplicate found in {amino_acid}')
    total_tRNAs.append(parse_kegg(amino_acid))
    tRNA_class_dict[amino_acid] = parse_kegg(amino_acid)
print(tRNA_class_dict)
######
count_files = ['auckland_small-RNA_5-indep_1.csv',
               'bristol_small-rna_5-indep_1.csv']

plot_data = []
for count_file in count_files:
    merged_genes = assign_tRNA_classes(count_file)
    plot_data.append([count_file, merged_genes])
normalised_plot_data = normalise_counts(plot_data)
plot_counts(normalised_plot_data)
