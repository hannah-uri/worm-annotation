import collections
import matplotlib.pyplot as plt
import numpy as np
import HTSeq


def parse_kegg(amino_acid='lys'):
    with open('../../output/ref/kegg/' + amino_acid + '-tRNA.txt') as f:
        return [field.strip().split('_')[1] for field in f.read().split('+')]


def merge_genes(genes):
    merged_genes = collections.Counter()
    for row in genes:
        if row[1] != '_no_feature' and row[1] != '_unmapped':
            print(row, row[1])
            merged_genes[row[0]] += int(row[2])
    return merged_genes


def assign_tRNA_classes(count_file, tRNA_class_dict):
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
    label_set = set()
    for count_file, trna_counts in plot_data:
        for amino_acid in trna_counts.keys():
            label_set.add(amino_acid)
    print('plotting...')
    labels = sorted(list(label_set))
    counts = []
    for i, (count_file, trna_counts) in enumerate(plot_data):
        counts.append([])
        for amino_acid in labels:
            counts[i].append(trna_counts[amino_acid])

    x = np.arange(len(labels))  # the label locations
    margin = 0.15
    width = (0.5 - margin) * 2 / len(plot_data)

    fig, ax = plt.subplots()
    rects = []
    for i, (count_file, merged_genes) in enumerate(plot_data):
        rects.append(ax.bar(x - 0.5 + margin + width * i, counts[i], width, label=count_file))

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
    for i, (count_file, trna_counts) in enumerate(plot_data):
        normalised_plot_data.append([count_file, collections.Counter()])
        print(count_file, trna_counts)
        tot = sum(trna_counts.values())
        for aa in trna_counts.keys():
            normalised_plot_data[i][1][aa] = round(trna_counts[aa] / tot * 100, 3)
    print('normalised_plot_data: ')
    print(normalised_plot_data)
    return normalised_plot_data


def parse_tRNAs_table(path):
    col_names = ['contig-name', 'tRNA #', 'begin', 'end', 'type', 'anti-codon', 'intron-begins', 'intron-ends',
                 'inf-score', 'note']
    tRNA_features = HTSeq.GenomicArrayOfSets('auto', stranded=False)
    with open(path) as f:
        for line in f.readlines()[3:]:
            fields = [field.strip() for field in line.split('\t')]
            if int(fields[col_names.index('begin')]) < int(fields[col_names.index('end')]):
                start = int(fields[col_names.index('begin')])
                end = int(fields[col_names.index('end')])
                strand = '+'
            else:
                start = int(fields[col_names.index('end')])
                end = int(fields[col_names.index('begin')])
                strand = '-'
            interval = HTSeq.GenomicInterval(chrom=fields[col_names.index('contig-name')],
                                             start=start,
                                             end=end,
                                             strand=strand)
            tRNA_features[interval] += fields[col_names.index('type')].lower()
    return tRNA_features


def parse_reference_annotation(annotation_path, tRNA_class_dict):
    print(tRNA_class_dict)
    inverted_tRNA_class_dict = {}
    for aa, gene_list in tRNA_class_dict.items():
        for gene in gene_list:
            inverted_tRNA_class_dict[gene] = aa
    print(inverted_tRNA_class_dict)
    annotations = HTSeq.GFF_Reader(annotation_path)
    tRNA_features = HTSeq.GenomicArrayOfSets('auto', stranded=True)
    for feature in annotations:
        if feature.type == 'tRNA':
            try:
                tRNA_features[feature.iv] += inverted_tRNA_class_dict[feature.name.split(':')[1]]
            except KeyError:
                tRNA_features[feature.iv] += 'unk'
            print(feature, f'\niv: {feature.iv}\nname: {feature.name}\nscore: {feature.score}\n'
                           f'source: {feature.source}\ntype: {feature.type}\nattr: {feature.attr}\n\n\n')
    print('tRNA list constructed')
    return tRNA_features


def generate_tRNA_class_dict():
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
    return tRNA_class_dict

#########
tRNA_table_paths = [
    '../../input/strains/auckland/tRNAs.txt',
    '../../input/strains/bristol/tRNAs.txt'
]
tRNA_features_list = []
for tRNA_table_path in tRNA_table_paths:
    tRNA_features_list.append(parse_tRNAs_table(tRNA_table_path))
print(tRNA_features_list)
tRNA_class_dict = generate_tRNA_class_dict()
reference_kegg_tRNA_features = parse_reference_annotation(
    '../../output/ref/transcriptome/PRJNA13758/c_elegans.PRJNA13758.WS270.annotations.gff3',
    tRNA_class_dict)

#####
# must be listed in same order as tRNA_table_paths, with stuff that
# uses the reference gff file at the end
sam_paths = [
    '../../output/bwa/auckland/auckland_small-rna_5-indep_1_contig-mapped.sam',
    '../../output/bwa/bristol/bristol_small-rna_5-indep_1_contig-mapped.sam',
    '../../output/bwa/auckland/auckland_small-rna_5-indep_1.sam',
    '../../output/bwa/bristol/bristol_small-rna_5-indep_1.sam'
]
#####
for i, sam_path in enumerate(sam_paths):
    if i < len(tRNA_features_list):
        tRNA_features = tRNA_features_list[i]
    else:
        tRNA_features = reference_kegg_tRNA_features

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

    print(i, sam_path, counts)

    with open(output_filename, 'w') as f:
        for key, val in counts.items():
            f.write(f'{key},{val}\n')
######

######
# count_files = ['auckland_small-rna_5-indep_1.csv',
#                'auckland_small-rna_5-indep_2.csv',
#                'bristol_small-rna_5-indep_1.csv',
#                'bristol_small-rna_5-indep_2.csv']
sam_filenames = [path.split('/')[len(path.split('/')) - 1] for path in sam_paths]
count_files = [filename.split('.')[0] + '.csv' for filename in sam_filenames]

plot_data = []
for i, count_file in enumerate(count_files):
    with open(count_file) as f:
        trna_counts = collections.Counter({line.split(',')[0]: int(line.split(',')[1].strip()) for line in f.readlines()
                                           if line.split(',')[0] != '_unmapped'
                                           and line.split(',')[0] != '_no_feature'})
    plot_data.append([sam_filenames[i], trna_counts])
normalised_plot_data = normalise_counts(plot_data)
plot_counts(normalised_plot_data)
