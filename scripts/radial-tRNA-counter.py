import collections
import matplotlib.pyplot as plt
import numpy as np
import HTSeq
from pprint import pprint as pp
import pickle
import os
##### snakemake stuff
paths = str(snakemake.input).split()
samples = [{'sam_path': '../../' + paths[i], 'sample': snakemake.params[0][i].capitalize()} for i in range(len(paths))]
output_path = '../../' + str(snakemake.output)  #output_filename = str(snakemake.output).split('/')[-1]
os.chdir('pycharm/worms')
#####

annotation_path = '../../output/ref/annotation/c_elegans.PRJNA13758.WS270.annotations.gff3'
# samples = [
#     {'sam_path': '../../output/bwa/bristol/bristol_small-rna_5-indep_1.sam',
#      'sample': 'Bristol'},
#     {'sam_path': '../../output/bwa/auckland/auckland_small-rna_5-indep_1.sam',
#      'sample': 'Auckland'}
# ]

def parse_reference_annotation(annotation_path, tRNA_class_dict):
    print('generating tRNA_features genomic array...')
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


def parse_kegg(amino_acid='lys'):
    with open('../../output/ref/kegg/' + amino_acid + '-tRNA.txt') as f:
        return [field.strip().split('_')[1] for field in f.read().split('+')]


def save_tRNA_features():
    tRNA_class_dict = generate_tRNA_class_dict()
    reference_kegg_tRNA_features = parse_reference_annotation(annotation_path, tRNA_class_dict)
    with open('tRNA_features.obj', 'wb') as f:
        pickle.dump(reference_kegg_tRNA_features, f)
    return reference_kegg_tRNA_features


def read_tRNA_features():
    with open('tRNA_features.obj', 'rb') as f:
        return pickle.load(f)


def count_tRNAs(samples, tRNA_features, output_path='tRNA-counts.tsv'):
    samples_with_counts = []
    with open(output_path, 'w') as f:
        f.write('sample\ttRNA\tcount\tfraction\n')
    for sample in samples:
        print(f'counting tRNAs in {sample["sample"]}')
        sample['counts'] = collections.Counter()
        tRNA_count_total = 0
        almnt_file = HTSeq.SAM_Reader(sample['sam_path'])
        for almnt in almnt_file:
            if not almnt.aligned:
                sample['counts']['_unmapped'] += 1
                continue
            aligned_tRNA_IDs = set()
            for iv, val in tRNA_features[almnt.iv].steps():
                aligned_tRNA_IDs |= val  # constructs a set of all tRNA IDs which the alignment could map to
            if len(aligned_tRNA_IDs) == 1:
                tRNA_ID = list(aligned_tRNA_IDs)[0]
                sample['counts'][tRNA_ID] += 1
                tRNA_count_total += 1
            elif len(aligned_tRNA_IDs) == 0:
                sample['counts']['_not_tRNA'] += 1
            else:
                sample['counts']['unk'] += 1  # ambiguous
                tRNA_count_total += 1

        sample['fractions'] = {}
        for tRNA, count in sorted(sample['counts'].items()):
            fraction = round(count / tRNA_count_total, 5)
            sample['fractions'][tRNA] = fraction
            if tRNA != '_unmapped' and tRNA != '_not_tRNA':
                with open(output_path, 'a') as f:
                    f.write(f'{sample["sample"]}\t{tRNA}\t{str(count)}\t{str(fraction)}\n')
        samples_with_counts.append(sample)
    return samples_with_counts


def read_samples_with_counts(filename='tRNA-counts.tsv'):
    samples_with_counts = []
    with open(filename, 'r') as f:
        for line in f.readlines()[1:]:
            fields = line.strip().split('\t')
            sample_name = fields[0]
            tRNA = fields[1]
            count = int(fields[2])
            fraction = float(fields[3])
            for i, sample in enumerate(samples_with_counts):
                if sample['sample'] == sample_name:
                    break
            else:
                samples_with_counts.append({'sample': sample_name, 'counts': {}, 'fractions': {}})
                i = len(samples_with_counts) - 1
            samples_with_counts[i]['counts'][tRNA] = count
            samples_with_counts[i]['fractions'][tRNA] = fraction
    return samples_with_counts


def plot_radial(samples_with_counts):
    tRNAs = list(sorted(samples_with_counts[0]['fractions'].keys()))
    angle_gap = 360 / len(tRNAs)
    theta = [angle_gap * i for i in range(0, len(tRNAs))]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='polar')
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_xticks([t/360 * 2 * 3.14159 for t in theta])
    ax.set_xticklabels(tRNAs)
    print(ax.get_yticks())

    for sample in samples_with_counts[0:1]:
        theta = [angle_gap * i for i in range(0, len(tRNAs))]
        radii = [sample['fractions'][tRNA] for tRNA in tRNAs]
        # for i in range(len(tRNAs)):
        #     print(i, tRNAs[i], theta[i], radii[i]*100)

        ax.scatter(theta, radii, alpha=0.75, label=sample['sample'])
    ax.legend()
    plt.show()


#tRNA_features = save_tRNA_features()
tRNA_features = read_tRNA_features()
samples_with_counts = count_tRNAs(samples, tRNA_features, output_path=output_path)
#samples_with_counts = read_samples_with_counts()
pp(samples_with_counts)
#plot_radial(samples_with_counts)
