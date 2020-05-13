STRAINS_WITH_5_INDEP_R1 = ['altadena', ':amares', 'auckland', 'bristol', 'hermanville', 'lake-forest-park', 'lisbon', 'madagascar', 'paolo-alto', 'roxel', 'salt-lake-city', 'san-fran', 'wuhan'] 
STRAINS_WITH_RNA_R1 = ['altadena', 'amares', 'auckland', 'bristol', 'hermanville', 'lake-forest-park', 'lisbon', 'madagascar', 'paolo-alto', 'roxel', 'salt-lake-city']
REPLICATE_DICT = {'altadena': 3, 'amares': 1, 'auckland': 1, 'bristol': 3, 'hermanville': 1, 'lake-forest-park': 3, 'lisbon': 3, 'madagascar': 3, 'roxel': 1, 'salt-lake-city': 3}
#REPLICATE_DICT = {'altadena': 3, 'amares': 1, 'auckland': 1, 'bristol': 1, 'hermanville': 1, 'lisbon': 3, 'madagascar': 1, 'roxel': 1, 'salt-lake-city': 1}
#REPLICATE_DICT = {'altadena': 3, 'amares': 1, 'auckland': 1, 'bristol': 3, 'hermanville': 1, 'lake-forest-park': 2, 'lisbon': 3, 'madagascar': 2, 'roxel': 1, 'salt-lake-city': 3}

rule all:
	input:
		'output/plots/trna/tRNA-counts_small-rna_5-indep_r1.png',
		'output/plots/trna/tRNA-counts_rna_r1.png',
		'output/summary/genes-reads.png',
		'output/summary/aeds.png'

def type_path(wildcards):
	if wildcards.type == 'rna':
		return 'input/strains/{strain}/r{replicate}/rna'
	elif wildcards.type == 'small-rna_5-dep':
		return 'input/strains/{strain}/r{replicate}/small-rna/5-dep'
	elif wildcards.type  == 'small-rna_5-indep':
		return 'input/strains/{strain}/r{replicate}/small-rna/5-indep'
	else:
		print('type: ' + wildcards.type)

def bonus_depth(wildcards):
	if wildcards.type =='rna':
		return '/*'
	else:
		return '/*/*'

rule merge:
	""" Merges the lanes from raw read files """
        input:
		type_path
	params:
		depth=bonus_depth
	output:
		temp('output/cat/{strain}/{strain}_r{replicate, [123]}_{type}.fastq')
	shell:
		'cat {input}{params.depth} | gunzip -c > {output}'

rule trim_small_rna:
	""" Trims TruSeq small RNA adaptors from 3' end, trims 1 from the start of each read, and sets the minimum read length to 10bp """
	input:
                'output/cat/{strain}/{strain}_r{replicate}_small-rna_5-{dep}.fastq'
	output:
		file='output/cutadapt/{strain}/{strain}_r{replicate}_small-rna_5-{dep}.fastq',
		log='output/cutadapt/{strain}/{strain}_r{replicate}_small-rna_5-{dep}.log'
	conda:
		'envs/cutadapt.yml'
	threads:
		4
	shell:
		'cutadapt -a TGGAATTCTCGGGTGCCAAGG {input} -u 1 -m 10 -j {threads} > {output.file} 2> {output.log}'

rule trim_rna:
	""" Trims 10 from the start & 2 from the end of each read, then trims TruSeq adaptors fom both ends, then discards read with over 50% N or under 25bp """
	input:
		'output/cat/{strain}/{strain}_r{replicate}_rna.fastq'
	output:
		file=temp('output/cutadapt/{strain}/{strain}_r{replicate}_rna.fastq'),
		read_count='output/cutadapt/{strain}/{strain}_r{replicate}.txt',
		log='output/cutadapt/{strain}/{strain}_r{replicate}_rna.log'
	conda:
		'envs/cutadapt.yml'
	threads:
		4
	shell:
		'mkdir -p output/cutadapt/{wildcards.strain}; '
		'cutadapt -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCA {input} -u 10 -u -2 -j {threads} --max-n 0.5 --minimum-length 25 > {output.file} 2> {output.log}; '
		"echo $(($(cat {output.file}| wc -l)/3)) > {output.read_count}"

rule prep_ref_genome:
	""" Downloads & indexes reference genome for alligning """
	output:
		'output/ref/genome/c_elegans.PRJNA13758.WS270.genomic.fa'
	conda:
		'envs/bwa.yml'
	shell:
		'wget -O {output}.gz ftp://ftp.wormbase.org/pub/wormbase/releases/WS270/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS270.genomic.fa.gz; '
		'gunzip {output}.gz; '
		'bwa index {output}'

rule allign:
	""" Alligns trimmed read files to reference genome """
	input:
		'output/ref/genome/c_elegans.PRJNA13758.WS270.genomic.fa',
		'output/cutadapt/{strain}/{strain}_r{replicate}_{type}.fastq'
	output:
		temp('output/bwa/{strain}/{strain}_r{replicate, [123]}_{type}.sam')
	conda:
		'envs/bwa.yml'
	threads:
		4
	shell:
		'bwa mem -t {threads} {input} > {output}'

def list_sams(wildcards):
	if wildcards.replicate == '1' and wildcards.type == 'small-rna_5-indep':
		return expand('output/bwa/{strain}/{strain}_r1_small-rna_5-indep.sam', strain=STRAINS_WITH_5_INDEP_R1)
	elif wildcards.replicate == '1' and wildcards.type =='rna':
		return expand('output/bwa/{strain}/{strain}_r1_rna.sam', strain=STRAINS_WITH_RNA_R1)
	else:
		print('replicate: ' + wildcards.replicate + ', type: ' + wildcards.type)

def return_strain_list(wildcards):
	if wildcards.replicate == '1' and wildcards.type == 'small-rna_5-indep':
		return STRAINS_WITH_5_INDEP_R1
	elif wildcards.replicate == '1' and wildcards.type =='rna':
		return STRAINS_WITH_RNA_R1


rule count_trnas:
	""" Classifies & counts the tRNAs """
	input:
		list_sams
	params:
		return_strain_list
	conda:
		'envs/worms.yml'
	output:
		'output/plots/trna/tRNA-counts_{type}_r{replicate}.tsv'
	script:
		'scripts/radial-tRNA-counter.py'

rule plot_trnas:
	""" Plots the distribution of tRNA types across the strains """ 
	input:
		'output/plots/trna/tRNA-counts_{type}_r{replicate}.tsv'
	output:
		'output/plots/trna/tRNA-counts_{type}_r{replicate}.png'
	params:
		transformation='log2'
	script:
		'scripts/radial_plot.R'

rule trinity:
	""" Assembles the transcriptome from trimmed RNA seq reads """
	input:
		'output/cutadapt/{strain}/{strain}_r{replicate}_rna.fastq'
	output:
		'output/trinity/{strain}/r{replicate}/trinity_out_dir/Trinity.fasta'
	conda:
		'envs/trinity.yml'
	threads:
		60
	log:
		'output/trinity/{strain}/r{replicate}/trinity.log'
	shell:
		'Trinity --seqType fq --max_memory 60G --single {input} --output output/trinity/{wildcards.strain}/r{wildcards.replicate}/trinity_out_dir --CPU {threads} > {log}'

rule count_repeats:
	""" Counts the repeats from the genome """
	input:
		'input/strains/{strain}/*.fa'
	output:
		temp('output/windowmasker/{strain}/{strain}_mask.count')
	conda:
		'envs/windowmasker.yml'
	shell:
		'windowmasker -in {input} -out {output} -mk_counts'

rule mask_repeats:
	""" Masks the repeats in the genome """
	input:
		genome='input/strains/{strain}/pilonN2n.fa',
		counts='output/windowmasker/{strain}/{strain}_mask.count'
	output:
		'output/windowmasker/{strain}/{strain}_masked.fa'
	conda:
		'envs/windowmasker.yml'
	shell:
		'windowmasker -in {input.genome} -ustat {input.counts} -out {output}.tmp -outfmt fasta -parse_seqids;'
		"cat {output}.tmp |  sed 's/lcl.//g' > {output};"
		'rm {output}.tmp'

rule download_ref_annotations:
	output:
		'output/ref/annotation/c_elegans.PRJNA13758.WS270.annotations.gff3'
	shell:
		'wget -O {output}.gz ftp://ftp.wormbase.org/pub/wormbase/releases/WS270/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS270.annotations.gff3.gz; '
		'gunzip {output}.gz'

rule maker_single:
	""" Creates a novel annotation, combining evidence from the transcriptome and ab initio predictors, using names from the reference annotation where possible """
	input:
		transcriptome='output/trinity/{strain}/r{replicate}/trinity_out_dir/Trinity.fasta',
		genome='output/windowmasker/{strain}/{strain}_masked.fa',
		reference_annotations='output/ref/annotation/c_elegans.PRJNA13758.WS270.annotations.gff3'
	output:
		temp(directory('output/maker/{strain}/r{replicate}'))
	conda:
		'envs/maker.yml'
	threads:
		60
	shell:
		'mkdir {output};'
		'cd {output};'
		'maker -CTL;'
		"cat maker_opts.ctl | "
		
		"sed 's/^model_gff=.*$/model_gff=\.\.\/\.\.\/\.\.\/ref\/annotation\/c_elegans.PRJNA13758.WS270.annotations.gff3 /' | " # uses gene models (names) from reference
		"sed 's/^trna=.*$/trna=1 /' | " # finds tRNAs using tRNAscan
		
		"sed 's/^genome=.*$/genome=\.\.\/\.\.\/\.\.\/windowmasker\/{wildcards.strain}\/{wildcards.strain}_masked.fa /' | " # genome to be annotated
		"sed 's/^est=.*$/est=\.\.\/\.\.\/\.\.\/trinity\/{wildcards.strain}\/r{wildcards.replicate}\/trinity_out_dir\/Trinity.fasta /' | " # transcriptome evidence
		"sed 's/^model_org=.*$/model_org= /' | " # disables RepeatMakser (was broken by licensing in 2019)
		"sed 's/^est2genome=.*$/est2genome=1 /' > maker_opts.ctl.2;"
		'mv maker_opts.ctl.2 maker_opts.ctl;'
		'maker -qq -cpus {threads}'

rule maker_double:
	""" Improves a previous annotation with data from replicates """
	input:
                premade='output/annotations/{strain}/r{replicate1}/{strain}.gff',
                transcriptome='output/trinity/{strain}/r{replicate2}/trinity_out_dir/Trinity.fasta',
                genome='output/windowmasker/{strain}/{strain}_masked.fa',
		reference_annotations='output/ref/annotation/c_elegans.PRJNA13758.WS270.annotations.gff3'
	output:
		temp(directory('output/maker/{strain}/r{replicate1}+{replicate2}'))
	conda:
                'envs/maker.yml'
	threads:
		60
	shell:
                'mkdir {output};'
                'cd {output};'
                'maker -CTL;'
                "cat maker_opts.ctl | "
		
		"sed 's/^maker_gff=.*$/maker_gff=\.\.\/\.\.\/\.\.\/annotations\/{wildcards.strain}\/r{wildcards.replicate1}\/{wildcards.strain}.gff /' | " # maker file to be improved with new data
                "sed 's/^est_pass=.*$/est_pass=1 /' | " # uses previous transcriptome evidence
                "sed 's/^pred_pass=.*/pred_pass=0 /' | " # re-generates ab initio predictions 
                "sed 's/^map_forward=.*/map_forward=1 /' | " # conserves names from previous run
                "sed 's/^model_pass=.*/model_pass=1 /' | "   # conserves names from previous run

		"sed 's/^model_gff=.*$/model_gff=\.\.\/\.\.\/\.\.\/ref\/annotation\/c_elegans.PRJNA13758.WS270.annotations.gff3 /' | " # uses gene models (names) from reference
                "sed 's/^trna=.*$/trna=1 /' | " # finds tRNAs using tRNAscan
                
		"sed 's/^genome=.*$/genome=\.\.\/\.\.\/\.\.\/windowmasker\/{wildcards.strain}\/{wildcards.strain}_masked.fa /' | " # genome to be annotated
		"sed 's/^est=.*$/est=\.\.\/\.\.\/\.\.\/trinity\/{wildcards.strain}\/r{wildcards.replicate2}\/trinity_out_dir\/Trinity.fasta /' | " # new transcritome evidence
                "sed 's/^model_org=.*$/model_org= /' | " # disables RepeatMakser (was broken by licensing in 2019)
                "sed 's/^est2genome=.*$/est2genome=1 /' > maker_opts.ctl.2;"
                'mv maker_opts.ctl.2 maker_opts.ctl;'
                'maker -qq -c {threads}'


rule collect_gffs:
	""" Collects the scattered outputs of MAKER into a single file and trims the fasta sequence from the ends """
	input:
		'output/maker/{strain}/r{replicate}'
	output:
		'output/annotations/{strain}/r{replicate}/{strain}.gff'
	shell:
		"echo '##gff-version 3' > {output};"
		"cat {input}/*/*/*/*/*/*chrI[^IXV]*.gff | sed '1,/^##FASTA/!d' | sed '/^#/d' > {output}.I.tmp;"
		"cat {input}/*/*/*/*/*/*chrII[^IXV]*.gff | sed '1,/^##FASTA/!d' | sed '/^#/d' > {output}.II.tmp;"
		"cat {input}/*/*/*/*/*/*chrIII[^IXV]*.gff | sed '1,/^##FASTA/!d' | sed '/^#/d' > {output}.III.tmp;"
		"cat {input}/*/*/*/*/*/*chrIV[^IXV]*.gff | sed '1,/^##FASTA/!d' | sed '/^#/d' > {output}.IV.tmp;"
		"cat {input}/*/*/*/*/*/*chrV[^IXV]*.gff | sed '1,/^##FASTA/!d' | sed '/^#/d' > {output}.V.tmp;"
		"cat {input}/*/*/*/*/*/*chrX[^IXV]*.gff | sed '1,/^##FASTA/!d' | sed '/^#/d' > {output}.X.tmp;"
		"cat {output}.I.tmp {output}.II.tmp {output}.III.tmp {output}.IV.tmp {output}.V.tmp {output}.X.tmp >> {output};"
		"rm {output}*.tmp"

rule summarise_individual:
	""" Counts the number of each feature and the reported evidence """
	input:
		'output/annotations/{strain}/r{replicate}/{strain}.gff'
	output:
		graph='output/summary/{strain}_r{replicate}.png',
		table='output/summary/{strain}_r{replicate}.tsv'
	script:
		'scripts/maker-analysis.R'


def return_reads_count_files(wildcards):
	read_count_files = []
	for strain, replicate_count in REPLICATE_DICT.items():
		for replicate in range(1, replicate_count + 1):
			read_count_files.append('output/cutadapt/' + strain + '/' + strain + '_r' + str(replicate) + '.txt')
	return read_count_files

def return_individual_summary_files(wildcards):
	individual_summary_files = []
	for strain, replicate_count in REPLICATE_DICT.items():
		for replicate in range(1, replicate_count + 1):
			individual_summary_files.append('output/summary/' + strain + '_r' + '+'.join([str(r) for r in range(1, replicate + 1)]) + '.tsv')
	return individual_summary_files

rule summarise_all:
	""" Collects all the summary files into a single table, along with a count of the number of trimmed reads used """
	input:
		return_reads_count_files,
		return_individual_summary_files
	output:
		'output/summary/total.csv'
	script:
		'scripts/group_summarise.py'

rule graph_genes_reads:
	""" Graphs the genes annotated vs the number of trimmed reads used """
	input:
		'output/summary/total.csv'
	output:
		'output/summary/genes-reads.png'
	script:
		'scripts/genes-reads.R'

def return_best_annotation_files(wildcards):
	annotation_files = []
	for strain, replicate_count in REPLICATE_DICT.items():
		annotation_files.append('output/annotations/' + strain + '/r' + '+'.join([str(r) for r in range(1, replicate_count + 1)]) + '/' + strain + '.gff')
	return annotation_files

rule count_aed:
	""" Counts the MAKER -reported AED scores from each transcript """
	input:
		return_best_annotation_files
	output:
		'output/summary/aeds.tsv'
	script:
		'scripts/aed_score.py'


rule graph_aed:
	""" Graphs the AED curves """
	input:
		'output/summary/aeds.tsv'
	output:
		'output/summary/aeds.png'
	script:
		'scripts/aed_graph.R'
