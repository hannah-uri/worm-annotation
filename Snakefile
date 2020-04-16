STRAINS_WITH_5_INDEP_R1 = ['altadena', 'amares', 'auckland', 'bristol', 'hermanville', 'lake-forest-park', 'lisbon', 'madagascar', 'paolo-alto', 'roxel', 'salt-lake-city', 'san-fran', 'wuhan'] 
STRAINS_WITH_RNA_R1 = ['altadena', 'amares', 'auckland', 'bristol', 'hermanville', 'lake-forest-park', 'lisbon', 'madagascar', 'paolo-alto', 'roxel', 'salt-lake-city']

rule all:
	input:
		#'output/plots/trna/tRNA-counts_small-rna_5-indep_r1.png',
		#'output/plots/trna/tRNA-counts_rna_r1.png',
		'output/summary/auckland.png',
		'output/summary/altadena.png',
		'output/summary/amares.png',
		'output/summary/bristol.png',
		#'output/maker/hermanville',
		#'output/maker/lake-forest-park',
		#'output/maker/lisbon',
		#'output/maker/madagascar',
		#'output/maker/roxel',
		#'output/maker/salt-lake-city'

rule prep_ref_genome:
	output:
		'output/ref/genome/c_elegans.PRJNA13758.WS270.genomic.fa'
	conda:
		'envs/bwa.yml'
	shell:
		'wget -O {output}.gz ftp://ftp.wormbase.org/pub/wormbase/releases/WS270/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS270.genomic.fa.gz;'
		'gunzip {output}.gz;'
		'bwa index {output}'

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
        input:  type_path
	params: depth=bonus_depth
        output:
                temp('output/cat/{strain}/{strain}_r{replicate, [123]}_{type}.fastq')
        shell:
                'cat {input}{params.depth} | gunzip -c > {output}'

rule trim_small_rna:
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
	input:
		'output/cat/{strain}/{strain}_r{replicate}_rna.fastq'
	output:
		file='output/cutadapt/{strain}/{strain}_r{replicate}_rna.fastq',
		log='output/cutadapt/{strain}/{strain}_r{replicate}_rna.log'
	conda:
		'envs/cutadapt.yml'
	threads:
		4
	shell:
		'mkdir -p output/cutadapt/{wildcards.strain};'
		'cutadapt -b file:input/adaptors/truseq.fa {input} -u 10 -u -2 -j {threads} --max-n 0.5 --minimum-length 25 > {output.file} 2> {output.log}'

rule allign:
	input:
		'output/ref/genome/c_elegans.PRJNA13758.WS270.genomic.fa',
		'output/cutadapt/{strain}/{strain}_r{replicate}_{type}.fastq'
	output:
		'output/bwa/{strain}/{strain}_r{replicate, [123]}_{type}.sam'
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
	input:  list_sams
	params: return_strain_list
	conda:
		'envs/worms.yml'
	output:
		'output/plots/trna/tRNA-counts_{type}_r{replicate}.tsv'
	script:
		'scripts/radial-tRNA-counter.py'

rule plot_trnas:
	input:
		'output/plots/trna/tRNA-counts_{type}_r{replicate}.tsv'
	output:
		'output/plots/trna/tRNA-counts_{type}_r{replicate}.png'
	params:
		transformation='log2'
	script:
		'scripts/radial_plot.R'

rule trinity:
	input:
		'output/cutadapt/{strain}/{strain}_r1_rna.fastq'
	output:
		'output/trinity/{strain}/trinity_out_dir/Trinity.fasta'
	conda:
		'envs/trinity.yml'
	threads:
		60
	shell:
		'Trinity --seqType fq --max_memory 60G --single {input} --output output/trinity/{wildcards.strain}/trinity_out_dir --CPU {threads}'

rule rep_count:
	input:
		'input/strains/{strain}/pilonN2n.fa'
	output:
		'output/windowmasker/{strain}/{strain}_mask.count'
	conda:
		'envs/windowmasker.yml'
	shell:
		'windowmasker -in {input} -out {output} -mk_counts'

rule rep_mask:
	input:
		genome='input/strains/{strain}/pilonN2n.fa',
		counts='output/windowmasker/{strain}/{strain}_mask.count'
	output:
		'output/windowmasker/{strain}/{strain}_masked.fa'
	conda:
		'envs/windowmasker.yml'
	shell:
		'windowmasker -in {input.genome} -ustat {input.counts} -out {output}.tmp -outfmt fasta -parse_seqids;'
		"cat {output}.tmp | sed 's/>lcl\|chrI_pilon/I' | sed 's/>lcl\|chrII_pilon/II' | sed 's/>lcl\|chrIII_pilon/III' | sed 's/>lcl\|chrIV_pilon/IV' | sed 's/>lcl\|chrV_pilon/V' | sed 's/>lcl\|chrX_pilon/X' > {output};"
		'rm {output}.tmp'

rule maker:
	input:
		transcriptome='output/trinity/{strain}/trinity_out_dir/Trinity.fasta',
		genome='output/windowmasker/{strain}/{strain}_masked.fa'
	output:
		directory('output/maker/{strain}')
	conda:
		'envs/maker.yml'
	shell:
		'mkdir output/maker/{wildcards.strain};'
		'cd output/maker/{wildcards.strain};'
		'maker -CTL;'
		"cat maker_opts.ctl | sed 's/^genome=.*$/\genome=\.\.\/\.\.\/\.\.\/output\/windowmasker\/{wildcards.strain}\/{wildcards.strain}_masked.fa /' | sed 's/^est=.*$/est=\.\.\/\.\.\/trinity\/{wildcards.strain}\/trinity_out_dir\/Trinity.fasta /' | sed 's/^model_org=.*$/model_org= /' | sed 's/^est2genome=.*$/est2genome=1 /' > maker_opts.ctl.2;"
		'mv maker_opts.ctl.2 maker_opts.ctl;'
		'maker -qq'

rule collect_gffs:
	# collects the scattered outputs of maker into a single file and trims the fasta sequence from the end 
	input:
		'output/maker/{strain}'
	output:
		'output/annotations/{strain}.gff'
	shell:
		"echo '##gff-version 3' > {output};"
		"cat {input}/*/*/*/*/*/*[^IXV]chrI[^IXV]*.gff | sed '1,/^##FASTA/!d' | sed '/^#/d' > {output}.I.tmp;"
		"cat {input}/*/*/*/*/*/*[^IXV]chrII[^IXV]*.gff | sed '1,/^##FASTA/!d' | sed '/^#/d' > {output}.II.tmp;"
		"cat {input}/*/*/*/*/*/*[^IXV]chrIII[^IXV]*.gff | sed '1,/^##FASTA/!d' | sed '/^#/d' > {output}.III.tmp;"
		"cat {input}/*/*/*/*/*/*[^IXV]chrIV[^IXV]*.gff | sed '1,/^##FASTA/!d' | sed '/^#/d' > {output}.IV.tmp;"
		"cat {input}/*/*/*/*/*/*[^IXV]chrV[^IXV]*.gff | sed '1,/^##FASTA/!d' | sed '/^#/d' > {output}.V.tmp;"
		"cat {input}/*/*/*/*/*/*[^IXV]chrX[^IXV]*.gff | sed '1,/^##FASTA/!d' | sed '/^#/d' > {output}.X.tmp;"
		"cat {output}.I.tmp {output}.II.tmp {output}.III.tmp {output}.IV.tmp {output}.V.tmp {output}.X.tmp >> {output};"
		"rm {output}*.tmp"

rule summarise:
	input:
		'output/annotations/{strain}.gff'
	output:
		'output/summary/{strain}.png'
	script:
		'scripts/maker-analysis.R'
