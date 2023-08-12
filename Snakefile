import re
import sys

run_files = []
sample_files = []
sample_setup_dict = {}
sample_run_dict = {}
sample_run_setup = config['sample_fastq_table']
# format of study_fq_tsv
# tab separated
# setup is either "paired" or "single"
# header of: study_id	sample_id	setup
# example line:
# zander	D3C_D0_1__HJLT7DSX3_19270591_S73_L002	paired
for line in open(sample_run_setup):
	if len(line.split('\t')) != 3:
		sys.exit(print('\n****************\nONLY ' + str(len(line.split())) + \
					' tab separated values in config[\'sample_fastq_table\']!!!!\n***************\n'))
	if 'sample_accession' in line:
		continue
	line = line.split('\t')
	if line[1] in run_files:
		sys.exit(print("Duplicated run: " + line[1] + "\nCheck sample input file for duplicates"))
	run_files.append(line[1])
	if line[0] not in sample_files:
		sample_files.append(line[0])
	sample_setup_dict[line[0]] = line[2].strip()
	if line[0] not in sample_run_dict.keys():
		sample_run_dict[line[0]] = [line[1]]
	else:
		oval = sample_run_dict[line[0]]
		oval.append(line[1])
		sample_run_dict[line[0]] = oval


fq_dir = config['fastq_dir']

def return_fq(prefix):
	setup = sample_setup_dict[prefix]
	out = []
	if setup == 'single':
		out.append(fq_dir + '/' + prefix + config['fqS_suffix'])
	else:
		out.append(fq_dir + '/' + prefix + config['fq1_suffix'])
		out.append(fq_dir + '/' + prefix + config['fq2_suffix'])
	return(out)

def salmon_input_maker(prefix):
	setup = sample_setup_dict[prefix]
	out = []
	if setup == 'single':
		out.append('-r ' + fq_dir + '/' + prefix + config['fqS_suffix'])
	else:
		out.append('-1 ' + fq_dir + '/' + prefix + config['fq1_suffix'])
		out.append('-2 ' + fq_dir + '/' + prefix + config['fq2_suffix'])
	return(out)

	

wildcard_constraints:
	run = '|'.join(run_files),
	sample = '|'.join(sample_files)

SALMON_QUANT_OUTPUT = ['salmon_quant/' + sample + '/quant.sf'  \
	for sample in sample_files] 
	

localrules: all, print_quants

rule all:
	input:
		'counts/gene_counts.csv.gz',
		'salmon_quant_output.tsv'

rule aggregate_fq:
	# merge run level files to the sample level
	input:
		lambda wildcards: expand(fq_dir + '/{run}{{suffix}}', run = sample_run_dict[wildcards.sample])
	output:
		fq_dir + '/{sample}{suffix}'
	shell:
		"""
		cat {input} > {output}
		"""

rule trim_fastq_se:
	input:
		fq_dir + '/{sample}' + config['fqS_suffix']
	output:
		temp('fastq_trimmed/{sample}' + config['fqS_suffix']),
	threads: 8
	conda: 'eiad_rna_quant.yaml'
	shell:
		"""
		fastp --thread {threads} -i {input} -o {output} 
		"""

rule trim_fastq_pe:
	input:
		r1 = ancient(fq_dir + '/{sample}' + config['fq1_suffix']),
		r2 = ancient(fq_dir + '/{sample}' + config['fq2_suffix'])
	output:
		r1 = temp(fq_dir + '/fastq_trimmed/{sample}' + config['fq1_suffix']),
		r2 = temp(fq_dir + '/fastq_trimmed/{sample}' + config['fq2_suffix']),
	threads: 8
	conda: 'eiad_rna_quant.yaml'
	shell:
		"""
		fastp --thread {threads} -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2}
		"""

rule download_salmon_index:
	input:
		url = config['refgenie_url']	
	output:
		directory('default')
	conda: 'eiad_rna_quant.yaml'
	shell:
		"""
		wget -O index.tgz {input}
		tar -xzvf index.tgz
		rm index.tgz
		"""

rule salmon_quant:
	input:
		fq_files = lambda wildcards: return_fq(wildcards.sample),
		index = 'default'
	output:
		'salmon_quant/{sample}/quant.sf'
	params:
		input_command = lambda wildcards: salmon_input_maker(wildcards.sample),
		dir = 'salmon_quant/{sample}'
	threads: 16
	conda: 'eiad_rna_quant.yaml'
	shell:
		"""
		salmon quant --index {input.index} --libType A --seqBias --gcBias \
			--validateMappings \
			-p {threads} \
			{params.input_command} \
			-o {params.dir}
		""" 

rule print_quants:
	input:
		SALMON_QUANT_OUTPUT
	output:
		'salmon_quant_output.tsv'
	run:
		with open(output[0], 'w') as f:
			f.write('\n'.join(str(quant) for quant in input))

rule make_counts:
	input:
		SALMON_QUANT_OUTPUT,
	output:
		'counts/gene_counts.csv.gz'
	params: config['R_make_counts_path']
	conda: 'eiad_rna_quant_r.yaml'
	shell:
		"""
		Rscript {params}
		"""
