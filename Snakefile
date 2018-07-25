import os
import subprocess as sp
import json
import re

'''
required modules: bedtools

config file:
    ids=text file with sraids
    bam_files=file with bam_locs
    bam_path= path/to/bams/
    fastq-location to fodler with other fastq's

'''
#adapted from https://www.biostars.org/p/139422/
#for paired end data, it will have 8 lines for 1 spot entry as per SRA format
def isPaired(filename):
    # this makes me sad but it works
    tmp=open('tmp.txt','w+')
    k= sp.run(["fastq-dump","-X","1","-Z","--split-spot", str(filename)],stdout=tmp)
    #        contents =
    tmp.close()
    tmp=open('tmp.txt','r')
    count=0
    for line in tmp:
        count+=1
    tmp.close()
    print(count)
    if count==4:

        return("-r")
    else:
        return("-2")
#new
def isFastqPaired(filename):
    with open(filename) as fastq:
        head=[next(fastq).split(' ')[0] for x in range(8)]
    print(head)
    first=head[0]
    head.pop(0)
    if first in head:
        # if the same spot is present twice must be paired end
        return('-2')
    else:
        return('-r')
#open file with all ids
SRA_IDS=[]
bam_files=[]# contains a list of bam names
salmonindex='ref/salmonindex'
with open(config['ids']) as file:
    for line in file:
        SRA_IDS.append(str(line).strip('"|\n'))
# witih open(config['bam_files']) as bams:
#     for line in bams:
#         tmp=str(line).strip('"|\n')
#         bam_files.append(tmp)
salmonindex='ref/salmonindex'
#define target files
rule all:
    input:
        expand("quant_files/{sraids}/quant.sf",sraids=SRA_IDS),
        'ref/gencodeRef.fa.gz',
        directory(salmonindex)
#download files using 
rule download_fastqs:
    output:'fastq_files/{srs}.fastq'
    run:
        shell("module load {config[sqlite_version]}")
        shell("module load {config[sratoolkit_version]}")
        with open(config['ids'],'r') as ids:
            for srs_id in ids:
                srs_id=srs_id.strip('\n')
                command='sqlite3 ' + config['sqlite_file'] + ' "SELECT run_accession FROM sra WHERE sample_accession=\'{}\'"'.format(srs_id)
                run_ids=sp.check_output(command, shell=True).decode('utf-8').strip("'|\n").split('\n')# outputs a list of run_ids
                for run in run_ids:# download all runs
                    sp.run("fastq-dump -X 5 -Z {}  > {}.{}.fqp ".format(run,srs_id,run),shell=True)#***REMOVE -X 5***
                sp.run('cat {}* > fastq_files/{}.fastq'.format(srs_id,srs_id),shell=True)
                sp.run('rm *.fqp',shell=True)
# rule bam_to_fastq:
#     input:[config['bam_path']+bamname for bamname in bam_files ]
#     output:['fastq_files/'+bamname+'fastq' for bamname in bam_files]
#     run:
#         for i in range(len(bam_files)):
#             command='bedtools bamtofastq -i {} -fq {}'.format(input[i],output[i])
#             sp.run(command,shell=True)
#download fasta reference and build salmon index_seq_hash
#might add webaddress and salmonindex same as config arguments
rule build_salmon_index:
    output:['ref/gencodeRef.fa.gz',directory(salmonindex)]
    run:
        print('running index download')
        sp.run('wget -O ref/gencodeRef.fa.gz -nc ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/gencode.v26.transcripts.fa.gz',shell=True)
        salmonindexcommand='salmon index -t {} -i {} --type quasi -k 31'.format(output[0],output[1])
        sp.run(salmonindexcommand, shell=True)

###############
#ALIGNMENT
#############


# snakemake requires multiple log files if you use wildcards
if not os.path.exists('logs'):
    log="logs/failed_mapping.log"
    sp.run('mkdir logs; touch logs/failed_mapping.log', shell=True)

rule run_salmon_SRA:
    input:"fastq_files/{sraid}.fastq",directory(salmonindex)
    output:"quant_files/{sraid}/quant.sf"
    run:
        pairedFlag=isFastqPaired(wildcards.sraid+'.fastq')
        salmon_command='salmon quant -i {} -l A  {} {} -o quant_files/{}'.format(input[1], pairedFlag,input[0],wildcards.sraid)
        sp.run(salmon_command,shell=True)
        with open('quant_files/{}/aux_info/meta_info.json'.format(wildcards.sraid)) as file:
            salmonLog=json.load(file)
            mappingscore=salmonLog["percent_mapped"]
        if mappingscore <= 50:
            with open(log,'a') as logFile:
                logFile.write('Sample {} failed QC mapping Percentage: {}'.format(wildcards.sraid,mappingscore))
#need data to test this on

# rule run_salmon_BAMS:
#     input:'{bam_fastqs}.fastq',directory(salmonindex)
#     output:'{bam_fastqs}/quant.sf'
#     run:
#         pairedFlag='-r'# need to change this
#         salmon_command='salmon quant -i {} -l A  {} {} -o {}'.format(input[1], pairedFlag,input[0],wildcards.bam_fastqs)
#         sp.run(salmon_command,shell=True)
#         with open('{}/aux_info/meta_info.json'.format(wildcards.bam_fastqs)) as file:
#             salmonLog=json.load(file)
#             mappingscore=salmonLog["percent_mapped"]
#         if mappingscore <= 50:
#             with open(log,'a') as logFile:
#                 logFile.write('Sample {} failed QC mapping Percentage: {}'.format(wildcards.bam_fastqs,mappingscore))

