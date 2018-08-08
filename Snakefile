'''
config file:
    sampleFile: tab seperated file  containing sample info
    refFasta_url: link to refernece fastq_path
    salmon_version:
    sratoolkit_version:
notes:
-the 5.2 version requires specifying directorys in output section of rule iwth directory(). Biowulf currently using 5.1
-need to make a rule to download all Gencode refs
'''
import subprocess as sp
def readSampleFile(samplefile):
    res={}
    with open(samplefile) as file:
        for line in file:
            info=line.strip('\n').split('\t')
            res[info[0]]={'files':info[1].split(','),'paired':True if info[2]=='y' else False}
    return(res)
def lookupRunfromID(card,sample_dict):
    id=card
    if '_' in id:
        i= '1' if id[-1]=='1' else '2'# check L/R file
        id=card[:-2]
    fqpfiles=sample_dict[id]['files']
    res=[]
    for file in fqpfiles:
        if sample_dict[id]['paired']:
            #PE
            res.append('fastqParts/{}_{}.fastq.gz'.format(file,i))
        else:
            #SE
            res.append('fastqParts/{}.fastq.gz'.format(file))
    return(res)
configfile:'config.yaml'
sample_dict=readSampleFile(config['sampleFile'])# sampleID:dict{path,paired,metadata}
sample_names=sample_dict.keys()
loadSRAtk="module load {} && ".format(config['sratoolkit_version'])
loadSalmon= "module load {} && ".format(config['salmon_version'])
salmonindex='ref/salmonindex'
STARindex='ref/STARindex'
ref_fasta='ref/gencodeRef.fa'
ref_GTF='ref/gencodeAno.gtf'
ref_PA='ref/gencodePA.fa'
badruns='badruns'
rule all:
    input: expand('quant_files/{sampleID}/quant.sf', sampleID=sample_names),
            ref_fasta,ref_GTF,ref_PA,
            salmonindex
rule downloadGencode:
    output:ref_fasta,ref_GTF,ref_PA
    shell:
        '''
        wget -O ref/gencodeRef.fa.gz {config[refFasta_url]}
        wget -O ref/gencodeAno.gtf.gz {config[refGTF_url]}
        wget -O ref/gencodePA.fa.gz {config[refPA_url]}
        gunzip ref/gencodeRef.fa.gz
        gunzip ref/gencodeAno.gtf.gz
        gunzip ref/gencodePA.fa.gz
        '''


rule getFQP:
    output: temp('fastqParts/{id}.fastq.gz')
    log: 'logs/{id}.dbg'
    run:
        id=wildcards.id
        id=id[:-2] if '_'in id else id #idididid
        try:
            sp.check_output(loadSRAtk + 'fastq-dump --gzip --split-3 -O fastqParts {}'.format(id),shell=True)
        except sp.CalledProcessError:
            with open('logs/{}.dbg'.format(wildcards.id),'a') as br:
                br.write(id)
rule aggFastqsPE:
    input:lambda wildcards:lookupRunfromID(wildcards.sampleID,sample_dict)
    output:temp('fastq_files/{sampleID}.fastq.gz')
    run:
        #this can use some cleaning up
        id=wildcards.sampleID
        fileParts=lookupRunfromID(id,sample_dict)
        i='1' if '_' in id and id[-1]=='1' else '2'# which strand
        id=id[:-2] if '_' in id else id
        for fqp in fileParts:
            if sample_dict[id]['paired']:
                sp.run('cat {fqp} >> fastq_files/{id}_{i}.fastq.gz '.format(fqp=fqp,i=i,id=id),shell=True)
            else:
                sp.run('cat {fqp} >> fastq_files/{id}.fastq.gz'.format(fqp=fqp,id=id),shell=True)
rule build_salmon_index:
    input: ref_fasta
    output:salmonindex
    run:
        salmonindexcommand=loadSalmon + 'salmon index -t {} --gencode -i {} --type quasi -k 31'.format(ref_fasta,salmonindex)
        sp.run(salmonindexcommand, shell=True)
# rule build_STARindex:
#     input: ref_PA, ref_GTF
#     output:STARindex,'ref/gencodeAno.gtf', 'ref/gencodePA.fa'
#     shell:
#         '''
#         module load STAR
#         wget -O ref/gencodeAno.gtf.gz -nc {config[refGTF_url]}
#         gunzip ref/gencodeAno.gtf.gz
#         wget -O ref/gencodePA.fa.gz -nc {config[refPA_url]}
#         gunzip ref/gencodePA.fa.gz
#         STAR --runThreadN 4 --runMode genomeGenerate --genomeDir {output[0]} --genomeFastaFiles {output[2]} --sjdbGTFfile {output[1]} --sjdbOverhang 100
#
#         '''



rule run_salmon:
    input: lambda wildcards: ['fastq_files/{}_1.fastq.gz'.format(wildcards.sampleID),'fastq_files/{}_2.fastq.gz'.format(wildcards.sampleID)] if sample_dict[wildcards.sampleID]['paired'] else 'fastq_files/{}.fastq.gz'.format(wildcards.sampleID),
            salmonindex
    output:'quant_files/{sampleID}/quant.sf'
    log: 'logs/{sampleID}.log'
    run:
        id=wildcards.sampleID
        paired=sample_dict[id]['paired']
        if paired:
            salmon_command=loadSalmon + 'salmon quant -i {} -l A --gcBias -1 {} -2 {} -o quant_files/{}'.format(input[2],input[0],input[1],id)
        else:
            salmon_command=loadSalmon + 'salmon quant -i {} -l A --gcBias -r {} -o quant_files/{}'.format(input[1],input[0],id)
        sp.run(salmon_command,shell=True)
        log1='logs/{}.log'.format(id)
        with open('quant_files/{}/aux_info/meta_info.json'.format(id)) as file:
            salmonLog=json.load(file)
            mappingscore=salmonLog["percent_mapped"]
        if mappingscore <= 50:
            with open(log1,'w+') as logFile:
                logFile.write('Sample {} failed QC mapping Percentage: {}'.format(id,mappingscore))
