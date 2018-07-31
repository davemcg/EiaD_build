'''
config file:
    sampleFile: tab seperated file  containing sample info
    refFasta_url: link to refernece fastq_path

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
agg_dict={}
for sample in sample_dict.keys():
    for file in sample_dict[sample]:
        agg_dict[file]=sample
fastq_parts=agg_dict.keys()
toDownload=[x for x in fastq_parts if '.' not in x ]# get SRRids
sample_names=sample_dict.keys()
samples_PE=[x for x in sample_dict.keys() if sample_dict[x]['paired'] ]
samples_SE=[x for x in sample_dict.keys() if not sample_dict[x]['paired'] ]
parts_PE=[x for x in fastq_parts if sample_dict[agg_dict[x]]['paired']]
parts_SE=[x for x in fastq_parts if sample_dict[agg_dict[x]]['paired']]
loadSRAtk="module load {} && ".format(config['sratoolkit_version'])
salmonindex='ref/salmonindex'
loadSalmon= "module load {} && ".format(config['salmon_version'])

rule all:
    input: expand('quant_files/{sampleID}/quant.sf', sampleID=sample_names),
            'ref/gencodeRef.fa.gz',
            salmonindex

rule getFQP:
    output: 'fastqParts/{id}.fastq.gz'
    run:
        id=wildcards.id
        id=id[:-2] if '_'in id else id 
        sp.run(loadSRAtk + 'fastq-dump -X 5 --gzip --split-3 -O fastqParts {}'.format(id),shell=True)


rule aggFastqsPE:
    input:lambda wildcards:lookupRunfromID(wildcards.sampleID,sample_dict)
    output:'fastq_files/{sampleID}.fastq.gz'
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
    output:['ref/gencodeRef.fa.gz',salmonindex]
    run:
        sp.run('wget -O ref/gencodeRef.fa.gz -nc {}'.format(config['refFasta_url']),shell=True)
        salmonindexcommand=loadSalmon + 'salmon index -t {} -i {} --type quasi -k 31'.format(output[0],output[1])
        sp.run(salmonindexcommand, shell=True)

if not os.path.exists('logs'):
    log="logs/failed_mapping.log"
    sp.run('mkdir -p logs; touch logs/failed_mapping.log', shell=True)

rule run_salmon:
    input: lambda wildcards: ['fastq_files/{}_1.fastq.gz'.format(wildcards.sampleID),'fastq_files/{}_2.fastq.gz'.format(wildcards.sampleID)] if sample_dict[wildcards.sampleID]['paired'] else 'fastq_files/{}.fastq.gz'.format(wildcards.sampleID),
            salmonindex
    output:'quant_files/{sampleID}/quant.sf'
    run:
        id=wildcards.sampleID
        paired=sample_dict[id]['paired']
        if paired:
            salmon_command=loadSalmon + 'salmon quant -i {} -l A -1 {} -2 {} -o quant_files/{}'.format(input[2],input[0],input[1],id)
        else:
            salmon_command=loadSalmon + 'salmon quant -i {} -l A -r {} -o quant_files/{}'.format(input[1],input[0],id)
        sp.run(salmon_command,shell=True)
        log='logs/failed_mapping.log'
        with open('quant_files/{}/aux_info/meta_info.json'.format(id)) as file:
            salmonLog=json.load(file)
            mappingscore=salmonLog["percent_mapped"]
        if mappingscore <= 50:
            with open(log,'a') as logFile:
                logFile.write('Sample {} failed QC mapping Percentage: {}'.format(id,mappingscore))
