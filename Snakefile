'''
config file:
    sampleFile: tab seperated file  containing sample info
    refFasta_url: link to refernece fastq_path
    salmon_version:
    sratoolkit_version:
notes:
-the 5.2 version requires specifying directorys in output section of rule iwth directory(). Biowulf currently using 5.1
-need to make a rule to download all Gencode refs

***IF YOU CHANGE A RULE NAME MAKE SURE TO CHECK cluster.json ****
'''
import subprocess as sp
import itertools as it

def readSampleFile(samplefile):
    # returns a dictionary of dictionaries where first dict key is sample id and second dict key are sample  properties
    res={}
    with open(samplefile) as file:
        for line in file:
            info=line.strip('\n').split('\t')
            res[info[0]]={'files':info[1].split(','),'paired':True if info[2]=='y' else False, 'tissue':info[3],'subtissue':info[4]}
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

def genrMATsinput(subtissue):
    all_combs= list(it.combinations(subtissue,2))
    res=['rmats_out/{}_VS_{}'.format(x[0],x[1]) for x in all_combs]
    return(res)

def fastq_for_rMATS(tissue,samp_dict ,gz=True):
    # given a tissue type, return the required fastqParts
    res=[]
    if gz:
        for key in samp_dict.keys():
            if samp_dict[key]['subtissue']==tissue and samp_dict[key]['paired']:
                res.append('{}_1.fastq.gz'.format(key))
                res.append('{}_2.fastq.gz'.format(key))
    else:
        for key in samp_dict.keys():
            if samp_dict[key]['subtissue']==tissue and samp_dict[key]['paired']:
                res.append('tmp/{}/{}_1.fastq'.format(tissue,key))
                res.append('tmp/{}/{}_2.fastq'.format(tissue,key))
    return(res)

def all_fastqs(samp_dict):
    res=[]
    for sample in samp_dict.keys():
        if samp_dict[sample]['paired']:
            res.append('fastq_files/{}_1.fastq.gz'.format(sample))
            res.append('fastq_files/{}_2.fastq.gz'.format(sample))
        else:
            res.append('fastq_files/{}.fastq.gz'.format(sample))
    return(res)

configfile:'config.yaml'
sample_dict=readSampleFile(config['sampleFile'])# sampleID:dict{path,paired,metadata}
# need to add something for defined for subtissues
subtissue=["Retina_Adult.Tissue",  "RPE_Cell.Line", "ESC_Stem.Cell.Line", "RPE_Adult.Tissue",'body']
sample_names=sample_dict.keys()
loadSRAtk="module load {} && ".format(config['sratoolkit_version'])
loadSalmon= "module load {} && ".format(config['salmon_version'])
salmonindex='ref/salmonindex'
salmonindex_trimmed='ref/salmonindex_trimmed'
STARindex='ref/STARindex'
ref_fasta='ref/gencodeRef.fa'
ref_GTF='ref/gencodeAno.gtf'
ref_GTF_basic='ref/gencodeAno_bsc.gtf'
ref_GTF_PA='ref/gencodeAno_pa.gtf'
ref_PA='ref/gencodePA.fa'
badruns='badruns'
ref_trimmed='ref/gencodeRef_trimmed.fa'

rule all:
    input: expand('RE_quant_files/{sampleID}',sampleID=sample_names),expand('ref/{st}.tmp',st=subtissue),genrMATsinput(subtissue)
    #,'smoothed_filtered_tpms.csv'

rule downloadGencode:
    output:ref_fasta,ref_GTF_basic,ref_GTF_PA,ref_PA
    shell:
        '''
        wget -O ref/gencodeRef.fa.gz {config[refFasta_url]}
        wget -O ref/gencodeAno_bsc.gtf.gz {config[refGTF_basic_url]}
        wget -O ref/gencodeAno_pa.gtf.gz {config[refGTF_pa_url]}
        wget -O ref/gencodePA.fa.gz {config[refPA_url]}
        gunzip ref/gencodeRef.fa.gz
        gunzip ref/gencodeAno_bsc.gtf.gz
        gunzip ref/gencodeAno_pa.gtf.gz
        gunzip ref/gencodePA.fa.gz
        '''

rule getFQP:
    output: temp('fastqParts/{id}.fastq.gz')
    run:
        id=wildcards.id
        id=id[:-2] if '_'in id else id #idididid
        try:
            sp.check_output(loadSRAtk + 'fastq-dump --gzip --split-3 -O fastqParts {}'.format(id),shell=True)
        except sp.CalledProcessError:
            with open('logs/{}.fqp'.format(wildcards.id)) as l:
                l.write('{} did not download'.format(wildcards.id))

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
    input: ref_fasta
    output:salmonindex
    run:
        salmonindexcommand=loadSalmon + 'salmon index -t {} --gencode -i {} --type quasi -k 31'.format(ref_fasta,salmonindex)
        sp.run(salmonindexcommand, shell=True)

rule build_STARindex:
    input: ref_PA, ref_GTF_basic
    output:STARindex
    shell:
        '''
        module load STAR
        mkdir -p ref/STARindex
        STAR --runThreadN 4 --runMode genomeGenerate --genomeDir {output[0]} --genomeFastaFiles {input[0]} --sjdbGTFfile {input[1]} --sjdbOverhang 100

        '''

rule preprMats_running:
    input: all_fastqs(sample_dict)
    output:expand('ref/{st}.rmats.txt',st=subtissue),expand('ref/{st}.tmp',st=subtissue)
    shell:
        '''
        module load R
        Rscript scripts/preprMATS.R
        '''

rule run_salmon:
    input: lambda wildcards: ['fastq_files/{}_1.fastq.gz'.format(wildcards.sampleID),'fastq_files/{}_2.fastq.gz'.format(wildcards.sampleID)] if sample_dict[wildcards.sampleID]['paired'] else 'fastq_files/{}.fastq.gz'.format(wildcards.sampleID),
            salmonindex
    output: 'quant_files/{sampleID}'
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
        salmon_info='quant_files/{}/aux_info/meta_info.json'.format(id)
        if os.path.exists(salmon_info):
            with open(salmon_info) as file:
                salmonLog=json.load(file)
                mappingscore=salmonLog["percent_mapped"]
            if mappingscore <= 50:
                with open(log1,'w+') as logFile:
                    logFile.write('Sample {} failed QC mapping Percentage: {}'.format(id,mappingscore))
        else:
            with open(log1,'w+') as logFile:
                logFile.write('Sample {} failed to align'.format(id))

rule unzip_fastqs:
    input:'ref/{tissue}.tmp'
    output: temp('tmp/{tissue}/{fastqname}')
    run:
        to_unzip= fastq_for_rMATS(wildcards.tissue,sample_dict,)
        for file in to_unzip:
            sp.run('gunzip fastq_files/{} -c > tmp/{}/{}'.format(file,wildcards.tissue,file.strip('.gz')),shell=True)

rule runrMATS:
    input: lambda wildcards: fastq_for_rMATS(wildcards.tissue1,sample_dict,gz=False)+fastq_for_rMATS(wildcards.tissue2,sample_dict,gz=False),
            STARindex
    output: 'rmats_out/{tissue1}_VS_{tissue2}',
    # might have to change read length to some sort of function
    shell:
        '''
        module load rmats
        rmats --s1 ref/{wildcards.tissue1}.rmats.txt --s2 ref/{wildcards.tissue2}.rmats.txt  -t paired --readLength 130 --gtf ref/gencodeAno_bsc.gtf --bi ref/STARindex --od {output[0]}
        '''

rule find_tx_low_usage:
    input: expand('quant_files/{sampleID}', sampleID=sample_names)
    output:'tx_for_removal'
    shell:
        '''
        module load R
        Rscript scripts/soneson_low_usage.R
        '''

rule remove_tx_low_usage:
    input:'tx_for_removal'
    output: ref_trimmed
    run:
        with open(ref_fasta) as infasta, open('tx_for_removal') as bad_tx, open(ref_trimmed,'w+') as outfasta:
            names=set()
            for line in bad_tx:
                names.add('>'+line.strip())
            oldline=infasta.readline().strip().split('|')[0]
            while oldline:
                if oldline not in names and '>' in oldline:
                    write=True
                elif oldline in names and '>' in oldline:
                    write=False
                if write:
                    outfasta.write(oldline+'\n')
                oldline=infasta.readline().strip().strip().split('|')[0]

rule rebuild_salmon_index:
    input:ref_trimmed
    output:salmonindex_trimmed
    run:
        salmonindexcommand=loadSalmon + 'salmon index -t {} --gencode -i {} --type quasi -k 31'.format(ref_trimmed,salmonindex_trimmed)
        sp.run(salmonindexcommand, shell=True)

rule reQuantify_Salmon:
    input: lambda wildcards: ['fastq_files/{}_1.fastq.gz'.format(wildcards.sampleID),'fastq_files/{}_2.fastq.gz'.format(wildcards.sampleID)] if sample_dict[wildcards.sampleID]['paired'] else 'fastq_files/{}.fastq.gz'.format(wildcards.sampleID),
            salmonindex_trimmed
    output:'RE_quant_files/{sampleID}'
    log: 'logs/{sampleID}.rq.log'
    run:
        id=wildcards.sampleID
        paired=sample_dict[id]['paired']
        if paired:
            salmon_command=loadSalmon + 'salmon quant -i {} -l A --gcBias -1 {} -2 {} -o RE_quant_files/{}'.format(input[2],input[0],input[1],id)
        else:
            salmon_command=loadSalmon + 'salmon quant -i {} -l A --gcBias -r {} -o RE_quant_files/{}'.format(input[1],input[0],id)
        sp.run(salmon_command,shell=True)
        log1='logs/{}.rq.log'.format(id)
        salmon_info='RE_quant_files/{}/aux_info/meta_info.json'.format(id)
        if os.path.exists(salmon_info):
            with open(salmon_info) as file:
                salmonLog=json.load(file)
                mappingscore=salmonLog["percent_mapped"]
            if mappingscore <= 50:
                with open(log1,'w+') as logFile:
                    logFile.write('Sample {} failed QC mapping Percentage: {}'.format(id,mappingscore))
        else:
            with open(log1,'w+') as logFile:
                logFile.write('Sample {} failed to align'.format(id))
rule collapse_logs:
    input:expand('logs/{sampleID}.rq.log', sampleID=sample_names)
    output:'logs/failed','logs/bad_mapping'
    shell:
        '''
        for i in logs/*.rq* ; do cat $i >> logs/allLogs && echo '.' >> logs/allLogs ; done
        for i in logs/*.fqp ; do cat $i >> logs/allLogs && echo '.' >> logs/allLogs ; done
        grep '^.$' logs/allLogs -v > logs/tmp
        grep align logs/tmp > logs/failed
        grep map logs/tmp > logs/bad_mapping
        '''

# rule quality_control:
#     input:expand('RE_quant_files/{sampleID}',sampleID=sample_names),'logs/failed','logs/bad_mapping'
#     output:'smoothed_filtered_tpms.csv'
#     shell:
#         '''
#         #module load R
#         #Rscript QC.R
#         touch 'smoothed_filtered_tpms.csv'
#         '''


