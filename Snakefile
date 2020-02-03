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
Things to do
-rewrite sonneson_low_usage
-rewrite script for removing_tx
'''
import subprocess as sp
import itertools as it

def readSampleFile(samplefile):
    # returns a dictionary of dictionaries where first dict key is sample id and second dict key are sample  properties
    res={}
    with open(samplefile) as file:
        for line in file:
            if line[0] == '#':
                continue
            info=line.strip('\n').split('\t')
            res[info[0]]={'files':info[1].split(','),'paired':True if info[2]=='y' else False, 'tissue':info[3],'subtissue':info[4]}
    return(res)

def lookupRunfromID(card,sample_dict):
    id=card
    if 'BLOOP_E-MTAB' in card: #not the best but it works
        return('bam_files/{}.bam'.format(id[:-2]))
    else:
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

#configfile:'config.yaml'
sample_dict=readSampleFile(config['sampleFile'])# sampleID:dict{path,paired,metadata}
sample_names=sample_dict.keys()

loadSRAtk="module load {} && ".format(config['sratoolkit_version'])
loadSalmon= "module load {} && ".format(config['salmon_version'])
salmonindex='ref/salmonindex'
salmonindex_trimmed='ref/salmonindex_trimmed'
STARindex='ref/STARindex'
ref_fasta='ref/gencodeRef.fa'
#ref_GTF='ref/gencodeAno.gtf'
ref_GTF_basic='ref/gencodeAno_bsc.gtf'
#ref_GTF_PA='ref/gencodeAno_pa.gtf'
#eref_PA='ref/gencodePA.fa'
badruns='badruns'
ref_trimmed='ref/gencodeRef_trimmed.fa'

rule all:
    input:
        expand('results/counts_{level}.csv.gz', level = ['gene', 'transcript']),
        'ref/salmonindexTrimmed',
        'results/word_clouds',
        config['EiaD_sqlite_file']
'''
****PART 1**** download files
'''
rule downloadGencode:
    output: ref_fasta,ref_GTF_basic
    shell:
        '''
        wget -O ref/gencodeRef.fa.gz {config[refFasta_url]}
        wget -O ref/gencodeAno_bsc.gtf.gz {config[refGTF_basic_url]}
        #wget -O ref/gencodePA.fa.gz {config[refPA_url]}
        gunzip ref/gencodeRef.fa.gz
        gunzip ref/gencodeAno_bsc.gtf.gz
        #gunzip ref/gencodePA.fa.gz
        '''

# set to 'NO' in config.yaml if you have fastq files
# you've added that aren't from SRA or you don't want
# to re-download them from SRA
# deposit fastq files in fastqParts/
if config['download_fastq'].upper() == 'YES':
    rule download_fastq:
        output: temp('fastqParts/{id}.fastq.gz')
        run:
            id=wildcards.id
            id=id[:-2] if '_'in id else id #idididid
            try:
                sp.check_output(loadSRAtk + 'fastq-dump --gzip --split-3 -O fastqParts {}'.format(id),shell=True)
            except sp.CalledProcessError:
                with open('logs/{}.fqp'.format(wildcards.id)) as l:
                    l.write('{} did not download'.format(wildcards.id))

# combine multiple fastq files together
rule aggregate_fastqs:
    input:ancient(lambda wildcards:lookupRunfromID(wildcards.sampleID,sample_dict))
    output:'fastq_files/{sampleID}.fastq.gz'
    run:
        #this can use some cleaning up - rule runs twice for paired
        id=wildcards.sampleID
        if '.bam' in input[0]:
            id=wildcards.sampleID[:-2]
            #need to collate a bam before you can convert, otherwise will lose many reads
            cmd='module load samtools &&  samtools collate -O {} | \
            samtools fastq -1 fastq_files/{}_1.fastq -2 fastq_files/{}_2.fastq -0 /dev/null -s /dev/null -n -F 0x900 -'.format(input[0],id,id)
            sp.run(cmd,shell=True)
            gunzip=' gunzip -c -f fastq_files/{}_1.fastq > fastq_files/{}_1.fastq.gz'.format(id,id)
            sp.run(gunzip,shell=True)
            gunzip=' gunzip -c -f fastq_files/{}_2.fastq > fastq_files/{}_2.fastq.gz'.format(id, id)
            sp.run(gunzip,shell=True)
        else:
            fileParts=lookupRunfromID(id,sample_dict)
            i='1' if '_' in id and id[-1]=='1' else '2'# which strand
            id=id[:-2] if '_' in id else id
            for fqp in fileParts:
                if sample_dict[id]['paired']:
                    sp.run('cat {fqp} >> fastq_files/{id}_{i}.fastq.gz '.format(fqp=fqp,i=i,id=id),shell=True)
                else:
                    sp.run('cat {fqp} >> fastq_files/{id}.fastq.gz'.format(fqp=fqp,id=id),shell=True)

# if set to yes in config.yaml, then salmon will be run twice
# after the first salmon run all transcripts across all samples
# will be evaluated and those transcripts with extremely low
# counts will be removed and the salmon index is rebuilt,
# and salmon is run again on the new index
if config['build_new_salmon_index'].upper() == 'YES':
    '''
    ****PART 2*** Initial quantification
    '''
    rule build_salmon_index:
        input:  ref_fasta
        output:'ref/salmonindex'
        run:
            salmonindexcommand=loadSalmon + 'salmon index -t {} --gencode -i {} --type quasi --perfectHash -k 31'.format(input[0],output[0])
            sp.run(salmonindexcommand, shell=True)



    rule run_salmon:
        input: ancient(lambda wildcards: ['fastq_files/{}_1.fastq.gz'.format(wildcards.sampleID),'fastq_files/{}_2.fastq.gz'.format(wildcards.sampleID)] if sample_dict[wildcards.sampleID]['paired'] else 'fastq_files/{}.fastq.gz'.format(wildcards.sampleID)),
            'ref/salmonindex'
        output: 'quant_files/{sampleID}/quant.sf'
        log: 'logs/{sampleID}.log'
        run:
            id=wildcards.sampleID
            #tissue=wildcards.tissue
            paired=sample_dict[id]['paired']
            if paired:
                salmon_command=loadSalmon + 'salmon quant -i {} -l A --validateMappings --gcBias --seqBias -p 4  -1 {} -2 {} -o {}'.format(input[2],input[0],input[1],'quant_files/{}'.format(id))
            else:
                salmon_command=loadSalmon + 'salmon quant -i {} -l A --validateMappings --gcBias --seqBias -p 4 -r {} -o {}'.format(input[1],input[0],'quant_files/{}'.format(id))
            sp.run(salmon_command,shell=True)
            log1='logs/{}.log'.format(id)
            salmon_info='quant_files/{}/aux_info/meta_info.json'.format(id)
            if os.path.exists(salmon_info):
                with open(salmon_info) as file:
                    salmonLog=json.load(file)
                    mappingscore=salmonLog["percent_mapped"]
                if mappingscore <= 35:
                    with open(log1,'w+') as logFile:
                        logFile.write('Sample {} failed QC mapping Percentage: {}'.format(id,mappingscore))
            else:
                with open(log1,'w+') as logFile:
                    logFile.write('Sample {} failed to align'.format(id))

    '''
    ****PART 3**** find and remove lowly used transcripts
    '''
    #problem with tximport; going to have to make a custom tsxdb from all the gtfs

    rule find_tx_low_usage:
        input: expand('quant_files/{sampleID}/quant.sf', sampleID=sample_names), 'ref/gencodeAno_bsc.gtf'
        output:'tx_for_removal.txt'
        shell:
            '''
            module load R
            Rscript {config[scripts_dir]}/soneson_low_usage.R {config[working_dir]} {ref_GTF_basic}
            '''

    rule remove_tx_low_usage:
        input:'tx_for_removal.txt', ref_fasta
        output: 'ref/gencodeRef_trimmed.fa'
        shell:
            '''
            python3 {config[scripts_dir]}/filterFasta.py {input[1]} {input[0]} {output[0]} Gencode
            '''


    '''
    ***PART 4*** requantify salmon
    Remove transcripts with no or very low usage (Sonenson method)
    '''

    rule rebuild_salmon_index:
        input:'ref/gencodeRef_trimmed.fa'
        output:'ref/salmonindexTrimmed'
        run:
            salmonindexcommand=loadSalmon + 'salmon index -t {} --gencode -i {} --type quasi --perfectHash -k 31'.format(input[0],output[0])
            sp.run(salmonindexcommand, shell=True)


    rule re_run_Salmon:
        input: lambda wildcards: ['fastq_files/{}_1.fastq.gz'.format(wildcards.sampleID),'fastq_files/{}_2.fastq.gz'.format(wildcards.sampleID)] if sample_dict[wildcards.sampleID]['paired'] else 'fastq_files/{}.fastq.gz'.format(wildcards.sampleID),
            'ref/salmonindexTrimmed'
        output:'RE_quant_files/{sampleID}/quant.sf'
        threads: 4
        log: 'logs/{sampleID}.rq.log'
        run:
            id=wildcards.sampleID
            #tissue=wildcards.tissue
            paired=sample_dict[id]['paired']
            if paired:
                salmon_command=loadSalmon + 'salmon quant -i {} -l A --validateMappings --gcBias --seqBias -p 8 -1 {} -2 {} -o {}'.format(input[2],input[0],input[1],'RE_quant_files/{}'.format(id))
            else:
                salmon_command=loadSalmon + 'salmon quant -i {} -l A --validateMappings --gcBias --seqBias -p 8 -r {} -o {}'.format(input[1],input[0],'RE_quant_files/{}'.format(id))
            sp.run(salmon_command,shell=True)
            log1='logs/{}.rq.log'.format(id)
            salmon_info='RE_quant_files/{}/aux_info/meta_info.json'.format( id)
            if os.path.exists(salmon_info):
                with open(salmon_info) as file:
                    salmonLog=json.load(file)
                    mappingscore=salmonLog["percent_mapped"]
                if mappingscore <= 35:
                    with open(log1,'w+') as logFile:
                        logFile.write('Sample {} failed QC mapping Percentage: {}'.format(id,mappingscore))
                else:
                    with open(log1,'w+') as logFile:
                        logFile.write('Sample {} passed mapping QC')
            else:
                with open(log1,'w+') as logFile:
                    logFile.write('Sample {} failed to align'.format(id))

# if set to no, then all samples will be run on the user given salmon index path
# which we imagine will be the latest trimmed salmon index from rebuild_salmon_index
# and/or the salmon_quant files downloaded from eyeIntegration
# the idea is not to do the fairly computationally expensive step of building salmon
# indices twice and running salmon twice, as adding a few new samples is not
# going to alter to the salon index much
else:
    rule re_run_Salmon:
        input: lambda wildcards: ['fastq_files/{}_1.fastq.gz'.format(wildcards.sampleID),'fastq_files/{}_2.fastq.gz'.format(wildcards.sampleID)] if sample_dict[wildcards.sampleID]['paired'] else 'fastq_files/{}.fastq.gz'.format(wildcards.sampleID),
            config['salmon_index_path']
        output:'RE_quant_files/{sampleID}/quant.sf'
        threads: 4
        log: 'logs/{sampleID}.rq.log'
        run:
            id=wildcards.sampleID
            #tissue=wildcards.tissue
            paired=sample_dict[id]['paired']
            if paired:
                salmon_command=loadSalmon + 'salmon quant -i {} -l A --validateMappings --gcBias --seqBias -p 8 -1 {} -2 {} -o {}'.format(input[2],input[0],input[1],'RE_quant_files/{}'.format(id))
            else:
                salmon_command=loadSalmon + 'salmon quant -i {} -l A --validateMappings --gcBias --seqBias -p 8 -r {} -o {}'.format(input[1],input[0],'RE_quant_files/{}'.format(id))
            sp.run(salmon_command,shell=True)
            log1='logs/{}.rq.log'.format(id)
            salmon_info='RE_quant_files/{}/aux_info/meta_info.json'.format(id)
            if os.path.exists(salmon_info):
                with open(salmon_info) as file:
                    salmonLog=json.load(file)
                    mappingscore=salmonLog["percent_mapped"]
                if mappingscore <= 35:
                    with open(log1,'w+') as logFile:
                        logFile.write('Sample {} failed QC mapping Percentage: {}'.format(id,mappingscore))
                else:
                    with open(log1,'w+') as logFile:
                        logFile.write('Sample {} passed mapping QC')
            else:
                with open(log1,'w+') as logFile:
                    logFile.write('Sample {} failed to align'.format(id))

rule process_poor_mapped_samples:
    input: expand('logs/{sampleID}.rq.log',sampleID=sample_names)
    output: 'ref/bad_mapping.txt'
    shell:
        '''
        for i in logs/*.rq.log; do cat $i && echo  ; done | grep failed - > {output}
        '''

localrules:extract_mapping_rate
rule extract_mapping_rate:
    input:
        expand('RE_quant_files/{sampleID}/quant.sf', sampleID=sample_names)
    output: 'results/mapping_rates.txt'
    shell:
        '''
        for i in RE_quant_files/*/logs/salmon_quant.log; 
            do echo "$i" | cut -f2 -d"/" | xargs printf ;
            printf " "; 
            grep "Mapping rate" $i | tail -n 1 | awk '{{print $NF}}'; 
        done > {output}
        '''

'''
***** Produce files for eyeIntegration web app
'''
# load all salmon quant with with tximport (lengthScaledTPM)
# run QSmooth normalization
# remove genes with near 0 counts for all samples
# remove samples with:
# 	low counts
#	normalize lengthScaled TPM by library size
#	run tsne, cluster with DBScan, remove samples that are more than 4 SD from cluster center
rule gene_quantification_and_normalization:
    input:
        tpms=expand('RE_quant_files/{sampleID}/quant.sf',sampleID=sample_names),
        gtf='ref/gencodeAno_bsc.gtf',
        bad_map='ref/bad_mapping.txt',
        mapping_rate='results/mapping_rates.txt'
    params:
        working_dir = config['working_dir']
    output:
        counts = 'results/counts_{level}.csv.gz',
        tpm = 'results/smoothed_filtered_tpms_{level}.csv',
        removed_samples = 'results/samples_removed_by_QC_{level}.tsv',
        cor_scores = 'results/cor_scores_{level}.tsv',
        batchCor_tpm = 'results/smoothed_filtered_tpms_batchCor_{level}.csv'
    shell:
        '''
        module load R
        Rscript {config[scripts_dir]}/QC.R {config[sampleFile]} {ref_GTF_basic} {params.working_dir} {wildcards.level} {input.bad_map} {input.mapping_rate} {output}
        '''

# output sample metadata and gene/tx lists for eyeIntegration
rule make_meta_info:
    input:
        expand('results/smoothed_filtered_tpms_{level}.csv',level = ['gene','transcript']),
        expand('results/samples_removed_by_QC_{level}.tsv', level = ['gene','transcript']),
        'results/mapping_rates.txt'
    params:
        working_dir = config['working_dir']
    output:
        metadata = 'results/core_tight.Rdata',
        tx_names = 'results/tx_names.Rdata',
        gene_names = 'results/gene_names.Rdata',
        gene_tx_info = 'results/gene_tx_gtf_info.Rdata'
    shell:
        '''
        module load R
        Rscript {config[scripts_dir]}/make_meta_info.R \
          {config[sampleFile]} \
          {ref_GTF_basic} \
          {config[sqlfile]} \
          {input} \
          {params.working_dir} \
          {output} \
          {config[scripts_dir]}
        '''

rule differential_expression:
    input: 
        'results/smoothed_filtered_tpms_{level}.csv',
        'results/core_tight.Rdata'
    params:
        working_dir = config['working_dir'], #'/data/swamyvs/autoRNAseq'
    output:
        comparisons = 'results/de_comparison_name_list_{level}.Rdata',
        limma_object = 'results/limma_DE_object_{level}.Rdata',
        list_of_dataframes = 'results/limma_DE_listDF_{level}.Rdata'
    shell:
        '''
        module load R
        Rscript {config[scripts_dir]}/diffExp.R {params.working_dir} {config[sampleFile]} {input} {output}
        '''

# for each gene/TX, by sub_tissue, calculate mean expression, rank, and decile
rule calculate_mean_rank_decile:
    input:
        lsTPM_file = 'results/smoothed_filtered_tpms_batchCor_{level}.csv',
        metadata_file = 'results/core_tight.Rdata',
    params:
        working_dir = config['working_dir']
    output:
        'results/mean_rank_decile_{level}.tsv'
    shell:
        '''
        module load R
        Rscript {config[scripts_dir]}/calculate_mean_rank_decile.R \
          {params.working_dir} \
          {input.lsTPM_file} \
          {input.metadata_file} \
          {output}
        '''

# build differential gene lists for all comparisons done in
# differential expression for GO enrichment
rule differential_gene_lists:
    input:
        limma_object = 'results/limma_DE_object_gene.Rdata',
    params:
        working_dir = config['working_dir']
    output:
        up = 'results/up_gene_comparison_lists.Rdata',
        down = 'results/down_gene_comparison_lists.Rdata',
        all_genes = 'results/all_genes.Rdata'
    shell:
        '''
        module load R
        Rscript {config[scripts_dir]}/find_diff_expressed_gene_lists.R \
          {params.working_dir} \
          {input} \
          {output}
        '''

# calculate GO term enrichments
rule GO_term_enrichment:
    input:
        limma_object = 'results/limma_DE_object_gene.Rdata'
    params:
        working_dir = config['working_dir']
    threads: 30
    output:
        all_vs_all_go = 'results/all_vs_all_GO.Rdata'
    shell:
       '''
       module load R
       Rscript {config[scripts_dir]}/calculate_GO_enrichment.R \
         {params.working_dir} \
         {threads} \
         {input} \
         {output}
       '''

# calculate t-SNE coordinates
rule tSNE:
    input:
        metadata = 'results/core_tight.Rdata',
        tpm = 'results/smoothed_filtered_tpms_batchCor_gene.csv',
        gtf = 'results/gene_tx_gtf_info.Rdata'
    params:
        working_dir = config['working_dir']
    output:
        'results/tSNE_coords.Rdata'
    shell:
        '''
        module load R
        Rscript {config[scripts_dir]}/calculate_tsne_5_to_50.R \
          {params.working_dir} \
          {input} \
          {output}
        '''

# make word cloud png for GO enrichment 
rule make_word_clouds:
    input:
        'results/all_vs_all_GO.Rdata'
    output:
        'results/word_clouds'
    params:
        working_dir = config['working_dir']
    shell:
        """
        module load {config[R_version]}
        mkdir -p {output}
        Rscript {config[scripts_dir]}/make_word_cloud_pngs.R \
          {params.working_dir} \
          {input} \
          {output}
        """

# create SQLite expression db
rule make_SQLite_db:
    input:
        tpms = expand('results/smoothed_filtered_tpms_batchCor_{level}.csv', level = ['gene', 'transcript']),
        tx_names = 'results/tx_names.Rdata',
        DE = expand('results/limma_DE_listDF_{level}.Rdata', level = ['gene', 'transcript']),
        DE_tests = 'results/de_comparison_name_list_gene.Rdata',
        GO = 'results/all_vs_all_GO.Rdata',
        mrd = expand('results/mean_rank_decile_{level}.tsv', level = ['gene', 'transcript']),
        metadata = 'results/core_tight.Rdata',
        gene_info = 'results/gene_tx_gtf_info.Rdata',
        tSNE = 'results/tSNE_coords.Rdata'
    params:
        working_dir = config['working_dir']
    output:
        config['EiaD_sqlite_file']
    shell:
       '''
       module load R
       Rscript {config[scripts_dir]}/make_sqlite_db.R \
         {params.working_dir} \
         {input} \
         {output}
       '''
