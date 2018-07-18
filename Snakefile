import os
import subprocess as sp

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

#open file with all ids
SRA_IDS=[]
with open(config['ids']) as file:
    for line in file:
        SRA_IDS.append(str(line).strip('"|\n'))

#define target files
rule all:
    input: expand("{sraids}/quant.sf",sraids=SRA_IDS)


# #download fastqs
# for id in SRA_IDS:
#     sp.run("fastq-dump -Z {}  > {}.fastq ".format(id,id),shell=True)



#this works
salmonindex='salmon.index'
# if  not os.path.exists(salmonindex):
#         sp.run('wget -nc ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/gencode.v26.transcripts.fa.gz')
#         salmonindexcommand='salmon index -t gencode.v26.transcripts.fa.gz -i' +salmonindex+ ' --type quasi -k 31'
#         sp.run(salmonindexcommand, shell=True)


#run salmon on 
rule run_salmon:
        input:"{sraid}.fastq"
        output:"{sraid}/quant.sf"
        run:
            print(wildcards.sraid)
            pairedFlag=isPaired(wildcards.sraid)
            salmon_command='salmon quant -i ' + salmonindex + ' -l A '+ pairedFlag+ ' {input} -o '+ wildcards.sraid
            sp.run(salmon_command,shell=True)  