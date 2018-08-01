import subprocess as sp
sample_table={}
bad_sample=[]
with open('srsids.txt') as ids , open('sampleTable.txt','w+') as outfile, open('badsample.txt','w+') as bs:
    for id in ids:
        id=id.strip('\n')
        command="sqlite3 ref/SRAmetadb_072418.sqlite \"select sample_accession,run_accession,library_layout from sra WHERE sample_accession=\'{}\' \"".format(id)
        sampleInfo=sp.check_output(command,shell=True).decode('utf-8').strip("\n")
        sampleInfo=sampleInfo.split('\n')
        sample_table[id]=[id, '','y']
        try:
            for entry in sampleInfo:            
                info=entry.split('|')
                sample_table[id][1]+=info[1]+','
                if 'PAIRED' not in info[2]:
                    sample_table[id][2]='n'
        except IndexError:
            bs.write( id+ '\n')
        outfile.write(sample_table[id][0] +'\t'+ sample_table[id][1].strip(',') +'\t'+ sample_table[id][2] +'\n')
