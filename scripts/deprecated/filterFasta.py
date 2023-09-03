import sys
if len(sys.argv)==5:
    Gencode=True
with open(sys.argv[1]) as infasta, open(sys.argv[2]) as bad_tx, open(sys.argv[3],'w+') as outfasta:
    names=set()
    for line in bad_tx:
        names.add('>'+line.strip())
    oldline=infasta.readline().strip()
    if Gencode:
        oldline=oldline.split('|')[0]
    while oldline:
        if Gencode:
            oldline=oldline.split('|')[0]
        if oldline not in names and '>' in oldline:
            write=True
        elif oldline in names and '>' in oldline:
            write=False
        if write:
            outfasta.write(oldline+'\n')
        oldline=infasta.readline().strip()
