with open('infile') as infasta, open('filter') as crit, open('out.fa','w+') as outfasta:
    names=[line.strip() for line in crit]
    oldline=infasta.readline().strip()
    while oldline:
        if oldline not in names and '>' in oldline:
            write=True
        elif oldline in names and '>' in oldline:
            write=False
        if write:
            outfasta.write(oldline+'\n')
        oldline=infasta.readline().strip()
