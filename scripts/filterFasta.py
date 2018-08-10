with open('ref/gencodeRef.fa') as infasta, open('tx_for_removal') as bad_tx, open('ref/gencodeRef_trimmed.fa','w+') as outfasta:
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
