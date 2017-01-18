#!/usr/bin/env python

import Bio
from Bio import SeqIO
from Bio import Entrez
from sys import argv

#def GetFASTAbyGI(GI):

script, GI, writeFile = argv

print ("This is the GI number: %s " % GI)
print ("This is the file to write too: %s " % writeFile)
print ("script: %s " % script)

entrezDbName = 'nucleotide'
Entrez.email = 'tod.p.stuber@usda.gov'

entryData = Entrez.efetch(db=entrezDbName, id=GI, retmode="text", rettype='fasta')

local_file=open(writeFile,"w")
local_file.write(entryData.read())
entryData.close()
local_file.close()

handle = open(writeFile, "r")
for record in SeqIO.parse(handle, "fasta"):
    print("FILE DOWNLOADED: *****", record.description,  "*****")
handle.close()
