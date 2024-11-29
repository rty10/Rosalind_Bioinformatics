from Bio.Seq import Entrez

with open("rosalind_gbk.txt", "r") as f:
    sample = f.readlines()
    sample = [s.rstrip('\n') for s in sample]

Entrez.email = "robb.young.2011@gmail.com"
handle = Entrez.esearch(
    db='nucleotide', 
    term=sample[0]+'[Organism] AND '+sample[1]+':'+sample[2]+'[Publication Date]'
)
record = Entrez.read(handle)

print(int(record['Count']))
del record, handle, sample