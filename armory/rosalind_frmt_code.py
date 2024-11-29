from Bio.Seq import Entrez, SeqIO

Entrez.email = "robb.young.2011@gmail.com"

with open("rosalind_frmt.txt", "r") as f:
    sample = f.read().rstrip('\n')

handle = Entrez.efetch(
    db='nucleotide', 
    id=[sample],
    rettype='fasta'
)
records = list(SeqIO.parse(handle, "fasta"))

min = 1e7
min_label = ''
for i in range(len(sample.split())):
    if len(records[i].seq) < min:
        min = len(records[i].seq)
        min_label = records[i].id
del handle, records, min


handle = Entrez.efetch(
    db='nucleotide', 
    id=[min_label],
    rettype='fasta'
)
records = handle.read()
print(records)
del records, handle, min_label, sample