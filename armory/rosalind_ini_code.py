from Bio.Seq import Seq

with open("rosalind_ini.txt", "r") as f:
    sample = f.read().rstrip('\n')

sequence = Seq(sample)

print(sequence.count('A'),
      sequence.count('C'),
      sequence.count('G'),
      sequence.count('T'))

del sequence, sample