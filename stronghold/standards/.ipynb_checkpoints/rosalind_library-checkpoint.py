# 

def dna_nt_count(string: str):
    '''
    -!- From Rosalind.info/Stronghold/dna
    Load in a string of DNA nucleotides and counts the number of A, C, G, and T bases.
    Returns a TUPLE of order ('A', 'C', 'G', 'T')
    '''
    a, t, g, c = 0,0,0,0
    for i in string:
        if i == 'A':
            a+=1
        elif i == 'C':
            c+=1
        elif i == 'G':
            g+=1
        elif i == 'T':
            t+=1
    return (a, c, g, t)


def dna_rna_transcribe(dna_string: str) -> str:
    '''
    -!- From Rosalind.info/Stronghold/rna
    Load in a string of DNA nucleotides and transcribe into RNA nucleotides,
       replacing thymine 'T' with uracil 'U'
    Result will be a output of the transcribed RNA string text
    '''
    return ''.join( [i for i in dna_string] ).replace('T','U')


def dna_complement(string: str) -> str:
    '''
    -!- From Rosalind.info/Stronghold/revc
    Load in string text of DNA nucleotides ('coding strand').
    Return a string of complementary bases ('complementary strand').
    The string text is read 5` to 3`, so the complement will also be 5`-3`.
    Notes: 
       A-T, T-A, C-G, and G-C
    '''
    nt_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    return ''.join( [nt_dict[i] for i in string] )


