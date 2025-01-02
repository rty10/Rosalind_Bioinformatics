# 


# --- Dependencies
import os, random, math

# -!- File loading, reading, and writing
def sequence_from_file(path_to_file) -> str:
    '''
    Load file with only one sequence
    return sequence
    '''
    with open(path_to_file, "r") as infile:
        return infile.read().rstrip('\n')


def fasta_dictionary(path_to_filename):
    '''
    A more robust function that in the `nb_gc` notebook. Here, open a FASTA file and keep only the identifier and the sequence (ignoring any additional information).
    Convert into a dictionary.
    Output: a dictionary where the key-value pairs are sequence IDs and sequences, respectively
    '''
    with open(path_to_filename, 'r') as f:
        lst  = f.readlines()
    f.close()
    for i in range(len(lst)):
        if lst[i].startswith('>'):
            lst[i] = lst[i].split(' ')[0]+'\n'
    lst  = [i.replace('\n', ' ') for i in lst]
    str1 = ''.join(lst)  
    lst2 = str1.split('>')
    lst2 = lst2[1:]
    seq_dict = {lst2[i].split(' ')[0]:''.join(lst2[i].split(' ')[1:]) for i in range(len(lst2))}
    del lst, lst2
    return seq_dict
# ---------------------------------------------------------------------------


# -!- Loading standard informatio (e.g. RNA codons, amino acid residues)
def mrna_codon_table(path_to_standards):
    '''
    Load standardf file that contains all unique 3-char mRNA codons for its AA-residue
    Return dictionary- key:value = RNA:AA pairings, including Stop codons
    '''
    with open(path_to_standards + "/rosalind_RNA_codon_table.txt", "r") as infile:
        rna_table      = infile.read().replace('\n', ' ').split()
        rna_codon_dict = {rna_table[i]:rna_table[i+1] for i in range(len(rna_table)) if i%2==0}
        del rna_table
    return rna_codon_dict


def protein_molar_mass(path_to_standards, aa_residue_seq: str) -> float:
    '''
    '''
    mm_water = 18.01056
    with open( path_to_standards + '/rosalind_monoisotopic_mass_table.txt', 'r') as infile:
        data = infile.read().replace('\n',' ').split()
        mass_dict = {data[i]:float(data[i+1]) for i in range(len(data)) if i%2==0}
        del data
    mm = 0
    for a in aa_residue_seq:
        mm += mass_dict[a]
    return round(mm,3)


def mrna_codon_frequency_table(dict_mrna_codon_table):
    '''
    Load dictionary that has each unique (3-char mRNA):(1-char AA-residue) pairing 
    Return frequency dictionary- key:value = AA-residue:freq, including Stop codons
    '''
    rna_codon_counter_dict = {i:0 for i in sorted(['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','Stop'])}
    for codon, counter in rna_codon_counter_dict.items():
        ct = 0
        for key, val in dict_mrna_codon_table.items():
            if val == codon:
                ct += 1
        rna_codon_counter_dict[codon] = ct
        del ct
    return rna_codon_counter_dict
# ---------------------------------------------------------------------------


# -!- General sequence details: counting, transcription, translation
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


def dna_rna_transcribe(dna_string: str) -> str:
    '''
    -!- From Rosalind.info/Stronghold/rna
    Load in a string of DNA nucleotides and transcribe into RNA nucleotides,
       replacing thymine 'T' with uracil 'U'
    Result will be a output of the transcribed RNA string text
    '''
    return ''.join( [i for i in dna_string] ).replace('T','U')


def rna_prot_translate(rna_seq: str) -> str:
    '''
    For the RNA nucleotide string, translate into an amino acid (AA) residue string.
    Step 1: load RNA codon table
        Using data found at https://rosalind.info/glossary/rna-codon-table/
        Result in a dictionary where the 3-nt RNA codon is the key and the 1-aa residue is the value 
    Step 2: begin translation 
        Make sure that is the RNA codon translate 'Stop' that one does not include 'Stop' in the residue string
    Result: string of AA residues
    '''
    with open(os.getcwd() + "/datasets/rosalind_RNA_codon_table.txt", "r") as f:
        text = f.read().replace('\n', ' ').split()
    codon_dict = {text[i]:text[i+1] for i in range(len(text)) if i%2==0}
    del text
    aa_seq = ''
    for i in range(len(rna_seq)):
        if (i%3 == 0):
            codon = rna_seq[i:i+3]
            if len(codon) == 3:
                if codon_dict[codon] == 'Stop':
                    break
                else:
                    aa_seq += codon_dict[ codon ]
            del codon
    del codon_dict
    return aa_seq
# ---------------------------------------------------------------------------


# -!- Functions relating to Graph Theory
def overlap_graph(sequence_dict, overlap_size:int):
    '''
    Load in a sequence dictionary from a FASTA-formatted file.
    Extract a suffix of size __overlap_size__ from one sequence and compare to all prefixes of the remaining sequences.
    Output: Print statement of the two keys in suffix-prefix order.
      Optional: this can be converted such that the output can be an adjacency list as a list-of-tuples
    '''
    adj_list = []
    sd_keys = list(sequence_dict.keys())
    for k1 in sd_keys:
        suff   = sequence_dict[k1][-overlap_size:]
        n_keys = [i for i in sd_keys if i != k1]
        for k2 in n_keys:
            if suff == sequence_dict[k2][:overlap_size]:
                adj_list.append(k1+' '+k2)
        del suff, n_keys
    del sd_keys
    for i in adj_list:
        print(i)
    del adj_list
    return
# ---------------------------------------------------------------------------


# -!- Functions relating to Motifs, Subsequences, and Sequence Comparables
def common_motif_checker(path_to_fasta):
    '''
    From Rosalind.info/stronghold/lcsm
    Load in the path to FASTA document and generate a sequence dictionary.
    First, determine all unique substrings found in every sequence of document, ranging in size from 2-nt to
     a maximum length of L-1, where L is the length of that sequence.
    Second, for each item in the unique substrings list, determine if found in every sequence in dictionary.
    Third, determine the largest length of substing in the common substrings list, and then keep substrings with only that length.
    Return: list of substrings of longest length
    '''
    seq_dict        = fasta_dictionary(path_to_fasta)
    full_substrings   = []
    for key in list(seq_dict.keys()):
        check_list=[]
        maxlength = len(seq_dict[key])
        if maxlength >= 250:
            maxlength = 250
        for i in range(2, maxlength):
            for j in range(0, len(seq_dict[key])):
                substring = seq_dict[key][j:j+i]
                if len(substring)==i:
                    check_list.append(substring)
                del substring
        check_list = list(set(check_list))
        if len(check_list) > len(full_substrings):
            full_substrings = check_list[:]
            std_key = key
        del maxlength, check_list
    print("All possible uniques substrings between n=2 & n=(L-1) with standard key "+std_key+":", len(full_substrings))
    common_substrings    = []
    max_substring_length = 2
    for substring in full_substrings:
        if len(substring) >= max_substring_length:
            counter_seqs = 1
            for key in [K for K in list(seq_dict.keys()) if K != std_key]:
                if substring in seq_dict[key]:
                    counter_seqs += 1
            if counter_seqs == (len(list(seq_dict.keys()))):
                common_substrings.append(substring)
                max_substring_length = len(substring)
            del counter_seqs     
    print("All possible uniques and common substrings in file:", len(common_substrings))
    # Prune substrings list, keeping only those who share the highest substring length
    max_substring_length = 0
    for item in common_substrings:
        if len(item) > max_substring_length:
            max_substring_length = len(item)  
    common_substrings = [i for i in common_substrings
                         if len(i) == max_substring_length]
    print("All possible common substrings at maximum substring length n="+str(max_substring_length)+":", len(common_substrings))
    del max_substring_length, full_substrings, std_key
    return common_substrings
# ---------------------------------------------------------------------------


# -!- Functions relating to Combinatorics
def permutation_integer_list(n: int):
    '''
    For some integer _n_, return a list of all possible strings using each of the numbers from 1 to _n_,
        without repeats.
    '''
    sequence = ''.join( [str(i+1) for i in range(n)] )
    perm_ct  = math.factorial(n)
    perm_lst = []
    val      = 0
    while val < perm_ct:
        gentext = ''.join( random.choice(sequence) for _ in range(n) )
        ct = 0
        for g in gentext:
            if gentext.count(g) == 1:
                ct += 1
        if (ct == n)&(' '.join([g for g in gentext]) not in perm_lst):
            perm_lst.append(' '.join([g for g in gentext]))
            val += 1
        del gentext, ct
    del val, perm_ct, sequence
    return sorted(perm_lst)


def protein_mrna_combinations(protein_sequence: str,
                              modulo: int,
                              dict_mrna_codon_frequency,
                              verbose=False) -> int:
    '''
    From Rosalind.info/stronghold/mrna
    Load in the protein sequence string, the modulo value for modular arithmetic, and the standard mRNA codon frequency dictionary
    Return is __verbose__-dependent:
        - default: return the protein_sequence_combination mod __modulo__ integer value
        - True: print outs of dictionaries, and the modular arithmetic value through each protein_sequence_combination iteration
    '''
    # Count frequency of each residue in sequence, compile into dictionary, then insert 'Stop':1    
    aa_res_counter = {i:0 for i in sorted(['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'])}

    for i in range(len(protein_sequence)):
        aa_res_counter[protein_sequence[i]] += 1
    
    # -!- to save space, remove anything where the value is 0
    aa_res_counter = {key:aa_res_counter[key]
                      for key in aa_res_counter.keys()
                      if aa_res_counter[key] != 0}
    
    # Compile dictionary of frequencies*(number of mRNA codon strings)
    mrna_seqs_counter = {i:0 for i in sorted(list(aa_res_counter.keys()))}
    mrna_seqs_counter['Stop'] = 3
    for key, val in aa_res_counter.items():
        mrna_seqs_counter[key] =  dict_mrna_codon_frequency[key] ** aa_res_counter[key]
    
    if verbose == True:
        print(str(len(protein_sequence))+' residues in length')
        print(protein_sequence)
        print()
        print(''.join(sorted(protein_sequence)))
        print()
        print("Frequency of each residue in sequence:"+'\t', aa_res_counter)
        print("Count of mRNA codons/aa-residue:"+'\t', {key:dict_mrna_codon_frequency[key]
                      for key in aa_res_counter.keys()
                      if aa_res_counter[key] != 0})
        print("Frequency of potential mRNA codons in sequence:"+'\t', mrna_seqs_counter)
        print()

    # Iteratively multiply the values in the previous dictionary to get total number of mRNA string combinations
    #rna_cts = 1
    #for key, val in mrna_seqs_counter.items():
    #    rna_cts *= val
    #    if verbose == True:
    #        print(str(rna_cts)+" mod "+str(modulo)+" =", str( divmod(rna_cts, modulo) ))

    # -!- try this for improved memory storage
    rna_cts = 1
    for key, val in mrna_seqs_counter.items():
         if verbose == True:
            print("("+str(rna_cts)+" x "+str(val)+") mod "+str(modulo)+" =", str( (rna_cts*val)%modulo ))
         rna_cts = (rna_cts*val)%modulo
            
    del mrna_seqs_counter, aa_res_counter
    
    # Return the rna_cts mod modulo value
    return (rna_cts % modulo)


