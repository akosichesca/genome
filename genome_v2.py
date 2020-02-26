#!/usr/bin/env python
# coding: utf-8

# # Read the Reference Genome

# In[1]:


from tqdm import tqdm
import sys
import math
import numpy as np
nrow = 64
ncol = 64

fname = sys.argv[1] #"GCF_000005845.2_ASM584v2_genomic.fna"

f = open(fname,'r')
lines = f.readlines()[1:]
f.close()

genome = ''.join(lines)
print(genome[0:200])


# ## Convert the genome code to bits

# In[2]:


g = 'G'
gbitlen = 3

# Converts each genome base into bits. The choice of this encoding is based on three different principles:
# 1. Since we are using TCAM we want an equal hamming distance with each substitution
# 2. Since transitions occur (2.1x) more often than transversions. 
#    Transversions will have 2x more penalty than transitions.
#    Transitions have a hamming distance of 1 and transversions have a hamming distance of 2.
# 3. Genomes are double stranded the reads can come from any of the two strands. 
#    The other strand will be a reverse complement of the other strand.
#    The encoding helps in doing this, in such a way that the reverse of the encoding is already its complement.
def base2bits(g):
    code = '' #'****'
    #print(g)
    if g == 'A':
        code = '0*01'
    elif g == 'C':
        code = '11*0'
    elif g == 'T':
        code = '10*0'
    elif g == 'G':
        code = '0*11'
    return code

# Converts a genome string into a genome bit string
def genome2bits(genome):
    gbits = '' #'x' * (len(genome) * gbitlen)
    for g in genome: #range(len(genome):
        gbits = gbits + base2bits(g)
    return gbits    
    
gbits = genome2bits(genome)        


# # Read the reads of the genome to compare

# In[3]:


# Read the reads on a FastQ file
def readfile(fname):
    readsg = ''
    with open(fname) as fp:
        line = fp.readline()
        while line:
            line = fp.readline()
            readsg = readsg + line
            line = fp.readline()
            line = fp.readline()
            line = fp.readline()
    return readsg

# Encode the genome string of the reads
def get_reads(readsg):
    reads = ''
    for r in readsg.split('\n'):
        gr = genome2bits(r)
        reads = reads + gr + '\n'
    reads = reads.split('\n')[:-2]
    return reads


# In[4]:


# Compute for the hamming distance between two ternary strings
def hamming_distance(str1, str2):
    n = len(str1)
    m = len(str2)
    #print(n,m)
    hamm = 0
    for i in range(n):
        if not (str1[i] == str2[i] or str1[i] == '*' or str2[i] == '*'):
            hamm = hamm + 1
    return hamm


# In[5]:


# Mimics the TCAM hardware
class TCAM_Arrays:
    
    def __init__(self,nrow, ncolumn, error_threshold=4):
        self.nrow = nrow
        self.ncolumn = ncolumn
        self.array = []
        self.tag = [] 
        self.error_threshold = error_threshold
    
    # Write on multiple TCAM arrays
    # TODO: Create a function that writes on a specific row. (not necessary for genome sequencing)
    def write(self,bits):
        nrc = self.nrow*self.ncolumn
        m = len(bits) % nrc 
        bits = bits + "0"*(nrc -m)
        bitsrow = [bits[i:i+self.nrow] for i in range(0, len(bits), self.nrow)]
        self.array = [bitsrow[i:i+self.ncolumn] for i in range(0, len(bitsrow), self.ncolumn)]
        self.narray = len(self.array)
        self.tag = [[0] * nrow for i in range(self.narray)]
    
    def reset_tag(self):
        self.tag = [[0] * nrow for i in range(self.narray)]
        
    def get_tag(self):
        return self.tag
    
    # Implement a TCAM search of the query in all the arrays
    # Returns:
    #    pos  = positions (array, column, hamming distance=0) of exact matches
    #    pos2 = positions (array, column, hamming distance) of partial matches 
    #              with hamming distance less than the error_threshold (default = 4)
    def search(self,query):
        na = 0
        pos = []
        pos2 = []
        for a in self.array:
            nc = 0
            for ac in a:
                if ac != '*'*64:
                    e = hamming_distance(ac,query)
                    if e == 0:
                        pos.append([na, nc, e])
                        self.tag[na][nc] = True
                    elif e <= self.error_threshold:
                        pos2.append([na, nc, e])
                        self.tag[na][nc] = 1
                nc = nc + 1
            na = na + 1
        return pos, pos2

    # Find the array with the most tag
    def best_array(self):
        bestarray = [sum(array) for array in tcam.get_tag()]
        max_bestarray = np.argmax(bestarray)
        return max_bestarray
        
    # Find the rows with matches in the best_array
    def best_match(self):
        res = [rows for rows, val in enumerate(self.tag[self.best_array()]) if val > 0] 
        best_match = [indexi[j], -1, int(max_bestarray), [res]]
        return best_match


# In[6]:


# Save positional information
def save_pos(pos, pos_fname, prefix=''):
    with open(pos_fname, "a") as f:
        for s in pos:
            f.write(prefix + str(s[0]) + ' ' + str(s[1]) + ' ' + str(s[2]) +"\n")

def findreads(tcam, reads, indexi, pos_fname):
    
    pos = ''
    pos2 = ''
    ismatch = [False] * len(reads)    
    best_matches = []

    #for j in  readsi: 
    f = open(pos_fname + '.best','w')

    for j in  tqdm(range(len(indexi))):
        #print(j, " of ", len(indexi))
        tcam.reset_tag()
        idx = indexi[j]
        
        best_matches = []
        best_match = [0] * 2
            
        for reverse in [False, True]:
            matchflag = False
            if reverse == False:
                query = reads[idx] 
            else:
                query = reads[idx][::-1]
            N = len(query) - 64
            for i in tqdm(range(0,N,4)):
                q = query[i:i+64]
                #print(i, " of ", N, " shifts")
                pos, pos2 = tcam.search(q)
                #save_pos(pos,pos_fname, 'Match ' + str(indexi[j]) + ' ' + str(i) + ' ')
                #save_pos(pos2,pos_fname + '.par', 'Partial Match ' + str(indexi[j]) + ' ' + str(i) + ' ')
                #if len(pos) > 0:
                #    best_matches.append([indexi[j], i, -1, pos]) 
            ismatch[j] = matchflag

            res = [rows for rows, val in enumerate(tcam.get_tag()[tcam.best_array()]) if val > 0] 
            #best_matches.append([indexi[j], -1, int(tcam.best_array()), [res]])
            
            if reverse == False:
                best_match[0] = [int(tcam.best_array()), res]
            else: 
                best_match[1] = [int(tcam.best_array()), res]
                

        
        if len(best_match[0][1]) >= len(best_match[1][1]):
            position = (int(best_match[0][0])) *nrow*ncol/4 + (int(best_match[0][1][0]))*nrow/4 
        else:
            position = (int(best_match[1][0])) *nrow*ncol/4 + (int(best_match[1][1][0]))*nrow/4 
            
        best_matches.append([idx, position])
    
        for s in best_matches:
            f.write(str(s[0]) + ',' + str(s[1]))
            f.write('\n')
        #clear_output()
            
    f.close()
    return pos, pos2, best_matches


fname = sys.argv[2] #'bp100/bp100'
readsg = readfile('reads/' + fname + ".bfast.fastq")
reads = get_reads(readsg)
tcam = TCAM_Arrays(nrow,ncol)
tcam.write(gbits)
indexi = range(0,2)
pos_100, pos2_100, bestmatch = findreads(tcam, reads, indexi, 'results/' + fname + '.pos')
#pos_100, pos2_100, bestmatch = findreads(tcam, reads[0:20], indexi, 'results/' + fname + '.posr', reverse = True)


# In[ ]:


print(tagarray)

