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
genome_start = int(sys.argv[2]) - 1
array_length = int(sys.argv[3])

f = open(fname,'r')

genome = []



with open(fname) as fp:
  fp.seek(genome_start*(int(nrow*ncol/4)+1))
  l = fp.readline()
#  for j in range(genome_start):
#    l = fp.readline()
  for j in range(array_length):
    #print(l)
    l = fp.readline()
    genome.append(l.strip())
    if not l:
      break
      
genome = ''.join(genome)

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


# In[3]:


#print(genome[0:10])

def base2mask(g):
    code = '' #'****'
    mask = ''
    #print(g)
    if g == 'A':
        code = '0001'
        mask = '0100'
    elif g == 'C':
        code = '1100'
        mask = '0010'
    elif g == 'T':
        code = '1000'
        mask = '0010'
    elif g == 'G':
        code = '0011'
        mask = '0100'
    elif g == 'N':
        code = '0000'
        mask = '1111'
    return code, mask

def genome2mask(genome):
    code = '' #'x' * (len(genome) * gbitlen)
    mask = ''
    codeint = []
    maskint = []
        
    if len(genome)%16 != 0:
        genome = genome + 'N'*(16-len(genome)%16)
    
    for g in genome: #range(len(genome):
        c,m = base2mask(g)
        code = code + c
        mask = mask + m
        
    #print(code,mask)
    for i in range(0,len(code),ncol):
        codeint.append(int(code[i:i+ncol],2))
        maskint.append(int(mask[i:i+ncol],2))
    return codeint, maskint

genc,genm = genome2mask(genome)


# # Read the reads of the genome to compare

# In[4]:


# Read the reads on a FastQ file
def readfile(fname, start, length):
    readsg = ''
    #print(start,length)
    with open(fname) as fp:
        line = fp.readline()
        for j in range(start):
            line = fp.readline()
            line = fp.readline()
            line = fp.readline()
            line = fp.readline()
        for j in range(length):
            line = fp.readline()
            readsg = readsg + line
            line = fp.readline()
            line = fp.readline()
            line = fp.readline()
            if not line:
                break
    return readsg

# Encode the genome string of the reads
def get_reads(readsg):
    reads = ''
    for r in readsg.split('\n'):
        #gr = genome2bits(r)
        reads = reads + r + '\n'
    reads = reads.split('\n')[:-2]
    return reads

# In[5]:


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


# In[6]:


# Mimics the TCAM hardware
class TCAM_Arrays:
    
    def __init__(self,nrow, ncolumn, error_threshold=4):
        self.nrow = nrow
        self.ncolumn = ncolumn
        self.narray = 0
        self.array = []
        self.array_code = []
        self.array_mask = []
        self.tag = [] 
        self.error_threshold = error_threshold
    
    # Write on multiple TCAM arrays
    # TODO: Create a function that writes on a specific row. (not necessary for genome sequencing)
    def write_masked(self,genc,genm):
        m = len(genc) % self.ncolumn 
        genc = genc + [0]*(self.ncolumn -m)
        genm = genm + [9223372036854775807]*(self.ncolumn -m)
        self.array_code = [genc[i:i+self.ncolumn] for i in range(0, len(genc), self.ncolumn)]
        self.array_mask = [genm[i:i+self.ncolumn] for i in range(0, len(genm), self.ncolumn)]
        self.narray = len(self.array_code)
        tag = [[0] * self.nrow for i in range(self.narray)]
    
    
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
                if ac != '*'*self.ncolumn:
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

    def search_bits(self,code,mask):
        na = 0
        pos = []
        pos2 = []
        #print(self.narray,self.ncolumn)
        for i in range(0,self.narray):
            nc = 0
            for j in range(0,self.ncolumn):
                if self.array_mask[i][j] != 9223372036854775807 and self.array_mask[i][j] != 18446744073709551615:
                    e = self.hamming_distance_mask(self.array_code[i][j],self.array_mask[i][j],code,mask)
                    if e == 0:
                        #print(i,j,self.array_code[i][j],self.array_mask[i][j],code,mask)
                        pos.append([na, nc, e])
                        pos2.append([na, nc, e])
                        self.tag[na][nc] = 1
                    elif e <= self.error_threshold:
                        pos2.append([na, nc, e])
                        #self.tag[na][nc] = 1
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
    

    def hamming_distance_mask(self,ac,am,bc,bm):
        hamm = bin(ac ^ bc & ~am & ~bm).count('1')
        #print(ac,bc,am,bm,hamm)
        return hamm


# In[ ]:





# In[ ]:





# In[10]:


# Save positional information
def save_pos(pos, pos_fname, prefix='',genome_start=0):
    with open(pos_fname, "a") as f:
        for s in pos:
            f.write(prefix + str(s[0]) + ' ' + str(s[1]) + ' ' + str(s[2]) + ' ' + str(genome_start) + "\n")

def findreads(tcam, reads, indexi, pos_fname, start,genome_start):
    
    pos = ''
    pos2 = ''
    ismatch = [False] * len(reads)    
    best_matches = []

    #for j in  readsi: 
    f = open(pos_fname + '.best','w')

    for j in  tqdm(range(len(indexi))):
    #for j in  range(len(indexi)):
        #print(j, " of ", len(indexi))
        tcam.reset_tag()
        idx = indexi[j]
        query = reads[idx]
        best_matches = []
        best_match = [0] * 2
            
        for reverse in [False, True]:
            matchflag = False
            N = len(query) - 16
            for i in tqdm(range(0,N,1)):
                q = query[i:i+16]
                #print(i, " of ", N, " shifts", q)
                qcode, qmask = genome2mask(q)
                if reverse == True:
                    qcode = int('{:064b}'.format(qcode[0])[::-1], 2)
                    qmask = int('{:064b}'.format(qmask[0])[::-1], 2)
                else:
                    qcode = qcode[0]
                    qmask = qmask[0]
                pos, pos2 = tcam.search_bits(qcode,qmask)
                #print(len(pos),len(pos2))
                #print(pos2)
                #save_pos(pos,pos_fname, 'Match ' + str(indexi[j]) + ' ' + str(i) + ' ')
                save_pos(pos2,pos_fname + '.par', 'Partial Match ' + str(indexi[j]) + ' ' + str(i) + ' ', genome_start)
                #if len(pos) > 0:
                #    best_matches.append([indexi[j], i, -1, pos]) 
            ismatch[j] = matchflag

            res = [rows for rows, val in enumerate(tcam.get_tag()[tcam.best_array()]) if val > 0] 
            #best_matches.append([indexi[j], -1, int(tcam.best_array()), [res]])
            
            if reverse == False:
                best_match[0] = [int(tcam.best_array()), res]
            else: 
                best_match[1] = [int(tcam.best_array()), res]
                

        #print(best_match)
        #print(tcam.narray,tcam.nrow,tcam.ncolumn)
        
        print(best_match)
        position = -1
        if best_match[0][1]:
                if len(best_match[0][1]) >= len(best_match[1][1]):
                    position = (int(best_match[0][0])) *nrow*ncol/4 + (int(best_match[0][1][0]))*nrow/4 + genome_start*nrow*ncol/4
                else:
                    position = (int(best_match[1][0])) *nrow*ncol/4 + (int(best_match[1][1][0]))*nrow/4 + genome_start*nrow*ncol/4
            
        best_matches.append([idx + start, position])
    
        for s in best_matches:
            f.write(str(s[0]) + ',' + str(s[1]) + ',' +  str(genome_start) + ' ' + str(res)[1:-1] )
            f.write('\n')
        #clear_output()
            
    f.close()
    return pos, pos2, best_matches

start = int(sys.argv[5])
fname = sys.argv[4] #'bp100/bp100'
readsg = readfile(fname + ".bfast.fastq",start,int(sys.argv[6]))
reads = get_reads(readsg)
tcam = TCAM_Arrays(nrow,ncol)
tcam.write_masked(genc,genm)
indexi = range(0,len(reads))


pos_100, pos2_100, bestmatch= findreads(tcam, reads, indexi, fname + sys.argv[5] + '_' + str(genome_start) + '.pos',start,genome_start)
#pos_100, pos2_100, bestmatch = findreads(tcam, reads, indexi, 'results/' + fname + '.posr', reverse = True)
#print(tcam.array_code[0])



