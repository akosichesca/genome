
# Get the reverse complement of each genome bits. Since
def reverse_complement(a):
    rstr =''.join(reversed(a))
    rstr = rstr.replace('A', 'a') 
    rstr = rstr.replace('C', 'c') 
    rstr = rstr.replace('G', 'C') 
    rstr = rstr.replace('c', 'G') 
    rstr = rstr.replace('T', 'A') 
    rstr = rstr.replace('a', 'T') 
    return rstr

def readfile(fname):
    readsg = ''
    with open(fname) as fp:
        line = fp.readline()
        while line:
            line = fp.readline()
            readsg = readsg + line
            #reads = reads + line.strip()
            line = fp.readline()
            line = fp.readline()
            line = fp.readline()
    return readsg
   #while line:
   #    print("Line {}: {}".format(cnt, line.strip()))
   #    line = fp.readline()
   #    cnt += 1

def get_reads(readsg, nrows=64):
    reads = ''
    for r in readsg.split('\n'):
        gr = genome2bits(r)
        reads = reads + gr + '\n'

    readsr = ''
    for r in readsg.split('\n'):
        readsr = readsr + genome2bits(reverse_complement(r)) + '\n'

    reads = reads.split('\n')[:-2]
    readsr = readsr.split('\n')[:-2]
    
    return reads, readsr

fname = 'bp8.bfast.fastq'
readsg = readfile(fname)
reads, readsr = get_reads(readsg)

print(reads[0])
print(len(reads))
print(len(readsr))
