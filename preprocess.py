import math
import numpy as np
nrow = 64
ncol = 64
nrowcol = nrow*ncol

fname = "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
fname_out = "human_64.fna"
f = open(fname_out,'w')


with open(fname) as fp:
    line = fp.readline().strip()
    while line:
        buf = []
        for i in range(int(nrowcol/4)):
            if ">" in line:
                print(line)
                i = i-1
            else:
                buf.append(line)
            line = fp.readline().strip()
            if not line:
                break
        buf = ''.join(buf)
        buf = [buf[i:i+int(nrowcol/4)] for i in range(0, len(buf), int(nrowcol/4))]
        for item in buf:
            f.write("".join(item.split()) + '\n')
f.close()

