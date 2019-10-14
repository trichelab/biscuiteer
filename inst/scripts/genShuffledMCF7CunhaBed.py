'''
genShuffledMCF7CunhaBed.py creates a BED file with shuffled coverages and
randomly selected methylation rates (chosen via randint(0,new_covg) and then
calculating the beta value.
'''

import sys
import numpy as np
import random
import collections

def usage():
    print('\ngenShuffledMCF7CunhaBed.py creates a BED file with shuffled')
    print('\tcoverages and randomly selected methylation rates (chosen via')
    print('\trandint(0,new_covg) and then calculating the beta value.')
    print('\n\npython genShuffledMCF7CunhaBed.py infile.bed outfile.bed')
    print('\tinfile.bed  --- original BED file')
    print('\toutfile.bed --- shuffled BED file')

if (len(sys.argv) != 3):
    usage()
    sys.exit()

inn = sys.argv[1]
owt = sys.argv[2]

# Open series data file
try:
    f = open(inn, 'r')
except IOError:
    print('File: {} does not exist'.format(inn))

out = open(owt, 'w')

col1 = []
col2 = []
col3 = []
col4 = []
col5 = []
for line in f:
    entry = line.rstrip().rsplit(sep='\t')

    # Column 1 is the chromosome
    # Same as input
    col1.append(entry[0])

    # Column 2 is the start of the CpG
    # Shifted with -1, 0, +1 from input
    # Shift is normal(mu = 0, sigma = 0.2), this keeps most of the shifts at 0,
    #     but provides a small amount of variation for example purposes
    col2.append(entry[1])

    # Column 3 is the end of the CpG
    col3.append(entry[2])

    # Column 4  is the methylation fraction ( M / (M+U) )
    # M is random integer between 0 and shuffled Cov value
    #col4.append(entry[3])

    # Column 5 is the coverage (i.e. read depth at that locus)
    # Will shuffle values around and randomly assign methylation levels
    col5.append(entry[4])

np.random.shuffle(col5)

#shifts = []
prev = 0
for i in range(0, len(col1)):
    # Column 1
    chr1 = col1[i]
    
    # Column 2
    start = int(col2[i])
    shift = round(np.random.normal(0, 0.2))
    #shifts.append(shift)
    new_start = start + shift
    if (new_start == prev):
        new_start += 1
    prev = new_start

    # Column 3
    new_end = new_start + (int(col3[i]) - int(col2[i]))

    # Column 4
    cov1 = int(col5[i])
    meth = np.random.randint(0, cov1)
    frac = ''
    if (cov1 != 0):
        temp = round(float(meth / cov1), 3)
        frac = '{:.3f}'.format(temp)
    elif (cov1 == 0):
        frac = '.'
    else:
        print('ERROR: COV IS WEIRD')

    lyst = [chr1, str(new_start), str(new_end), frac, str(cov1)]

    strings = '\t'.join(lyst)

    out.write(strings + '\n')

#print(collections.Counter(shifts))

out.close()
f.close()
