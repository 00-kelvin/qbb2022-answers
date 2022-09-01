#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt

vcf = sys.argv[1]
fs = open( vcf )

ac = []
for i, line in enumerate( fs ):
    if "#" in line:
        continue
    fields = line.split()
    info = fields[7].split(";")
    ac.append( int(info[0].replace("AC=","")) )

fig, ax = plt.subplots()

# changed to log scale
ax.hist( ac, density=True, log=True)

# added consistent y-axis limits
ax.set_ylim(1e-6, 1e-2)

# add title
plt.title(vcf + " allele count frequency")

# add axis labels
ax.set_xlabel("Total number of alternate alleles in called genotypes")
ax.set_ylabel("Frequency")

fig.savefig( vcf + ".png" )


fs.close()

