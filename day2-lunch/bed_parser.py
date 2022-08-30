#!/usr/bin/env python3

import sys

def parse_bed(fname):
    try:
        fs = open(fname, 'r')
    except:
        raise FileNotFoundError("That file doesn’t appear to exist")
    bed = []

    # start a counter for malformed lines
    malformed = 0

    # added field types for columns 7 and 8
    # also changed field type 5 to int as per the BED documentation
    field_types = [str, int, int, str, int, str, int, int]
    for i, line in enumerate(fs):
        if line.startswith("#"):
            continue
        fields = line.rstrip().split()
        fieldN = len(fields)

        # Include only allowed lengths
        if not (3 <= fieldN <= 9 or fieldN == 12):
            
            # Increment malformed line counter
            malformed += 1
            continue

        try:

            # Fixing types for any allowed bed format (3-9 and 12)
            for j in range(fieldN):
                if j < 8:
                    fields[j] = field_types[j](fields[j])

                # Splitting the RGB entries and converting to ints
                elif j == 8:
                    fields[j] = fields[j].split(',')
                    fields[j] = [int(k) for k in fields[j]]

                    # check that the RGB entry is length 3
                    assert len(fields[j]) == 3

                # Fix type for the blockCount
                elif j == 9:
                    fields[j] = int(fields[j])

                # Split the blockSize and blockStart int lists of ints
                else:

                    # make sure to strip trailing commas
                    fields[j] = fields[j].rstrip(',').split(',')
                    fields[j] = [int(k) for k in fields[j]]

                    # Check that these lists are the length of the block count
                    assert len(fields[j]) == fields[9]
    
            bed.append(fields)
        except:
            
            # increment malformed if any of the above failed
            malformed += 1

    # print count of malformed lines if nonzero
    if malformed > 0:
        print(f"{malformed} lines are malformed", file=sys.stderr)
    fs.close()
    return bed

if __name__ == "__main__":
    fname = sys.argv[1]
    bed = parse_bed(fname)
