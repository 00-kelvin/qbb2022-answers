#!/usr/bin/env python3

#import sys module to interact with command line
import sys

#start the function to read the vcf
def parse_vcf(fname):
    
    #initialize variables that we will use later

    #vcf will be a list containing each line of the file
    vcf = []

    #info_description will be a dictionary with headers in the info section
    info_description = {}

    #info_type is a dictionary that will tell us what type each field in the 
    # info column will be
    info_type = {}

    #tells us formatting of the genotype information
    format_description = {}

    #to translate the types in the header to python types
    type_map = {
        "Float": float,
        "Integer": int,
        "String": str
        }

    #start counter of malformed lines
    malformed = 0

    #defensive coding to open the file or raise error if it doesn't exist
    try:
        fs = open(fname)
    except:
        raise FileNotFoundError(f"{fname} does not appear to exist", file=sys.stderr)

    #go through each line in the file (including values and indices)
    for h, line in enumerate(fs):
        
        #for the header lines (don't need to annotate)
        if line.startswith("#"):     
            try:                
                if line.startswith("##FORMAT"):
                    fields = line.split("=<")[1].rstrip(">\r\n") + ","
                    i = 0
                    start = 0
                    in_string = False
                     
                    while i < len(fields):
                        if fields[i] == "," and not in_string:
                            name, value = fields[start:i].split('=')
                            if name == "ID":
                                ID = value
                            elif name == "Description":
                                desc = value
                            start = i + 1
                        elif fields[i] == '"':
                            in_string = not in_string
                        i += 1
                    format_description[ID] = desc.strip('"')
                elif line.startswith("##INFO"):
                    fields = line.split("=<")[1].rstrip(">\r\n") + ","
                    i = 0
                    start = 0
                    in_string = False
                    while i < len(fields):
                        if fields[i] == "," and not in_string:
                            name, value = fields[start:i].split('=')
                            if name == "ID":
                                ID = value
                            elif name == "Description":
                                desc = value
                            elif name == "Type":
                                Type = value
                            start = i + 1
                        elif fields[i] == '"':
                            in_string = not in_string
                        i += 1
                    info_description[ID] = desc.strip('"')
                    info_type[ID] = Type
                elif line.startswith('#CHROM'):
                    fields = line.lstrip("#").rstrip().split("\t")
                    vcf.append(fields)
            except:
                raise RuntimeError("Malformed header")
        #now looking at the actual variant lines
        else:
            #making sure the variants are properly formatted
            try:
                
                #making a list of the entries in the columns by stripping 
                #trailing characters and then splitting by tabs
                fields = line.rstrip().split("\t")
                
                #making the entry in the position field into an int
                fields[1] = int(fields[1])

                #change the entry in qual to a float, unless it's a '.'
                if fields[5] != ".":
                    fields[5] = float(fields[5])

                #initialize a dictionary for the data in the info column
                info = {}

                #split info column into list of key/value pairs by ;
                #then looping through those entries
                for entry in fields[7].split(";"):

                    #creates a list of 2 items (e.g. AC, 91) for each entry
                    temp = entry.split("=")

                    #for fields that have no equal sign (just a flag)
                    if len(temp) == 1:

                        #set key as the value, value as None
                        info[temp[0]] = None
                    
                    #for the entries that have an equal sign
                    else:

                        #assigning variable names to the info before and after
                        #the equal sign
                        name, value = temp

                        #use the info_type dictonary to determine what type the
                        #item with that name has and store it to variable Type
                        Type = info_type[name]

                        #set the key for name to be its associated value, but
                        #converted to its appropriate type
                        info[name] = type_map[Type](value)
                
                #set the list for data for the info field to the dictionary we 
                #just made
                fields[7] = info

                #if there are more than 8 fields (the file has genotype info)...
                if len(fields) > 8:

                    #split the format column by colons into a list and store it
                    #back in the 9th column
                    fields[8] = fields[8].split(":")

                    #if there was a colon...
                    if len(fields[8]) > 1:

                        #for all the individuals thereafter...
                        for i in range(9, len(fields)):

                            #we also need to split their entries by colons
                            fields[i] = fields[i].split(':')
                    else:
                        #if there is only one piece of info in the format column,
                        #changes it from a list containing 1 item to a string
                        #consisting of the first item in the list
                        fields[8] = fields[8][0]
                
                #append all our formatted lists to the parsed vcf list
                vcf.append(fields)
            
            #if the variants were formatted wrong, increment malformed counter
            except:
                malformed += 1

    #replace the INFO header in the #CHROM line with info_description dict 
    vcf[0][7] = info_description

    #if the length of the first line is > 8 (there is genotype info)
    if len(vcf[0]) > 8:

        #set the FORMAT header in the #CHROM line to format_description dict
        vcf[0][8] = format_description

    #if some of the lines were malformed, print how many to the standard error
    if malformed > 0:
        print(f"There were {malformed} malformed entries", file=sys.stderr)
    
    #return formatted list!
    return vcf

#allows you to call this script as the main function in Terminal
if __name__ == "__main__":
    fname = sys.argv[1]
    vcf = parse_vcf(fname)

    #then prints the first 10 lines of the formatted list
    for i in range(10):
        print(vcf[i])
