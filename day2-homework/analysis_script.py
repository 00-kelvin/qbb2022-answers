#import vcf parser
from vcfParser import parse_vcf

#parse the two vcf files
random = parse_vcf("random_snippet.vcf")
dbSNP = parse_vcf("dbSNP_snippet.vcf")

#initialize empty position-ID dictionary
SNP_ID_dict = {}

#loop through all the lines in dbSNP
for i in range(len(dbSNP)):

	#set the key in dict to positions and value to IDs
	SNP_ID_dict[dbSNP[i][1]] = dbSNP[i][2]

#initialize counter of non corresponding IDs
noID = 0

#loop through every line of random snippet
for i in range(len(random)):
	
	#try to set the ID in random to the ID in SNP dictionary
	#with corresponding position value
	try:
		random[i][2] = SNP_ID_dict[random[i][1]]
	
	#if there is no matching position, increment noID counter
	except:
		noID += 1

#print non-matching record count
print(f"There are {noID} records that do not have a corresponding ID in dbSNP")

#print first 100 lines (first 80 characters)
print("The first 100 records in random_snippet (80 character limit):")
for record in random[:100]: 
	print(str(record)[:80])

