#USAGE: python scriptname.py input_filename [number_lines_to_display]
import sys #IMPORT module
filename = sys.argv[1] #SET the input file 
if len(sys.argv) > 2: #IF the user specified a desired number of lines to display 
    n_lines = int(sys.argv[2]) #SET the "desired number" of displayed lines 
#END IF 
else: #OTHERWISE 
    n_lines = 10 #SET the "desired number" of displayed lines to a default 
#END OTHERWISE 
list_lines = [] #INITIALIZE a storage list for lines in the file
for line in open(filename): #FOR every line in the open file 
   list_lines.append(line) #ADD the line to the storage list
#END FOR
tail = list_lines[(len(list_lines) - n_lines):]#SET a subset of the storage list to be the last "desired number" of items in the storage list 
for line in tail: #FOR every line in the subset
   print(line.strip("\r\n")) #PRINT the line 
#END FOR 
