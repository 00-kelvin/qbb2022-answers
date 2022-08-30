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

# This is a beautiful script. I really like how your comments are the pseudo-code.
# There are only two things that I will add to help make you code even better.
# First, adding blank lines to separate your code into different functional blocks
# can help make it more readble. Second, putting longer comments on lines of 
# their own will also help with readability. As far as the logic goes, you're
# spot on. Great job! - Mike

# Notes from small group discussion:
# Can use .reverse to reverse the order of a list
# Can use list[-n_lines:] to start from n_lines from the end
# Could have done some "defensive programming" -- check if the inputs are valid


