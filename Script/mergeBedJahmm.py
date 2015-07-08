#!/usr/bin/env python
# -*- coding:utf-8 -*-

#date created: 031012
#date modified : 081012
 # change the lines 46-48 for checking if all the files have equal rows
 # these lines are removed by checking if all the lines of the files are read 

############################################
# merge bed files into one big file        #
#                                          #
# !Note: before merging make sure!         #
# - all files have the same rows           #
# - all the rows are sorted                #
#                                          #
# The big file will contain:               #
# - each colname: the id of the filename   #
# - col: number of reads                   #
############################################

import sys
import re
import os

def EOF_or_raise(f):
   """
   Check if all the rows of the files are read 
   """
   try:
      f.next()
   except StopIteration:
      return
   else:
      raise Exception(str(f))

#put all the files that you wish to merge into one list
f_names = sys.argv[1:]
f_list = []
for fname in f_names:
   f_list.append(open(fname))

#grep only the IDs of the file
#e.g hg19-H1-003.bed, the id would be 003
#IDs = [re.search('-(\d{3}[a-z]*).', fname).group(1) for fname in f_names]

# fname = 300bin-Input.bed => Input
colNames = [os.path.splitext(fname)[0].split('-')[1] for fname in f_names]

#write the colnames
sys.stdout.write('\t'.join(colNames) + '\n')

try:
   #get the list containing number of reads from bed file
   while True:
      try:
         #number of reads at index 3, or column 4 (chr, start, end, read)
         #read for each file line by line (for each loop)
         items = [f.readline().split()[3] for f in f_list]
      #when the list of files is empty then stop
      except IndexError:
         break
      sys.stdout.write('\t'.join(items) + '\n')

   for f in f_list:
      EOF_or_raise(f)

finally:
   for f in f_list: f.close()
