'''
#!/usr/bin/env python3

import sys
import re
import os

fafile = sys.argv[1]
idfile = sys.argv[2]
os.path.abspath(fafile)
dir_name=os.path.dirname(fafile)
# store all ids
idfh = open(idfile,'r')
idSet = set()
for line in idfh:
    idSet.add(line.rstrip())
idfh.close()

base_filename=''
curId=''
filename_suffix='fa'
path=os.path.join(dir_name, base_filename + "." + filename_suffix)
f=open(path,'w')
pattern = re.compile(">([\w.]+)")
fafh = open(fafile,'r')
for line in fafh:
    line = line.strip()
    if line[0] == '>':
        curId = pattern.match(line).group(1)
        if curId in idSet:
            base_filename=curId
            path=os.path.join(dir_name, base_filename + "." + filename_suffix)
            f.close()
            if os.path.exists(path):
                f=open(path,'a')
                f.write(line)
                f.write("\n")
            else:
                f=open(path,'w')
                print line
                f.write(line)
                f.write("\n")
    else:
        f.write(line)
        f.write("\n")
fafh.close()
'''
