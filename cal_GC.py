'''
#!/usr/bin/env python3
import sys
import os
import re

def getGCNum(inSeq):
    count = 0
    for i in inSeq.upper():
        if i == 'G' or i=='C':
            count+=1
    return count

def getLowerBase(inSeq):
    count = 0
    for i in inSeq:
        if not i.isupper():
            count+=1
    return count

def getNCount(inSeq):
    count = 0
    for i in inSeq:
        if i == 'N' or i=='n':
            count+=1
    return count

inFaFile = sys.argv[1]

totalLen = 0
lenArr = []
GCcount = 0
seqCounts = 0
nBaseCounts = 0
lowerBaseCounts = 0

curId = ""
curSeq = ""
firstRecordflag = True
faFh = open(inFaFile,'r')
pattern = re.compile(">([\w.]+)")
Id=''
for line in faFh:
    line = line.rstrip()
    if line[0] == '>':
        Id = pattern.match(line).group(1)
        if not firstRecordflag:
            totalLen += len(curSeq)
            lenArr.append(len(curSeq))
            GCcount += getGCNum(curSeq)
            nBaseCounts += getNCount(curSeq)
            lowerBaseCounts += getLowerBase(curSeq)
        #
        firstRecordflag = False
        seqCounts += 1
        curSeq=""
    else:
        curSeq = curSeq + line
# collect the last seq
# print(totalLen)
totalLen += len(curSeq)
lenArr.append(len(curSeq))
GCcount += getGCNum(curSeq)
nBaseCounts += getNCount(curSeq)
lowerBaseCounts += getLowerBase(curSeq)
# 

faFh.close()
'''
# getMaxLen MinLen MeanLen N50
lenArr.sort()
# 
maxLen = lenArr[-1]
minLen = lenArr[0]
meanLen = totalLen/seqCounts
if seqCounts%2:
    medianLen = (lenArr[seqCounts//2]+lenArr[seqCounts//2+1])/2
else:
    medianLen = lenArr[seqCounts//2]

lenArr.reverse()
halfTotalLen = totalLen/2
cumLen = 0
for curLen in lenArr:
    cumLen += curLen
    if cumLen >= halfTotalLen:
        N50 = curLen
        break
'''
print Id,"\t",seqCounts,"\t",float(GCcount)/float(totalLen)
#print "Total Len is:","\t",totalLen
#print "Total Seq Num is:","\t",seqCounts
#print "Total N count is:","\t",nBaseCounts
#print "Total lowBase is:","\t",lowerBaseCounts
#print "Total GC contene is:","\t",float(GCcount)/float(totalLen)
#print GCcount
#print totalLen
#print("maxLen is:""\t",maxLen)
#print("minLen is:""\t",minLen)
#print("meanLen is:""\t",meanLen)
#print("medianLen is:""\t",medianLen)
#print("N50 is:""\t",N50)
'''











