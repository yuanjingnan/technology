'''
#!/usr/bin/python

import sys
import re


def read(l):
    i=open(l,'rb')
    arr=[]
    col={}
    col2={}
    for line in i:
        sp=line.split( )
        gene_name=sp[0]
        scaf_id=sp[5]
        scaf_start=sp[7]
        scaf_end=sp[8]
        align_rate=sp[10]
        scaf_block_start_end=sp[12]
        scaf_block_identify=sp[13]
        block_sp=scaf_block_start_end.split(";")
        block_sp_identify=scaf_block_identify.split(";")
        j=0
        sum_sub_len=0
        sum_block_sp_identify=0
        sum_exact_map=0
        while j<int(sp[9]):
            sub_sp=block_sp[j].split(",")
            sub_len=int(sub_sp[1])-int(sub_sp[0])+1
            arr.append([scaf_id,int(sub_sp[0]),int(sub_sp[1]),block_sp_identify[j],gene_name,scaf_start,scaf_end,sub_len,sub_len*float(block_sp_identify[j]),align_rate])
            s="\t"
            n=(scaf_id,gene_name,scaf_start,scaf_end)
            name=s.join(n)
            sum_sub_len=sum_sub_len+int(sub_len)
            sum_block_sp_identify=sum_block_sp_identify+float(block_sp_identify[j])/int(sp[9])
            sum_exact_map=sum_exact_map+float(sub_len)*float(block_sp_identify[j])
            col[name]=(sum_sub_len,float(sum_exact_map/100),sum_block_sp_identify,align_rate,scaf_start,scaf_end,)
            col2[name]=(sum_sub_len,float(sum_sub_len)/(int(scaf_end)-int(scaf_start)+1)*100,float(sum_exact_map)/float(sum_sub_len),sum_block_sp_identify,align_rate)
            j=j+1
    arr1=sorted(arr,key=lambda x: (x[0], x[1]),reverse=False)
    return arr1,col,col2



##########arr1:scaf_id block_start block_end block_identify gene_name scaf_start scaf_end sub_len exact_map align_rate
##########col:key:scaf_id gene_name scaf_start scaf_end
##########col:value:all_block_len map_exact_len ave_rate align_rate scaf_start scaf_end


def overlap(arr,col):
    filter_arr=[]
    for i in arr:
        scaf_id=i[0]
        gene_name=i[4]
        scaf_start=i[5]
        scaf_end=i[6]
        s="\t"
        n=(scaf_id,gene_name,scaf_start,scaf_end)
        name=s.join(n)
        if int(col[name][2])>70 :
            filter_arr.append([scaf_id,scaf_start,scaf_end,gene_name])
#            print scaf_id,"\t",gene_name
    arr1=sorted(filter_arr,key=lambda x: (x[0],x[1],x[2],x[3]),reverse=False)

####arr1:scaf_id scaf_start scaf_end gene_name

    scaf_id=""
    scaf_end=0
    scaf_start=0
    gene_name=""
    overlap_arr=[]
    a=0
    overlap_len=0
    overlap_start=0
    overlap_end=0
    same=[]
    t="\t"
    for l in arr1:
        if overlap_arr==[]:
#            print "start"
            start_n=(l[0],l[3],l[1],l[2])
            start_name=t.join(start_n)
            overlap_arr.append([l[0],l[1],l[2],l[3],col[start_name][0],col[start_name][1]])
#            print "result",[l[0],l[1],l[2],l[3],col[start_name][0]]
        if overlap_arr!=[]:
            if l==same:
                continue
            else:
                same=l
#            print "####",overlap_arr
            old2=overlap_arr[-1]
#            print old2,"old"
            if l[0]==old2[0]:
                if l[1]>=old2[2]:
                    new_n=(l[0],l[3],l[1],l[2])
                    new_name=t.join(new_n)
#                    print "new scaf_start is more than old sacf_end"
                    overlap_arr.append([l[0],l[1],l[2],l[3],col[new_name][0],col[new_name][1]])
#                    print "result",[l[0],l[1],l[2],l[3],col[new_name][1]]
                else:
#                    print "new sacf has overlap with old scaf"
                    overlap_arr.pop()
                    overlap_start=l[1]
                    old_n=(old2[0],old2[3],old2[1],old2[2])
                    old_name=t.join(old_n)
                    new_n=(l[0],l[3],l[1],l[2])
                    new_name=t.join(new_n)
#                    print new_name,old_name
                    if l[2]>=old2[2]:
#                        print "new scaf_end is more than old scaf_end"
                        overlap_end=old2[2]
                        overlap_len=float(overlap_end)-float(overlap_start)+1
                        scaf_len=float(scaf_end)-float(scaf_start)+1
#                        print "overlap_len=",overlap_len
#                        print "scaf_len=",scaf_len
                        overlap_rate=float(overlap_len/scaf_len)
#                        print overlap_rate
                        if float(overlap_rate)>0.4:
#                            print "more than 0.4"
                            if int(col[old_name][1])> int(col[new_name][1]):
#                                print "former longer"
                                overlap_arr.append([old2[0],old2[1],old2[2],old2[3],col[old_name][0],col[old_name][1]])
#                                print "result",[old2[0],old2[3],old2[1],old2[2],col[old_name][0],col[old_name][1]]
                            else:
                                overlap_arr.append([l[0],l[1],l[2],l[3],col[new_name][0],col[new_name][1]])
#                                print "result",[l[0],l[1],l[2],l[3],col[new_name][0],col[new_name][1]]
                        else:
#                            print overlap_len,float(overlap_rate),"less than 0.4"
                            overlap_arr.append([old2[0],old2[1],old2[2],old2[3],col[old_name][0],col[old_name][1]])
#                            print "result",[old2[0],old2[3],old2[1],old2[2],col[old_name][0],col[old_name][1]]
                            overlap_arr.append([l[0],l[1],l[2],l[3],col[new_name][0],col[new_name][1]])
#                            print "result",[l[0],l[1],l[2],l[3],col[new_name][0],col[new_name][1]]
                    else:
#                        print "new scaf_end is less old scaf_end"
                        if int(col[old_name][1])> int(col[old_name][1]):
                            overlap_arr.append([old2[0],old2[1],old2[2],old2[3],col[old_name][0],col[old_name][1]])
                        else:
                            overlap_arr.append([l[0],l[1],l[2],l[3],col[new_name][0],col[new_name][1]])
            else:
#                print "newnoe"
                n1=(l[0],l[3],l[1],l[2])
                name1=t.join(n1)
                overlap_arr.append([l[0],l[1],l[2],l[3],col[name1][0],col[name1][1]])
#                print "result",[l[0],l[1],l[2],l[3],col[name1][0],col[name1][1]]
    return overlap_arr

############overlap_arr:scaf_id,scaf_start,scaf_end,gene_name,scaf_sum_block_len,ave_iden

##mian

(store,integet,test)=read(sys.argv[1])
#for i in store:
#    print i
#for key in test:
#    p=(key,str(test[key][0]),str(test[key][1]),str(test[key][2]),str(test[key][3]),str(test[key][4]))
#    t="\t"
#    check=t.join(p)
#    print check
score=overlap(store,integet)
for i in score:
    a=(i[0],str(i[1]),str(i[2]),i[3],str(i[5]))
    t="\t"
    long=t.join(a)
    print long
'''
