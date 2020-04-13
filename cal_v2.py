'''
#!/usr/bin/python

import sys
import re


def read(l):
    i=open(l,'rb')
    arr=[]
    col={}
    for line in i:
        sp=line.split( )
        gene_name=re.match(r'([\da-zA-Z]*)_*',sp[0]).group(1)
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
            col[name]=(sum_sub_len,float(sum_exact_map)/float(sum_sub_len),sum_block_sp_identify,align_rate,scaf_start,scaf_end)
            j=j+1
    arr1=sorted(arr,key=lambda x: (x[0], x[1]),reverse=False)
    return arr1,col



##########arr1:scaf_id block_start block_end block_identify gene_name scaf_start scaf_end sub_len exact_map align_rate
##########col:key:scaf_id gene_name scaf_start scaf_end
##########col:value:all_block_len map_rate ave_rate align_rate scaf_start scaf_end


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
        if int(col[name][1])>70 and int(col[name][2])>70 and float(col[name][3])>70:
            filter_arr.append([scaf_id,scaf_start,scaf_end,gene_name])
    arr1=sorted(filter_arr,key=lambda x: (x[0],x[1],x[2],x[3]),reverse=False)

####arr1:scaf_id scaf_start scaf_end gene_name

    scaf_id=""
    scaf_end=0
    scaf_start=0
    gene_name=""
    overlap_col=[]
    a=0
    overlap_len=0
    overlap_start=0
    overlap_end=0
    same=[]
    for l in arr1:
#        print a,l,"######################################"
#        print scaf_id
#        print scaf_start
#        print scaf_end
#        print gene_name
#        print overlap_start
#        print overlap_end
        if same==l:
            continue
        else:
            same=l

        if l[0]==scaf_id:
            overlap_col.pop()
            if l[1]>=scaf_end:
#                print "new scaf_start is more than old sacf_end"
                overlap_col.append([scaf_id,scaf_start,scaf_end,gene_name,int(scaf_start)+int(overlap_len),scaf_end])
#                print "result",[scaf_id,scaf_start,scaf_end,gene_name,int(scaf_start)+int(overlap_len),scaf_end]
                overlap_start=1
                overlap_end=0
                overlap_len=overlap_end-overlap_start+1
            else:
#                print "new sacf has overlap with old scaf"
                overlap_start=l[1]
                past_overlap_len=overlap_len
                if l[2]>=scaf_end:
#                    print "new scaf_end is more than old scaf_end"
                    overlap_end=scaf_end
                    t2="\t"
                    n2=(scaf_id,gene_name,scaf_start,scaf_end)
                    n3=(l[0],l[3],l[1],l[2])
                    name2=t2.join(n2)
                    name3=t2.join(n3)
                    overlap_len=int(overlap_end)-int(overlap_start)+1
#                    print "overlap_len=",overlap_end
                    overlap_rate=float(overlap_len/(int(scaf_end)-int(scaf_start)+1))
#                    print overlap_rate
                    if overlap_rate>0.4:
#                        print "more than 0.4"
                        if col[name2][0]> col[name3][0]:
#                            print "former longer"
                            overlap_col.append([scaf_id,scaf_start,scaf_end,gene_name,int(scaf_start)+int(past_overlap_len),overlap_end])
#                            print "result",[scaf_id,scaf_start,scaf_end,gene_name,int(scaf_start)+int(past_overlap_len),overlap_end]

                        else:
                            overlap_col.append([scaf_id,scaf_start,scaf_end,gene_name,int(scaf_start)+int(past_overlap_len),int(overlap_start)-1])
                            e=int(scaf_start)+int(past_overlap_len)
                            f=int(overlap_start)
                            if e == f:
                                overlap_col.pop()
#                            print "result",[scaf_id,scaf_start,scaf_end,gene_name,int(scaf_start)+int(past_overlap_len),int(overlap_start)-1]
                            overlap_col.append([scaf_id,scaf_start,scaf_end,l[3],overlap_start,overlap_end])
#                            print "latter longer"
#                            print "result",[scaf_id,scaf_start,scaf_end,l[3],overlap_start,overlap_end]
                    else:
#                        print overlap_len,float(overlap_rate),"less than 0.4"
                        overlap_col.append([scaf_id,scaf_start,scaf_end,gene_name,int(scaf_start)+int(past_overlap_len),overlap_end])
#                        print "result",[scaf_id,scaf_start,scaf_end,gene_name,int(scaf_start)+int(past_overlap_len),overlap_end]
                        overlap_col.append([scaf_id,scaf_start,scaf_end,l[3],overlap_start,overlap_end])
#                        print "result",[scaf_id,scaf_start,scaf_end,l[3],overlap_start,overlap_end]
                else:
#                    print "new scaf_end is less more scaf_end"
                    overlap_col.append([scaf_id,scaf_start,scaf_end,gene_name,int(scaf_start)+int(overlap_len),scaf_end])
                    continue
        else:
            overlap_col.append([scaf_id,str(scaf_start),str(scaf_end),gene_name,int(scaf_start)+int(overlap_len),scaf_end])
            b=int(scaf_end)+1
            c=int(overlap_len)+int(scaf_start)

            if c == b:
                overlap_col.pop()
#            print overlap_len
#            print "result",[scaf_id,str(scaf_start),str(scaf_end),gene_name,int(scaf_start)+int(overlap_len),scaf_end]
#            print "newnoe"
            overlap_end=0
            overlap_start=1
            overlap_len=int(overlap_end)-int(overlap_start)+1
            overlap_col.append([l[0],l[1],l[2],l[3],l[1],l[2]])
#            print "result",[l[0],l[1],l[2],l[3],l[1],l[2]]

        scaf_id=l[0]
        gene_name=l[3]
        scaf_start=l[1]
        scaf_end=l[2]
#        a=a+1
    return overlap_col

############overlap_col:scaf_id,scaf_start,scaf_end,gene_name,map_start,map_end

##mian

(store,integet)=read(sys.argv[1])
#for i in store:
#    print i
#for key in integet:
#    print integet[key]
score=overlap(store,integet)
for i in score:
    a=int(i[5])+1
    b=int(i[4])
    if a != b:
        a=(i[0],str(i[1]),str(i[2]),str(i[3]),str(i[4]),str(i[5]))
        t="\t"
        long=t.join(a)
        print long
'''
