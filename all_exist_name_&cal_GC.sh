'''
cat XXX.id |while read s;do python all_exist_name.py $s/XXX.cds XXX.re ;done
cat XXX.id |while read s ;do echo "XXX.re |while read n;do python cal_GC.py $s/\$n.fa>$s/\$n.fa.gc;done" >$s/$s.sh;done
'''
