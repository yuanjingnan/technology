'''
cat XXX_list |while read s;do echo "python cal.v3.py "$s"/XXX.tab|awk -v OFS=\"\\t\" '{if(\$5>30){print \$4,\$5}}'|cat|awk '{if(\$1~/bb/){gsub(/bb/,\"aa\")}print \$0}'>$s/partial_candidate">$s/$s.partial.sh;done
'''
