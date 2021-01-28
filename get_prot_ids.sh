for i in `ls *.fa`; do grep ">" $i | awk -F"protein_id:" '{print $2}' | awk -F"|" '{print $1}' > $i.ids; done
