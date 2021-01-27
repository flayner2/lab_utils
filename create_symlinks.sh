for i in `ls *.txt`; do ln -s $i ${i:0:13}; done 
