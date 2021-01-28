# Creates a dictionary file for the specified database terms 
# from a set of Interproscan output files.
# Change the grep search term e.g. "\sPfam\s" to the desired database.
# Change the output name after the redirect e.g. > dicionario.tsv
grep "\sPfam\s" *tsv | awk -F"\t" '{OFS="\t"}{print $5, $6}' | sort | uniq >  dicionario.tsv