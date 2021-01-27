grep "\sPfam\s" *tsv | awk -F"\t" '{OFS="\t"}{print $5, $6}' | sort | uniq >  dicionario.tsv
