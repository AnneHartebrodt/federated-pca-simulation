for line in $(cat $1); do
	result=$(echo "$line" | awk 'NR>1{print PREV} {PREV=$0} END{printf("%s",$0)}')
result=${result//"samples"/""}
result=${result//".tsv"/""}
result=${result/\\n/""}
echo -e "$result\t$line" >> $2
done;