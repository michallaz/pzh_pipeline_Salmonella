#!/bin/bash
# ver 0.1

GENOM=$1

run_blastn () {
	if [ ${1} == 'STY1774' ]; then
		blastn -query ${2} -db /cgMLST2_entero/${1}.fasta -out blastn_tmp_${1}.tab -outfmt 6 -max_target_seqs 3 -word_size 5 
	else
		 blastn -query ${2} -db /cgMLST2_entero/${1}.fasta -out blastn_tmp_${1}.tab -outfmt 6 -max_target_seqs 3
	fi
	return 0
}

export -f run_blastn
cat /cgMLST2_entero/all_3002_allele.txt | xargs -I {} --max-procs=50 bash -c "run_blastn {} ${GENOM}"

for K in `cat /cgMLST2_entero/all_3002_allele.txt`
do
	mv blastn_tmp_${K}.tab blastn_tmp.tab
	# 1 dostajemy jeden 100 % wynik pod katem seq id, sprawdzamy tylko czy pokrywa sie on w 100% z dlugoscia allelu (nie jest trimowany na koncach)
	#echo ${K}
	HITY=`cat blastn_tmp.tab | awk '{if($3 == 100) print $0}' | wc -l`
	if [ ${HITY} -eq 1 ]; then
		ALLEL_ID=`cat blastn_tmp.tab | awk '{if($3 == 100) print $2}' |  sed -r  s'/.*_(.+)/\1/'`
		ALLEL_LEN=`cat blastn_tmp.tab | awk '{if($3 == 100) print $4}'`
		LINIA=`grep -nw ${K}_${ALLEL_ID} /cgMLST2_entero/${K}.fasta | cut -d ":" -f1`
		LINIA=`echo "$LINIA + 1" | bc -l`
		ALLEL_TRUE_LEN=`awk -v l=$LINIA 'NR == l { print length} '  /cgMLST2_entero/${K}.fasta`
		if [ ${ALLEL_TRUE_LEN} -eq ${ALLEL_LEN} ]; then
			echo -e "${K}\t${ALLEL_ID}\t100\tnormal" >> log.log
		elif [ ${ALLEL_TRUE_LEN} -gt ${ALLEL_LEN} ];then
			echo -e "${K}\t${ALLEL_ID}\t100\tshort" >> log.log
		fi
	elif [ ${HITY} -gt 1 ]; then
		ALLEL_ID=(`cat blastn_tmp.tab  | awk '{if($3 == 100) print $2}' | sed -r  s'/.*_(.+)/\1/' |sort -nu`)
		
		FINAL_ID=`echo -1`
		FINAL_TRUE_LEN=`echo -1` # wybieramy na koniec taki allel sposorod idealnie mapowanych ktory jest najdluzszy ! 
		for ID in ${ALLEL_ID[@]}; do
			#petla od najnizszego numeru allelu, w kolejnej iteracji bedzie trzeba zwracac wiecej cyfr jesli jest kontaminacja
			#ALLEL_LEN=`cat blastn_tmp.tab | awk '{if($3 == 100) print $0}' | grep -P '_`echo ${ID}`\s' | cut -f4`
			# ponizej wyciagamy najdluzsze mapowanie , czasami zdarzalo sie ze wiele contigow mapuje 
			# sie na 100% z tym allelem, ale sa to mapowania krotsze, wiec bierzemy najdluzsze
			ALLEL_LEN=`cat blastn_tmp.tab | awk -v id=$ID '{ola="_"id"t"; $2=$2"t"} {if($2 ~ ola && $3 == 100) print $4}' | sort -rn | head -1`
			#echo $ALLEL_LEN
			LINIA=`grep -nw ${K}_${ID} /cgMLST2_entero/${K}.fasta | cut -d ":" -f1`
			LINIA=`echo "$LINIA + 1" | bc -l`
			#ALLEL_TRUE_LEN=`head -n ${LINIA}  /enterobase/blastdb/${K}.fasta | tail -1 | wc -c`
			ALLEL_TRUE_LEN=`awk -v l=$LINIA 'NR == l { print length} '  /cgMLST2_entero/${K}.fasta`
			if [ ${ALLEL_TRUE_LEN} -eq ${ALLEL_LEN} ]; then
				if [ ${ALLEL_TRUE_LEN} -gt ${FINAL_TRUE_LEN} ]; then
					FINAL_ID=`echo ${ID}`
					FINAL_TRUE_LEN=`echo ${ALLEL_TRUE_LEN}`
				fi

			fi
		
		done
		echo -e "${K}\t${FINAL_ID}\t100\tmulti" >> log.log
	elif  [ ${HITY} == 0 ]; then
		#poki co tylko allel o najwyzszym scorez bez analizy pokrycia
		ALLEL_ID=`cat blastn_tmp.tab | sort -rnk3 | head -1 | cut -f2 |  sed -r  s'/.*_(.+)/\1/' `
		echo -e "${K}\t-1\t<100\tunk" >> log.log
	fi

	rm blastn_tmp.tab
done
