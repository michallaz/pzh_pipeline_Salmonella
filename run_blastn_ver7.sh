#!/bin/bash
# ver 7
# skrypt zwraca plik log.log
# w ktrym allele podzielone sa na 4 grupy
# normal - jest tylko jeden allel ze 100% zgodnoscia w genomie, i ten allel jest obecny w genomie w 'pelnej' dlugosci
# short - jest tylko jeden allel ze 100% zgodnoscia w genomie, ale ten allel nie jest caly w genomie /jest uciety/
# multi - jest wiele alleli ze 100 seq ident, wybieramy allel ktory jest pelny i najdluzszy
# unk - nie ma allelu ze 100 seq ident, wybieramy allel ktory ma najwieksza wartosc score = dlugosc allelu w genommie - ilosc mismatch - ilosc gap

# Dla kazdej z tych klas log.log zawiera troche inne dodatkowe kolumny 5+
# normal - dlugosc allelu w cgmlst i w genomie (de facto musza to byc te same liczby)
# multi - lista wszystkich alleli ze 100% seq id
# short dlugosc allelu w cgmlst i w genomie (de facto musza to byc te same liczby, wtedy widac o ile za krtki jest allel w genomie)
# unk -  dlugosc allelu w cgmlst i score 




GENOM=$1
MAX_PROC=$2
run_blastn () {
	if [ ${1} == 'STY1774' ]; then
		# przerazlwie wolny ale ten allel jest ekstremalnie krotki
		/blast/bin/blastn -query ${2} -db /cgMLST2_entero/${1}.fasta -out blastn_tmp_${1}.tab -outfmt "6 saccver slen qstart qend sstart send pident evalue mismatch gaps" -max_target_seqs 5 -word_size 5 
	else
		/blast/bin/blastn -query ${2} -db /cgMLST2_entero/${1}.fasta -out blastn_tmp_${1}.tab -outfmt "6 saccver slen qstart qend sstart send pident evalue mismatch gaps" -max_target_seqs 5
	fi
	return 0
}

export -f run_blastn
cat /cgMLST2_entero/all_3002_allele.txt | xargs -I {} --max-procs=${MAX_PROC} bash -c "run_blastn {} ${GENOM}"

for K in `cat /cgMLST2_entero/all_3002_allele.txt`
do

	mv blastn_tmp_${K}.tab blastn_tmp.tab
	# 1 dostajemy jeden 100 % wynik pod katem seq id, sprawdzamy tylko czy pokrywa sie on w 100% z dlugoscia allelu (nie jest trimowany na koncach)
	#echo ${K}
	HITY=`cat blastn_tmp.tab | awk '{if($7 == 100) print $0}' | wc -l`
	
	# gdyby cos nie poszlo z tymi if-ami nizej 
	# to sa defaultowe wartosci dal zmiennych 
	ALLEL_ID="-1"
        ALLEL_TRUE_LEN="-1"
	
	if [ ${HITY} -eq 1 ]; then
		# Mamy jeden hit mapujacy sie w 100% na allel z cgMLST
		# ale sprawdzamy jeszcze czy mapowanie obejmuje cala dlugosc allelu z cgMLST, czy blast ucial alignment gdzies po drodze
		# id allelu z cgMLST
		ALLEL_ID=`cat blastn_tmp.tab | awk '{if($7 == 100) print $1}' |  sed -r  s'/.*_(.+)/\1/'`
		# dlugosc tego allelu
		ALLEL_LEN=`cat blastn_tmp.tab | awk '{if($7 == 100) print $2}'`
		
		# dlugosc mapowania query
		QUERY_LEN=`cat blastn_tmp.tab | awk '{if($7 == 100) print $4 - $3 + 1}'`


		# Zastanowic sie co jesli mam za krotka sekwencje w stosunku do allelu z cgMLST
		# czy nalezy braz odczyty z nie 100% seq id ale ktore sa krotsze i "lepiej" mapuja sie na alle z cgmlst 
		if [ ${QUERY_LEN} -eq ${ALLEL_LEN} ]; then
			echo -e "${K}\t${ALLEL_ID}\t100\tnormal\t${ALLEL_LEN}\t${QUERY_LEN}" >> log.log
		
		elif [ ${ALLEL_LEN} -gt ${QUERY_LEN} ];then
			echo -e "${K}\t${ALLEL_ID}\t100\tshort\t${ALLEL_LEN}\t${QUERY_LEN}" >> log.log
		fi


	elif [ ${HITY} -gt 1 ]; then
		# ALLEL_ID=(`cat blastn_tmp.tab  | awk '{if($7 == 100) print $1}' | sed -r  s'/.*_(.+)/\1/' |sort -nu`)
		ALLELS=(`cat blastn_tmp.tab  | awk '{if($7 == 100) print $1}'`)
		
		# Te zienne podmieniam jesli natrafiam na dluzszy allel
		FINAL_ID=`echo -1`
		FINAL_TRUE_LEN=`echo -1`
		FINAL_QUERY_LEN=`echo -1`

		for ID in ${ALLELS[@]}; do
			# wybieramy te allele ktore maja 100% pident, a mapowania obejmuja pelna dlugosc sekwencji tego allelu z cgMLST
			# uwaga ID moze byc w wynikach blasta n-razy wybieramy to ktore mapuje sie na pident 100
			
			#uwaga 2 moze byz tak ze ten sam allel ma 2 razy pident 100, ale rozne dlugosci , wybeiramy ten z najdluzszym query len

			ALLEL_ID=`cat blastn_tmp.tab |  awk -v id="${ID}" '{if($7 == 100 && $1 == id)  print $1}' |  sed -r  s'/.*_(.+)/\1/' | head -1 `
			ALLEL_LEN=`cat blastn_tmp.tab | awk -v id="${ID}" '{if($7 == 100 && $1 == id)  print $2}' | head -1`
			QUERY_LEN=`cat blastn_tmp.tab | awk -v id="${ID}" '{if($7 == 100 && $1 == id) print $4 - $3  + 1}' | sort -rnk1 | head -1`

			# zidentyfikowany allel jest dluzszy niz allel aktualnie traktowany jako "final"
			# oraz alignment obejmuje cala dlugosc sekwencji tego allelu
			if [ ${QUERY_LEN} -eq  ${ALLEL_LEN} ] || [ ${QUERY_LEN} -gt ${FINAL_TRUE_LEN} ]; then
				FINAL_ID=`echo ${ALLEL_ID}`
				FINAL_TRUE_LEN=`echo ${ALLEL_LEN}`
				FINAL_QUERY_LEN=`echo ${QUERY_LEN}`

			fi
		done
		echo -e "${K}\t${FINAL_ID}\t100\tmulti\t${FINAL_TRUE_LEN}\t${FINAL_QUERY_LEN}" >> log.log

	elif  [ ${HITY} == 0 ]; then
		# Nie ma allelu o 100% pident z allelami w cgMLST
		# Iterujemy po calym pliku i szukamy takiego allelu o najnizszej sumie dlugosc mapowania - (mismatch + gaps)
		# Defaultowe wartosci gdyby analizowany plik byl pusty
		FINAL_ID=`echo -1`
                FINAL_TRUE_LEN=`echo -1`
                FINAL_QUERY_LEN=`echo -1`
		FINAL_PIDENT="-1"	
		FINAL_SCORE="-1"
		while read LINE; do
		
			ALLEL_ID=`echo ${LINE} | cut -d " " -f1 |  sed -r  s'/.*_(.+)/\1/' `
			ALLEL_LEN=`echo ${LINE} | cut -d " " -f2`
                	QUERY_LEN=`echo ${LINE} | awk '{print $4 - $3 + 1}'`
			NO_MISMATCH=`echo ${LINE} | awk '{print $9}'`
			NO_GAPS=`echo ${LINE} | awk '{print $10}'`
			PIDENT=`echo ${LINE} | awk '{print $7}'`
			SCORE=`echo "${QUERY_LEN} - ${NO_MISMATCH} - ${NO_GAPS}" | bc -l | awk '{print int($1)}'`
			if [ ${SCORE} -gt ${FINAL_SCORE} ]; then
				FINAL_ID=`echo ${ALLEL_ID}`
				FINAL_TRUE_LEN=`echo ${ALLEL_LEN}`
				FINAL_QUERY_LEN=`echo ${QUERY_LEN}`
				FINAL_SCORE=`echo ${SCORE}`
				FINAL_PIDENT=`echo ${PIDENT}`

			fi


		done < blastn_tmp.tab 
		
		echo -e "${K}\t${FINAL_ID}\t${FINAL_PIDENT}\tunk\t${FINAL_TRUE_LEN}\t${FINAL_QUERY_LEN}\t${FINAL_SCORE}" >> log.log
	fi

	rm blastn_tmp.tab
done
