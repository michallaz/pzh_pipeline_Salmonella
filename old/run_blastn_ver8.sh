#!/bin/bash
# ver 8

# W wersji 8 dalej opracowywanej na podstawie przykladu ST8 zmieniono if-a przy multimapowani
# Nie liczy sie dlugosc a tylko numer allelu, im nizszy tym lepszy

# chyba entero tez sprawdza czy kolejny allel zaczyba sie ZA poprzednim

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
	PREVIOUS_POS="-1"
	# zamienic na mv po testach
	cp blastn_tmp_${K}.tab blastn_tmp.tab
	# 1 dostajemy jeden 100 % wynik pod katem seq id, sprawdzamy tylko czy pokrywa sie on w 100% z dlugoscia allelu (nie jest trimowany na koncach)
	#echo ${K}
	HITY=`cat blastn_tmp.tab | awk '{if($7 == 100) print $0}' | wc -l`
	
	# Sa hity ze 100 seq id 
	if [ ${HITY} -gt 0 ]; then
		ALLELS=(`cat blastn_tmp.tab  | awk '{if($7 == 100) print $1}'`)
		
		# Te zienne podmieniam jesli natrafiam na dobry allel
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
			HIT_START=`cat blastn_tmp.tab |  awk '{if($4 > $3) {print $3} else {print $4}}'`
			HIT_END=`cat blastn_tmp awk | '{if($4 > $3) {print $4} else {print $3}}'`
			
			# znalazlem pelny hit i final id -1
			if [ ${QUERY_LEN} -eq  ${ALLEL_LEN} ] && [ ${QUERY_LEN} -ge ${FINAL_TRUE_LEN} ] && [ ${FINAL_ID} -eq -1 ]; then
				FINAL_ID=`echo ${ALLEL_ID}`
				FINAL_TRUE_LEN=`echo ${ALLEL_LEN}`
				FINAL_QUERY_LEN=`echo ${QUERY_LEN}`
			#elif [ ${QUERY_LEN} -eq  ${ALLEL_LEN} ] && [ ${QUERY_LEN} -gt ${FINAL_TRUE_LEN} ]; then
				# nowy allel jest dluzszy niz ten znaleziony wczesniej, to ma pierwszenstwo
			#	FINAL_ID=`echo ${ALLEL_ID}`
                        #       FINAL_TRUE_LEN=`echo ${ALLEL_LEN}`
                        #        FINAL_QUERY_LEN=`echo ${QUERY_LEN}`
			elif [ ${QUERY_LEN} -eq  ${ALLEL_LEN} ] && [ ${ALLEL_ID} -lt ${FINAL_ID} ] && [ ${HIT_START} -gt ${FINAL_ID} -gt ${PREVIOUS_POS} ];then
				# allel jest dluzszy lub tej samej dlugosci, ale numer allelu jest nizszy niz ten wczesniejszy
				FINAL_ID=`echo ${ALLEL_ID}`
                                FINAL_TRUE_LEN=`echo ${ALLEL_LEN}`
                                FINAL_QUERY_LEN=`echo ${QUERY_LEN}`
				PREVIOUS_POS="${HIT_END}"

			fi
		done
		echo -e "${K}\t${FINAL_ID}\t100\tmulti\t${FINAL_TRUE_LEN}\t${FINAL_QUERY_LEN}" >> log.log

	elif  [ ${HITY} == 0 ]; then
		# Nie ma allelu o 100% pident z allelami w cgMLST
		# Iterujemy po calym pliku i szukamy takiego allelu o najnizszej sumie dlugosc mapowania - (mismatch + gaps)
		# Defaultowe wartosci gdyby analizowany plik byl pusty
		# Na koniec zwrocimy i tak -1, ale zostawiam sobie do testow dodatkowe informacje
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
		
		# echo -e "${K}\t${FINAL_ID}\t${FINAL_PIDENT}\tunk\t${FINAL_TRUE_LEN}\t${FINAL_QUERY_LEN}\t${FINAL_SCORE}" >> log.log
		# docelowo entero i tak daje tu brak allelu wiec aby byc spojnym damy -1
		# mozna dac kod jakis inny ...
		echo -e "${K}\t-1\t${FINAL_PIDENT}\tunk\t${FINAL_TRUE_LEN}\t${FINAL_QUERY_LEN}\t${FINAL_SCORE}" >> log.log
	fi

	rm blastn_tmp.tab
done
