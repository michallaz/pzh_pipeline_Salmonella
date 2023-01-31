#!/bin/bash

#na razie tylko komendy z etoki do dostania MLST, cgMLST, antygenow OH

F1=`ls *_R1*.fastq.gz`
F2=`ls *_R2*.fastq.gz`

FILENAME=`basename -s .fastq.gz ${F1}`
NANOPOER=`echo $1` # 1 - tak, 0 - nie

#metaphlan i kraken2 do analizy zanieczysczen
echo 'Metaphlan'
metaphlan ${F1},${F2} --bowtie2out metagenome.bowtie2.bz2 --nproc 24 --input_type fastq --unclassified_estimation --bowtie2db /bowtie_db/ --bt2_ps very-sensitive --min_mapq_val 30 --tax_lev 'a' -o metaphlan_output.txt 

echo 'Kraken2'
/opt/docker/EToKi/externals/kraken2 -db /kraken2_sdb --report raport_kraken2.txt --threads 24 --gzip-compressed --minimum-base-quality 30 --paired  --use-names ${F1} ${F2} >> raport_kraken2_individualreads.txt 2>&1

### Teraz jakis parser wynikow tych softow

### skladanie
echo 'Etoki budowanie genomu'
EToKi.py prepare --pe ${F1},${F2} -p prep_out
EToKi.py assemble --pe prep_out_L1_R1.fastq.gz,prep_out_L1_R2.fastq.gz --se prep_out_L1_SE.fastq.gz -p asm_out

### ocenianie genomu na podstawie tego co opisano w grancie i jest z defaultu w enterobase, czy genom spelnia zalozone kryteria

echo 'Etoki MLST 7 genomwy'
EToKi.py MLSType -i asm_out/etoki.fasta -r /Achtman7GeneMLST_entero/refset.fasta -k test_sample -o MLSTout.txt -d /enterobase/MLST_etoki/MLST_database.tab

# cgMLST
# sciezki do bazy cgMLST ustawiono wewnatrz podskryptu run_blastn_ver5.sh
echo 'cgMLST'
./run_blastn_ver5.sh asm_out/etoki.fasta

#EToKi.py MLSType -i asm_out/etoki.mapping.reference.fasta -r /enterobase/cgMLSTv2_etoki/refset.fasta -k test_sample -o cgMLSTout.txt -d /enterobase/cgMLSTv2_etoki/cgMLST_database.tab

#seqsero od zlozonego genomu
echo 'antygeny OH'
SeqSero2_package.py -m k -t 4 -p 4 -i asm_out/etoki.fasta -d seqsero_out
#sistr
sistr --qc -vv --alleles-output allele-results.json --novel-alleles novel-alleles.fasta --cgmlst-profiles cgmlst-profiles.csv -f tab -o sistr-output.tab asm_out/etoki.fasta

# odpalanie hiercc, w skrypcie prep_hierCC.py ustawiono sciezki na sztywno
echo 'okreslenie hiercc'
ST=`python prep_hierCC.py`
pHierCC -p /cgMLST2_entero/profiles.list -a /cgMLST2_entero/profile_output.npy -z to_hiercc.txt -o out_hierCC

#pointfinder
echo 'Pointfinder' ten wymaga pandas 1.0.5
python -m resfinder -o resfinder_out -s "s.enterica"  -l 0.6 -t 0.8 --acquired --point -k /opt/docker/kma/kma -db_disinf /opt/docker/disinfinder_db/ -db_res /opt/docker/resfinder_db/ -db_point /opt/docker/pointfinder_db/ -ifq ${FILENAME}*
#virulence
echo 'Analyzing virulence islands for Salmonella'
mkdir tmp
mkdir wyniki_spifinder
spifinder.py -i ${F1} ${F2} -o wyniki_spifinder -tmp tmp -mp kma -p /opt/docker/spifinder_db/ -x -t 0.95

# wyciaganie plasmidow 

# szukanie circular contigs

#skladanie z mapowaniem na genom referencyjny  i powtorzenie cgMLST ?
#TO DO to jest dosc proste potrzebny jest tylko genom referencyjny

# wyciaganie sekwencji bialek, budowanie struktur wybranych bialek ? tylko tych z mutacjami w genomie referencyjnym ?
#TO DO to tez jest proste wystarczy po mapowaniu na genom referencyjny wyciagnac nextstrainem 

