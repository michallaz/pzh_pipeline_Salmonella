params.cpus = 25
params.min_coverage_ratio = 0.1 
// Etoki odrzuca contigi ktorych srednie pokrycie jest mniejsze niz 0.2 globalnego pokrycia
// my nie jestesmy tak ostrzy dajemy 0.1 co tlumaczymy uzyciem innego alignera
// params.min_coverage_ratio_nanopore = 0.1 // Ten parametr nie jest moze  posluzy do nanopre'a

params.species = 's.enterica' // inne mozliwosci to c.jejuni c.coli e.coli . Gatungi do resfindera
params.reads = '/mnt/sda1/michall/Salmonella/*_R{1,2}_001.fastq.gz'

// eqluowane w main worflow jesli jest cos innego niz illumina
// to wchozi sciezka nanoporowe'a
params.machine = 'Illumina'

params.kraken2_db_absolute_path_on_host = "/home/michall/kraken2/kraken2_db/kraken2_sdb/" // na sztywno bo nie chce mi sie tego ustawiac za kazdym razem
params.quality_initial = 5 // Parametr stosowany aktualnie tylko przez krakena
params.Achtman7GeneMLST_db_absolute_path_on_host = "/home/michall/git/pzh_pipeline_Salmonella/Achtman7GeneMLST_entero" //ponownie na sztywno do poprawy docelowo 
params.cgMLST_db_absolute_path_on_host = "/mnt/sda1/michall/db/cgMLST_30042024" // sciezka do alleli z cgMLST + informacja o profilacj hierCC
params.enterobase_api_token = "eyJhbGciOiJIUzI1NiIsImlhdCI6MTcyMDQzNjQxMSwiZXhwIjoxNzM2MjA0NDExfQ.eyJfIjoibUsyNFlZSHd4SyIsInVzZXJuYW1lIjoiTWljaGFsX0xhem5pZXdza2kiLCJpZCI6ODg4MCwiYWRtaW5pc3RyYXRvciI6bnVsbCwiZW1haWwiOiJtbGF6bmlld3NraUBwemguZ292LnBsIiwiYXBpX2FjY2Vzc19jbG9zdHJpZGl1bSI6IlRydWUiLCJhcGlfYWNjZXNzX2Vjb2xpIjoiVHJ1ZSIsImFwaV9hY2Nlc3Nfc2VudGVyaWNhIjoiVHJ1ZSJ9.VEsyVPv8sn1zG7d3uFqEjfk6XFS2qP8P5Y5mh9VPE9w" // klucz api nadawany przez ENTEROBASE po rejestracji na stronie + wystapieniu o klucz 


// Kontenery uzywane w tym skrypcie 
// salmonella_illumina:2.0 - bazowy kontener z programami o kodem
// staphb/prokka:latest - kontener z prokka, w notatkach mam ze budowanie programu od 0 jest meczace bo to kod sprzed ponad 4 lat 
// infl_nanopore_eqa:2.18 - tu jest medaka



process check_etoki {
  // Testowa funkcja
  container  = 'salmonella_illumina:2.0'
  tag "Fixing fastq dla sample $x"
  input:
  tuple val(x), path(reads)
  output:
  stdout
  script:
  read_1 = reads[0]
  read_2 = reads[1]
  """
  python /opt/docker/EToKi/EToKi.py -h
  """


}

process view_output {
  // Testowa funkcja
  input:
  val(x)
  output:
  stdout
  script:
  """
  echo $x
  """
}

process run_fastqc {
  // Lekko zmodyfkowany modul z sars-a
  tag "fastqc for sample ${x}"
  container  = 'salmonella_illumina:2.0'
  publishDir "pipeline_wyniki/${x}/QC", mode: 'copy'
  input:
  tuple val(x), path(reads)
  output:
  tuple path("*forward_fastqc.txt"), path("*reverse_fastqc.txt")

  script:
  """
  QUALITY=5
  fastqc --format fastq \
         --threads ${params.cpus} \
         --memory 2024 \
         --extract \
         --delete \
         --outdir . \
         ${reads[0]} ${reads[1]}
    r1=\$(basename ${reads[0]} .fastq.gz)
    r2=\$(basename ${reads[1]} .fastq.gz)
    /data/parse_fastqc_output.py \${r1}_fastqc/fastqc_data.txt \${QUALITY} >> ${x}_forward_fastqc.txt
    /data/parse_fastqc_output.py \${r2}_fastqc/fastqc_data.txt \${QUALITY} >> ${x}_reverse_fastqc.txt
  """
}

process clean_fastq {
  // Prosta funkcja zaimplementowa w etoki do czyszczenia plikow fastq, ma na celu rozdzielenie odczytow ktore sa sparowane
  // od tych ktore pary nie maja, trimmowanie odczytow na podstawie jakosci, zmiana nazwy odczytow .
  // Ta funkcja generalnie mimikuje dzialanie trimmomatic-a. Wiec teoretycznie mozna uzyc jego
  // Default na trimming to quality 6

  container  = 'salmonella_illumina:2.0'
  tag "Fixing fastq dla sample $x"
  input:
  tuple val(x), path(reads)
  output:
  tuple val(x), path('prep_out_L1_R1.fastq.gz'), path('prep_out_L1_R2.fastq.gz'), emit: PE_path
  tuple val(x), path('prep_out_L1_SE.fastq.gz'), emit: SE_path
  tuple val(x), path('prep_out_L1_R1.fastq.gz'), path('prep_out_L1_R2.fastq.gz'), path('prep_out_L1_SE.fastq.gz'), emit: All_path
  script:
  read_1 = reads[0]
  read_2 = reads[1]
  """
  python /opt/docker/EToKi/EToKi.py prepare --pe ${read_1},${read_2} -p prep_out -c ${params.cpus}
  """
}


process spades {
  // Funkcja do odpalania spadesa
  // Powstaly plik fasta jest poprawiany bo 
  // Hapo-G nie akceptuje "." w plikach fasta
  cpus params.cpus
  container  = 'salmonella_illumina:2.0'
  tag "Spades dla sample $x"
  input:
  tuple val(x), path('R1.fastq.gz'), path('R2.fastq.gz'), path('SE.fastq.gz')
  output:
  tuple val(x),  path('scaffolds_fix.fasta')
  script:
  """
  #  for 150bp reads SPAdes uses k-mer sizes 21, 33, 55, 77
  #  We strongly recommend not to change -k parameter unless you are clearly aware about the effect.
  # --isolate This flag is highly recommended for high-coverage isolate and multi-cell Illumina data;
  # --careful Tries to reduce the number of mismatches and short indels. ale nie dziala z isolate
  #  wywolanie spades w ramach etoki to tylko definicja inputu + liczby procesorow
  # /opt/docker/EToKi/externals/spades.py -t 8 --pe-1 1 /Salomenlla/test_etoki/id_151_novel_etoki_with_comments/prep_out_L1_R1.fastq.gz --pe-2 1 /Salomenlla/test_etoki/id_151_novel_etoki_with_comments/prep_out_L1_R2.fastq.gz --pe-s 2 /Salomenlla/test_etoki/id_151_novel_etoki_with_comments/prep_out_L1_SE.fastq.gz -o spades

  python /opt/docker/EToKi/externals/spades.py --isolate -t ${task.cpus}  --pe-1 1 R1.fastq.gz --pe-2 1 R2.fastq.gz --pe-s 2 SE.fastq.gz  -o spades_manual
  cat spades_manual/scaffolds.fasta  | awk '{if(\$0 ~ />/) {split(\$0, ala, "_"); print ">NODE"ala[2]} else print \$0}' >> scaffolds_fix.fasta
  """


}


process bwa_paired {
  // Funkcja do mapowania odczytow Pair-end
  cpus params.cpus
  container  = 'salmonella_illumina:2.0'
  tag "RE-mapowanie PE dla sample $x"
  input:
  tuple val(x), path('genomic_fasta.fasta'),  path(read_1),  path(read_2)
  output:
  tuple val(x), path('mapowanie_bwa_PE.bam')
  script:
  """
  /opt/docker/EToKi/externals/bwa index genomic_fasta.fasta
  /opt/docker/EToKi/externals/bwa mem -t ${task.cpus} -T 30 genomic_fasta.fasta ${read_1} ${read_2} |  /opt/docker/EToKi/externals/samtools fixmate -m -@ ${task.cpus} - - | /opt/docker/EToKi/externals/samtools sort -@ ${task.cpus} -  |  /opt/docker/EToKi/externals/samtools markdup -r -@ ${task.cpus}  -O BAM - mapowanie_bwa_PE.bam
  """
}


process bwa_single {
  // Funkcja do mapowania odczytow Single-end na genom
  cpus params.cpus
  container  = 'salmonella_illumina:2.0'
  tag "RE-mapowanie SE dla sample $x"
  input:
  tuple val(x), path('genomic_fasta.fasta'), path(reads)
  output:
  tuple val(x), path('mapowanie_bwa_SE.bam')
  script:
  """
  /opt/docker/EToKi/externals/bwa index genomic_fasta.fasta
  /opt/docker/EToKi/externals/bwa mem -t ${task.cpus} -T 30 genomic_fasta.fasta ${reads} |  /opt/docker/EToKi/externals/samtools sort -@ ${task.cpus} -  |  /opt/docker/EToKi/externals/samtools markdup -r -@ ${task.cpus} -O BAM - mapowanie_bwa_SE.bam
  """
}

process merge_bams {
  // process do mergowania bam-ow
  // potrzebna jesli uzyje jednak uzyje Hapo-G
  
  container  = 'salmonella_illumina:2.0'
  tag "Merging bam files for sample $x"
  input:
  tuple val(x), path(bam1), path(bam2)
  output:
  tuple val(x), path('merged_bams.bam')
  script:
  """
  /opt/docker/EToKi/externals/samtools merge -f merged_bams.bam $bam1 $bam2
  """ 
}

process _run_pilon {
  // Poki co zastepujemy hapo-g pilonem 
  // MaxPilonRun params.pilon_run
  // To bez zsensu bo bamy trzeba remapowac na nowy genom a ni brac mapowania
  // na poprzedni scaffold  
  // Funkcja zostaje ale docelowo jest do USUNIECIA

  container  = 'salmonella_illumina:2.0' 
  // w kontenerze 2.0 dodalem pilona ale nie chce kasowac 1.0
  tag "Pilon for sample $x"
  publishDir "pilon_wyniki/${x}", mode: 'copy'
  input:
  tuple val(x), path(bam1), path(bam2), path('genomic_fasta.fasta')
  output:
  tuple val(x), path('last_pilon.fasta'), path('pilon*changes')
  script:
  """ 
  # To trzeba umiescic w petli i polishowac tyle razy ile jest parametr lub nie ma zmian
  /opt/docker/EToKi/externals/samtools index  $bam1
  /opt/docker/EToKi/externals/samtools index  $bam2
  RUN=1  
  NO_CHANGES=100
  while [[ \${RUN} -le ${params.pilon_run} && \${NO_CHANGES} -gt 1 ]]
  # warunek jest prosty przerwij petle gdy osiagniesz maksymalna ilosc iteracji polishowania
  # albo poprzednia iteracja zwrocila nie wiecej niz jedna poprawge do scaffold-u
  do
  mkdir pilon_bwa_\${RUN}
  if [ -e last_pilon.fasta ]; then
        # robie to tylko po to by nie nadpisywac genomic_fasta.fasta
   	java -jar /opt/docker/EToKi/externals/pilon.jar  --genome last_pilon.fasta --frags $bam1 --unpaired $bam2 --output pilon_polish --vcf --changes --outdir pilon_bwa_\${RUN}
  else
  	java -jar /opt/docker/EToKi/externals/pilon.jar  --genome genomic_fasta.fasta --frags $bam1 --unpaired $bam2 --output pilon_polish --vcf --changes --outdir pilon_bwa_\${RUN}
  fi

  cp  pilon_bwa_\${RUN}/pilon_polish.fasta last_pilon.fasta
  cp  pilon_bwa_\${RUN}/pilon_polish.changes pilon_polish_\${RUN}.changes
  NO_CHANGES=`cat pilon_bwa_\${RUN}/pilon_polish.changes | wc -l`
  RUN=`echo "\${RUN} + 1" | bc -l`
  done
  """ 

}
process run_pilon {
  // Poki co zastepujemy hapo-g pilonem
  // Dokladne uzycie pilona jest tu https://github.com/broadinstitute/pilon/wiki/Requirements-&-Usage
  container  = 'salmonella_illumina:2.0'
  // w kontenerze 2.0 dodalem pilona ale nie chce kasowac 1.0
  tag "Pilon for sample $x"
  // publishDir "pilon_wyniki/${x}", mode: 'copy'
  input:
  tuple val(x), path(bam1), path(bam2), path('genomic_fasta.fasta')
  output:
  tuple val(x), path('latest_pilon.fasta'), path('latest_pilon.changes'), emit: ALL
  tuple val(x), path('latest_pilon.fasta'), emit: ONLY_GENOME
  script:
  """
  # indeksacja bam-ow
  /opt/docker/EToKi/externals/samtools index  $bam1
  /opt/docker/EToKi/externals/samtools index  $bam2

  # wywolanie pilon-a
  java -jar /opt/docker/EToKi/externals/pilon.jar --threads ${params.cpus} --genome genomic_fasta.fasta --frags $bam1 --unpaired $bam2 --output pilon_polish --vcf --changes --outdir pilon_bwa
 
  # przygotowanie outputu
  # cat pilon_bwa/pilon_polish.fasta | awk  -v ALA=${x} 'BEGIN{OFS=""}; {if(\$0 ~ />/) print \$0,"_", ALA; else print \$0}' >> last_pilon.fasta
  cat pilon_bwa/pilon_polish.fasta >> latest_pilon.fasta
  cp  pilon_bwa/pilon_polish.changes latest_pilon.changes
  """
}


process run_coverage {
  // Liczenie pokrycia dla kazdego contiga polaczona z filtorwanie odczytow
  // Przy uzyciu NASZEGO skryptu
  container  = 'salmonella_illumina:2.0'
  tag "Coverage-based filtering for sample $x"
  publishDir "pipeline_wyniki/${x}", mode: 'copy'
  input:
  tuple val(x), path(bam1), path('genomic_fasta.fasta')
  output:
  tuple val(x), path('final_scaffold_filtered.fa'), emit: ONLY_GENOME
  // final_scaffold_filtered.fa to nazwa ustawiona NA SZTYWNO w skrypcie coverage_filter.py
  script:
  """
  /opt/docker/EToKi/externals/samtools index  $bam1
  python  /opt/docker/EToKi/externals/coverage_filter.py genomic_fasta.fasta $bam1 ${params.min_coverage_ratio}
  """
}



process extract_final_stats {
// skopiowac z etoki
  container  = 'salmonella_illumina:2.0'
  tag "Calculating basic statistics for sample $x"
  publishDir "pipeline_wyniki/${x}", mode: 'copy'
  input:
  tuple val(x), path(fasta)
  output:
  tuple val(x), path('Summary_statistics.txt')
  // Dla zabawy napiszemy kod tutaj a nie w skrypcie
  // Uwzglednienie wciec to morederstwo, przerzucam do skryptu
  script:
  """
  python  /opt/docker/EToKi/externals/calculate_stats.py $fasta
  """
}

process run_7MLST {
  // wykorzystujemy bezposrednio Etoki, 7-genomy MLST w tym nie da sie pomylic
  container  = 'salmonella_illumina:2.0'
  containerOptions "--volume ${params.Achtman7GeneMLST_db_absolute_path_on_host}:/Achtman7GeneMLST_entero"
  tag "Predicting MLST for sample $x"
  publishDir "pipeline_wyniki/${x}", mode: 'copy'
  input:
  tuple val(x), path(fasta)
  output:
  tuple val(x), path('MLSTout.txt')
  script:
  """
  /opt/docker/EToKi/EToKi.py MLSType -i $fasta -r /Achtman7GeneMLST_entero/references.fasta -k ${x} -o MLSTout.txt -d /Achtman7GeneMLST_entero/MLST_database.tab
  """
}


process parse_7MLST {
  // Funkcja do parsowania wynikow run_7MLST
  // Zwraca plik w ktorym mamy 9 kolumn
  // Pierwsza to Sequnce type zgodnie z plikiem z profilem
  // Kolejen 7 kolumn to wersje alleli dla 7 genow
  // Kolumna 9 to odleglosc znalezionych wersji alleli do wskazanego ST, jesli mamy pelna zgodnosc to wypiujemy 0 i jest jedno-wierszowy wpis
  // Jesli nie ma profilu z pelna zgodnoasci to wpisujemy wszystkie ST ktore maja wskazana odleglosc / pewnie bedzie to max odleglosc 1, czyli roznica w jednym z alleli /
  // proces nie jest czescie run_7MLST bo korzystam z moich pythonowych skryptow
  
  container  = 'salmonella_illumina:2.0'
  containerOptions "--volume ${params.Achtman7GeneMLST_db_absolute_path_on_host}:/Achtman7GeneMLST_entero"
  tag "Pasring MLST for sample $x"
  publishDir "pipeline_wyniki/${x}", mode: 'copy', pattern: 'parsed_7MLST.txt'

  input:
  tuple val(x), path('MLSTout.txt')
  output:
  tuple val(x), path('parsed_7MLST.txt')
  script:
"""
#!/usr/bin/python
import  sys
sys.path.append('/data')
from all_functions_salmonella import *

known_profiles, klucze = create_profile('/Achtman7GeneMLST_entero/MLST_profiles.txt')
identified_profile = parse_MLST_fasta('MLSTout.txt')
matching_profile = getST(MLSTout = identified_profile, \
                         profile_dict = known_profiles, \
                         lista_kluczy = klucze)

  

if isinstance(matching_profile, str) and int(matching_profile) != 0:
	# we identified only one matching profile 
	formatted_string = '{aroC}\\t{dnaN}\\t{hemD}\\t{hisD}\\t{purE}\\t{sucA}\\t{thrA}'.format(**identified_profile)
	with open('parsed_7MLST.txt', 'w') as f:
		f.write(f'ST\\taroC\\tdnaN\\themD\\thisD\\tpurE\\tsucA\\tthrA\\tDistance\\n')
		f.write(f'{matching_profile}\\t{formatted_string}\\t0\\n')
elif isinstance(matching_profile, list):
	# no matching profile in the database we print all profiles with distance 1
	formatted_string = '{aroC}\\t{dnaN}\\t{hemD}\\t{hisD}\\t{purE}\\t{sucA}\\t{thrA}'.format(**identified_profile)
	with open('parsed_7MLST.txt', 'w') as f:
		f.write(f'ST\\taroC\\tdnaN\\themD\\thisD\\tpurE\\tsucA\\tthrA\\tDistance\\n')
		f.write(f'Novel\\t{formatted_string}\\t0\\n')
		for hit in matching_profile:
			formatted_string = '{aroC}\\t{dnaN}\\t{hemD}\\t{hisD}\\t{purE}\\t{sucA}\\t{thrA}'.format(**known_profiles[hit])
			f.write(f'{hit}\\t{formatted_string}\\t1\\n')
		
else:
	#No matching profiles with distance 1
	with open('parsed_7MLST.txt', 'w') as f:
		f.write(f'ST\\taroC\\tdnaN\\themD\\thisD\\tpurE\\tsucA\\tthrA\\tDistance\\n')
		f.write(f'Novel\\t{formatted_string}\\t0\\n')
		f.write(f'No matching profiles with at most 1 allelic difference')

"""
}

process run_Seqsero {
  container  = 'salmonella_illumina:2.0'
  tag "Predicting OH for sample $x with Seqsero"
  publishDir "pipeline_wyniki/${x}", mode: 'copy'
  cpus params.cpus
  input:
  tuple val(x), path(fasta)
  output:
  tuple val(x), path('seqsero_out/SeqSero_result.txt')
  script:
  """
  # -m to rodzaj algorytmu -m to chyba opart o k-mery
  # -t 4 to informacja ze inputem sa contigi z genomem
  # -p to procki 
  python /opt/docker/SeqSero2/bin/SeqSero2_package.py -m k -t 4 -p ${task.cpus} -i $fasta -d seqsero_out
  """
}

process run_sistr {
  container  = 'salmonella_illumina:2.0'
  tag "Predicting OH for sample $x with Sistr"
  publishDir "pipeline_wyniki/${x}", mode: 'copy'
  cpus params.cpus
  input:
  tuple val(x), path(fasta)
  output:
  tuple val(x), path('sistr-output.tab')
  script:
  // Uwaga sistr korzysta z WLASNEJ BAZY do pobrania z 
  // SISTR_DB_URL = 'https://sairidapublic.blob.core.windows.net/downloads/sistr/database/SISTR_V_1.1_db.tar.gz'
  // adres ustawiony jest na sztywno w /usr/local/lib/python3.8/dist-packages/sistr/src/serovar_prediction/constants.py
  // Nie wiem co bedzie w przyszlych wersjachh sistr i nie wiem jak to jest aktualizowane
  // baza rozpakowywana jest do 
  // /usr/local/lib/python3.8/dist-packages/sistr/data
  // NIE ZNALAZLEM opcji podania innej sciezki do bazy niz ten default w opcjach programu sistr
  // Baza sciagana jest na poziomie dockera
  // jesli baza nie jest sciagnieta program robi to sam przy wywolaniu sistr-a
  // baza ma jakies 50 Mb /przed rozpakowanie/ wiecjej uzywanie to nie jest jaki wielki problem 
  // ponadto w katalogu /usr/local/lib/python3.8/dist-packages/sistr/ musi byc plik dbstatus.txt
  // bo jak go nie ma to sistr_cmd dalej wymusza sciaganie bazy
  """
  /usr/local/bin/sistr --qc -vv --alleles-output allele-results.json --novel-alleles novel-alleles.fasta --cgmlst-profiles cgmlst-profiles.csv -f tab -t ${task.cpus} -o sistr-output.tab $fasta
  """
}

process run_pointfinder {
  container  = 'salmonella_illumina:2.0'
  tag "Predicting microbial resistance for sample $x"
  publishDir "pipeline_wyniki/${x}/pointfinder", mode: 'copy', pattern: "resfinder_out/ResFinder_results_table.txt"
  publishDir "pipeline_wyniki/${x}/pointfinder", mode: 'copy', pattern: "resfinder_out/pheno_table_*.txt"
  publishDir "pipeline_wyniki/${x}/pointfinder", mode: 'copy', pattern: "resfinder_out/PointFinder_results.txt"
  // pattern mowi co kopiowac
  input:
  tuple val(x), path(fasta)
  output:
  // tak naprawde pheno_table jest istotne ale pozostale pliki tlumacza jakie geny niosace opornosci znaleziono i jakie mutacje
  // w genach wywolujace opornosc znaleziono
  tuple val(x), path('resfinder_out/pheno_table*.txt'), path('resfinder_out/ResFinder_results_table.txt'), path('resfinder_out/PointFinder_results.txt')
  script:
  // resfinder rozumie 4 organizmy
  // campylobacter jejuni
  // campylobacter coli
  // escherichia coli
  // salmonella enterica
  // inne opcje to
  // -l 0.6 Minimum (breadth-of) coverage of ResFinder within the range 0-1.
  // -r 0.8 Threshold for identity of ResFinder within the range 0-1
  // oba parametry decyduja jaki procent alignmentu i z jaka identycznoscia musi miec nasz sample w stosunku 
  // do odczytow z bazy aby uznac ze dany gen wystepuje w probce 
  // --acquired Run resfinder for acquired resistance genes, czyli szukaj znanych genow wywolujacych opornosci
  // --point Run pointfinder for chromosomal mutations, szukamy w konkretnych genach mutacji punktowych
  // odpowiedzialnych za nabycie konkretnych opornosci
  // pozostale opcje to sciezki do baz /instalowanych wraz z tworzeniem obrazu/ 
  """
  # resfinder-owi mozna tez podac pliki fastq (-ifq) jako input
  python -m resfinder -o resfinder_out/ -s ${params.species}  -l 0.6 -t 0.8 --acquired --point -k /opt/docker/kma/kma -db_disinf /opt/docker/disinfinder_db/ -db_res /opt/docker/resfinder_db/ -db_point /opt/docker/pointfinder_db/ -ifa ${fasta}
  #virulence
  # echo 'Analyzing virulence islands for Salmonella'
  # mkdir tmp
  # mkdir wyniki_spifinder
  # wybralem opcje jak na stronie (inne niz default w skrypcoe pythonnowym)
  # spifinder.py -i \${F1} \${F2} -o wyniki_spifinder -tmp tmp -mp kma -p /opt/docker/spifinder_db/ -x -t 0.95 -l 0.6
 """
}

process run_cgMLST {
  container  = 'salmonella_illumina:2.0'
  // podmontujemy z zewnatrz cgMLST, montujemy na /cgMLST2_entero bo skrypt run_blastn_ver6.sh
  // ma tak za hard-kodowane i nie chce tego zmieniac poki co
  containerOptions "--volume ${params.cgMLST_db_absolute_path_on_host}:/cgMLST2_entero"
  tag "Predicting cgMLST for sample $x"
  publishDir "pipeline_wyniki/${x}", mode: 'copy'
  cpus params.cpus
  input:
  tuple val(x), path(fasta)
  output: 
  tuple val(x), path('cgMLST.txt')
  // output to na razie plik ktory zwraca informacja jaka wersja allelu jest z informacja czy 
  // jaki byl procent identycznosci sekwencyjne j mapowania na ten allel i czy ewentualnie mapowanie obejmowala caly allel, jego czesc czy moze bylo wiele mapowan na ten allel 
  // output nie jest wiec kompatybilny z pHierCC, potem to POPRAWIC

  // Aha run_blastn_ver6.sh wykorzystuje pod spodem xargs zeby rownolegle puszczac max ${task.cpus} 1-procesowych blastow
  script:
  """
  /data/run_blastn_ver10.sh $fasta ${task.cpus}
  cat log.log | cut -f1,2 > cgMLST.txt
  """
}

process parse_cgMLST {
  // Funkcja do parsowania wynikow run_cgMLST
  // Zwraca plik parsed_cgMLST.txt ktory albo zawiera pojedyncza cyfre z infromrmacja o pasujacym profilu
  // "Multiple_distance1" jesli jest wiele ST z odlegloscia 1
  // lub wiele linijek w trypi ST:odleglosc jesli nie bylo ST z odlegloscia 0 lub 1
  // UWAGA wczytanie wszystkich 400k profili wymaga ponad 50 Gb Ram-u poki co
  // wiec lepiej nie puszczac tego dla wiecej niz 4 probek na raz

  container  = 'salmonella_illumina:2.0'
  containerOptions '--volume /mnt/sda1/michall/db/cgMLST_30042024:/cgMLST2_entero'
  tag "Parsing cgMLST for sample $x"
  publishDir "pipeline_wyniki/${x}", mode: 'copy', pattern: 'parsed_cgMLST.txt'
  maxForks 4
  input:
  tuple val(x), path('cgMLST.txt')
  output:
  tuple val(x), path('parsed_cgMLST.txt')
  script:
"""
#!/usr/bin/python
import  sys
sys.path.append('/data')
from all_functions_salmonella import *

known_profiles, klucze = create_profile('/cgMLST2_entero/profiles.list')
identified_profile = parse_MLST_blastn('cgMLST.txt')
matching_profile = getST(MLSTout = identified_profile, \
			profile_dict = known_profiles, \
			lista_kluczy = klucze)



if isinstance(matching_profile, str) and int(matching_profile) != 0:
	# we identified only one matching profile
	with open('parsed_cgMLST.txt', 'w') as f:
		f.write(f'{matching_profile}\\t0\\n')
elif isinstance(matching_profile, list):
	# no matching profile in the database we print all profiles with distance 1
	with open('parsed_cgMLST.txt', 'w') as f:
		f.write(f'Multiple_distance1\\t1\\n')

elif isinstance(matching_profile, dict):
	#No matching profiles with distance 1
	with open('parsed_cgMLST.txt', 'w') as f:
		for klucz,wartosc in matching_profile.items():
			f.write(f'{klucz}:{wartosc}\\n')

"""
}


process run_prokka {
  // Propkka to soft do przewidywania genow w genomie (zwraca tez sekwencje bialek i annotace)
  // Ktory pod spodem korzysta z prodigal to przewiddywnia genow
    // Komentarz czemu uzywa prodigal 2.6.3 a nie 3.0 
    // Prodigal sluzy do predykcji genow w genomie
    // Uzywamy wersji 2.6.3 choc na ich wiki caly czas mowia o 3.0, ktorej nie moge znalezc
    // Moze ma to zwiazek z issues w ktorych facet pisze ze 3.0 to work in progress 
    // https://github.com/hyattpd/Prodigal/issues/45
  
  // instalacja prokka jest meczaca ale korzystamy z wystawionego przez autorow obrazu 


  // opcja --metagenome w prokka jest wywolaniem prodigal z opcja -p meta
  // opcja meta in which Prodigal applies pre-calculated training files to the provided input sequence and predicts genes based on the best results.
  //  i winnym miejscu
  // Isolated, short sequences (<100kbp) such as plasmids, phages, and viruses should generally be analyzed using meta mode
  // Testowe puszczenie z opcja meta i bez niej nawet na dosc stabilnym genomie jak sample 151 daje troche inne wyniki

  // opcja --compliant Force Genbank/ENA/DDJB compliance: --addgenes --mincontiglen 200 --centre XXX (default OFF)
  // opcja --kingdom Bacteria jest defaultem zaklada jaki typ genomu analizujemy

  container  = 'staphb/prokka:latest'
  tag "Predicting genes for sample $x"
  publishDir "pipeline_wyniki/${x}", mode: 'copy'
  cpus params.cpus
  input:
  tuple val(x), path(fasta)
  output:
  tuple val(x), path('prokka_out/prokka_out*gff'), path('prokka_out/*faa'), path('prokka_out/*.ffn'), path('prokka_out/*.tsv')
  script:
  """
  prokka --metagenome --cpus ${params.cpus} --outdir prokka_out --prefix prokka_out --compliant --kingdom Bacteria $fasta 
  """
}

process run_VFDB {
  // Baza z czynnikami wirulencji w roznych bakteriach w tym salmonelli
  // Przygotowana baza z instrukcja jak ja przygotowac jest w /mnt/sda1/michall/db/VFDB/README
  container  = 'salmonella_illumina:2.0'
  containerOptions '--volume /mnt/sda1/michall/db/VFDB:/db'
  // Ponownie montujemy na sztyno do /db bo taka lokalizacje na sztywno ma wpisany moj skrypt
  tag "Predicting VirulenceFactors for sample $x"
  publishDir "pipeline_wyniki/${x}", mode: 'copy'
  cpus params.cpus
  input:
  // inputem jest output procesu run_prokka
  tuple val(x), path(gff), path(faa), path(ffn), path(tsv)
  output:
  tuple val(x), path('VFDB_summary.txt')
  script:
  """
  SPEC='Salmonella' # wpradzie jest juz paramter params.species, ale baza VFDB go nie zrozumie.
  PIDENT=80 # minimalna identycznosc sekwencyjna aby stwierdzic ze jest hit 
  COV=80 # minimalne pokrycie query i hitu aby stwierdzic ze jest hit 
  EVAL=0.01  # maksymalne e-value
  /opt/docker/EToKi/externals/run_VFDB.sh $ffn ${task.cpus} \${SPEC} \${PIDENT} \${EVAL} \${COV}
  """

}



process run_spifinder {
  container  = 'salmonella_illumina:2.0'
  tag "Predicting virulence islands for sample $x"
  publishDir "pipeline_wyniki/${x}", mode: 'copy'
  cpus params.cpus
  input:
  tuple val(x), path(fasta)
  output:
  tuple val(x), path('spifinder_results/*')

  script:
  """
  mkdir spifinder_results # program wymaga tworzenia katalogu samodzielnie
  # -mp opcja jak szukamy wysp kma jest gdy inputem sa dane surowe, blastn gdy podajemy zlozony genom/contigi
  # -p sciezka do bazy, sciagana w trakcie budowy kontenera
  # -l i -t to parametrty na alignment coverage i seq id
  # -x to rozszerzony output
  python /opt/docker/spifinder/spifinder.py -i $fasta -o spifinder_results -mp blastn -p /opt/docker/spifinder_db/ -l 0.6 -t 0.9 -x
  """
}


process run_kraken2_illumina {
  // modul wziety z pipeline do SARS
  // poprawilem tylko pare rzeczy zwiazane z parametrami i kontenerami ktore u mnie sa inne
  // kraken2 instalowany jest przez ETOKI i jest w path kontenera do salmonelli
  tag "kraken2:${x}"
  container  = 'salmonella_illumina:2.0'
  publishDir "pipeline_wyniki/${x}/kraken2", mode: 'copy', pattern: "report_kraken2_individualreads.txt"
  publishDir "pipeline_wyniki/${x}/kraken2", mode: 'copy', pattern: "report_kraken2.txt"
  publishDir "pipeline_wyniki/${x}", mode: 'copy', pattern: "summary_kraken.txt"
  containerOptions "--volume ${params.kraken2_db_absolute_path_on_host}:/home/external_databases/kraken2"
  maxForks 5

  input:
  tuple val(x), path(reads)

  output:
  tuple val(x), path('report_kraken2.txt'), path('report_kraken2_individualreads.txt'), path('summary_kraken.txt')

  script:
  """
  kraken2 --db /home/external_databases/kraken2 \
          --report report_kraken2.txt \
          --threads ${params.cpus} \
          --gzip-compressed \
          --minimum-base-quality ${params.quality_initial} \
          --use-names ${reads[0]} ${reads[1]} >> report_kraken2_individualreads.txt 2>&1
  # parse kraken extract two most abundant FAMILIES
  LEVEL="G" # G to chyba rodzzaj (genus ?)  S to pewnie gatunek (SPECIES)
  SPEC1=`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -1 | tr -s " " | cut -f6 | tr -d "="`
  SPEC2=`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f6 | tr -d "="`
  ILE1==`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}'| grep -w \${LEVEL} | sort -rnk 1 | head -1 | tr -s " " | cut -f1 | tr -d " "`
  ILE2==`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}'| grep -w \${LEVEL} | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f1 | tr -d " "`

  echo -e "${x}\t\${SPEC1}\${ILE1}%\t\${SPEC2}\${ILE2}%" >> summary_kraken.txt
  """
}

process run_pHierCC {
  
  // Funkcja odpytuje API Enterobase w celu wyciagniecia profuili z tej bazu
  // Ponadto Funkcja odpytuje dwa pliki przygotowane przeze mnie 
  // Jeden zawierajacy klastrowanie SINGLE linkage zbudowane na 430k profili z enterobase
  // Drugi zawierajacy klastrowanie COMPLETE linkage zbudowane na 430k profili z enterobae
  // Jesli ST jest nowy NIE obecny w zadnej z baz program nic nie zwraca ?
  
  container  = 'salmonella_illumina:2.0'
  containerOptions "--volume ${params.cgMLST_db_absolute_path_on_host}:/cgMLST2_entero"
  tag "Predicting hierCC levels for sample $x"
  publishDir "pipeline_wyniki/${x}/phiercc", mode: 'copy'
  input:
  tuple val(x), path('parsed_cgMLST.txt')
  output:
  tuple val(x), path('parsed_phiercc_enterobase.txt'), path('parsed_phiercc_minimum_spanning_tree.txt'), path('parsed_phiercc_maximum_spanning_tree.txt')
  script:
"""
#!/usr/bin/python
### kod pobrany ze strony https://enterobase.readthedocs.io/en/latest/api/api-getting-started.html ###

import  sys
from urllib.request import urlopen
from urllib.error import HTTPError
import urllib
import base64
import json
import gzip

API_TOKEN = "${params.enterobase_api_token}"

def __create_request(request_str):
    base64string = base64.b64encode('{0}: '.format(API_TOKEN).encode('utf-8'))
    headers = {"Authorization": "Basic {0}".format(base64string.decode())}
    request = urllib.request.Request(request_str, None, headers)
    return request

def getST(my_file):
    # prosta funkcja zwraca wszystkie numery ST z odlegloscia 0 (patrz funkcja parse_cgMLST)
    my_ST = []
    with open(my_file) as f:
        for line in f:
            line = line.rsplit()
            if int(line[1]) == 0:
                my_ST.append(line[0])
    return my_ST

my_ST = getST('parsed_cgMLST.txt')
if len(my_ST) == 0:
    # nie znaleziono ST inicjalizuje puste pliki
    with open('parsed_phiercc_enterobase.txt', 'w') as f:
        f.write('ST\\tHC0\\tHC2\\tHC5\\tHC10\\tHC20\\tHC50\\tHC100\\tHC200\\tHC400\\tHC900(ceBG)\\tHC2000\\tHC2600\\tHC2850(subsp.)\\n')
        f.write(f'UNK\\tUNK\\n')

    with open('parsed_phiercc_minimum_spanning_tree.txt', 'w') as f:
        f.write('ST\\tHC0\\tHC2\\tHC5\\tHC10\\tHC20\\tHC50\\tHC100\\tHC200\\tHC400\\tHC900(ceBG)\\tHC2000\\tHC2600\\tHC2850(subsp.)\\n')
        f.write(f'UNK\\tUNK\\n')
 
    with open('parsed_phiercc_maximum_spanning_tree.txt', 'w') as f:
        f.write('ST\\tHC0\\tHC2\\tHC5\\tHC10\\tHC20\\tHC50\\tHC100\\tHC200\\tHC400\\tHC900(ceBG)\\tHC2000\\tHC2600\\tHC2850(subsp.)\\n')
        f.write(f'UNK\\tUNK\\n')
    
# wyciaganie danych z enterobase	
with open('parsed_phiercc_enterobase.txt', 'w') as f:
    f.write('ST\\tHC0\\tHC2\\tHC5\\tHC10\\tHC20\\tHC50\\tHC100\\tHC200\\tHC400\\tHC900(ceBG)\\tHC2000\\tHC2600\\tHC2850(subsp.)\\n')
    for ST in my_ST:
        address = f"https://enterobase.warwick.ac.uk/api/v2.0/senterica/cgMLST_v2/sts?st_id={ST}&scheme=cgMLST_v2&limit=5"
        try:
            response = urlopen(__create_request(address))
            data = json.load(response)
            #data['STs'][0]['info']['hierCC']
            formatted_string = '{d0}\\t{d2}\\t{d5}\\t{d10}\\t{d20}\\t{d50}\\t{d100}\\t{d200}\\t{d400}\\t{d900}\\t{d2000}\\t{d2600}\\t{d2850}'.format(**data['STs'][0]['info']['hierCC'])
            f.write(f'{ST}\\t{formatted_string}\\n')
        except:
            f.write(f'{ST}\\tUNK\\n')
            

# wyciaganie danych z moich przeliczonych wczesniej plikow
with open('parsed_phiercc_minimum_spanning_tree.txt', 'w') as f:
    f.write('ST\\tHC0\\tHC2\\tHC5\\tHC10\\tHC20\\tHC50\\tHC100\\tHC200\\tHC400\\tHC900(ceBG)\\tHC2000\\tHC2600\\tHC2850(subsp.)\\n')
    for ST in my_ST:
        with gzip.open('/cgMLST2_entero/profile_single_linkage.HierCC.gz') as f2:
            for line in f2:
                line = list(map(lambda x: x.decode('utf-8', errors='replace'), line.split()))
                if line[0] == ST:
                    f.write(f'{line[0]}\\t{line[1]}\\t{line[3]}\\t{line[6]}\\t{line[11]}\\t{line[21]}\\t{line[51]}\\t{line[101]}\\t{line[201]}\\t{line[401]}\\t{line[901]}\\t{line[2001]}\\t{line[2601]}\\t{line[2851]}\\n')
                    break # konczymy iterowac pl pliku
        
with open('parsed_phiercc_maximum_spanning_tree.txt', 'w') as f:
    f.write('ST\\tHC0\\tHC2\\tHC5\\tHC10\\tHC20\\tHC50\\tHC100\\tHC200\\tHC400\\tHC900(ceBG)\\tHC2000\\tHC2600\\tHC2850(subsp.)\\n')
    for ST in my_ST:
        with gzip.open('/cgMLST2_entero/profile_complete_linkage.HierCC.gz') as f2:
            for line in f2:
                line = list(map(lambda x: x.decode('utf-8', errors='replace'), line.split()))
                if line[0] == ST:
                    f.write(f'{line[0]}\\t{line[1]}\\t{line[3]}\\t{line[6]}\\t{line[11]}\\t{line[21]}\\t{line[51]}\\t{line[101]}\\t{line[201]}\\t{line[401]}\\t{line[901]}\\t{line[2001]}\\t{line[2601]}\\t{line[2851]}\\n')
                    break # konczymy iterowac pl pliku

"""
}

process extract_historical_data {
  // W przypadku znalezienia ST ktory jest  bazie enterobase
  // Wyciagamy informacje gdzie i kiedy byly zebrane probki z tym ST
  // Generalnie rozszerzenie o pytania o ST na konkretnym poziomie hierCC
  // mozna latwo doimplementowac, choc wymagaja dodatkowe odpytania bazy w tablei 'STs'
  // tabela straindata zawierajaca informacje o ST nie zawiera informacji o phierCC 
  container  = 'salmonella_illumina:2.0'
  tag "Extracting historical data for sample $x"
  publishDir "pipeline_wyniki/${x}/phiercc", mode: 'copy'
  input:
  tuple val(x), path('parsed_cgMLST.txt')
  output:
  tuple val(x), path('enterobase_historical_data.txt')
  script:
"""
#!/usr/bin/python

import  sys
from urllib.request import urlopen
from urllib.error import HTTPError
import urllib
import base64
import json


API_TOKEN = "${params.enterobase_api_token}"

def __create_request(request_str):
    base64string = base64.b64encode('{0}: '.format(API_TOKEN).encode('utf-8'))
    headers = {"Authorization": "Basic {0}".format(base64string.decode())}
    request = urllib.request.Request(request_str, None, headers)
    return request

def getST(my_file):
    # prosta funkcja zwraca wszystkie numery ST z odlegloscia 0 (patrz funkcja parse_cgMLST)
    my_ST = []
    with open(my_file) as f:
        for line in f:
            line = line.rsplit()
            if int(line[1]) == 0:
                my_ST.append(line[0])
    return my_ST

my_ST = getST('parsed_cgMLST.txt')


# downloading entrire strain table 
address = 'https://enterobase.warwick.ac.uk/api/v2.0/senterica/strains?my_strains=false&sortorder=asc&return_all=true&offset=0'
response = urlopen(__create_request(address))
data = json.load(response)

# in dict - keys - strain name; values - 2- element list with countr and collection year 
strain_info= {}
for strain in data['Strains']:
    strain_info[strain['strain_barcode']] = []
    strain_info[strain['strain_barcode']].append(strain['country'])
    strain_info[strain['strain_barcode']].append(strain['collection_year'])

# querying "strain" table in f'{step}'-sized chunks 
# 150 seems to be optimal number , larger quieries are rejected
start = 0
step = 150
orignal_list = list(strain_info.keys())
for end in range(step,len(orignal_list), step):
    # create chunls
    #print(end) # for tests
    list_of_ids = orignal_list[start:end] # extract barcodes in chunks
    # build url link
    address2 = f"https://enterobase.warwick.ac.uk/api/v2.0/senterica/straindata?limit={step}&sortorder=asc&" + ('barcode={}&'*step).format(*list_of_ids) + "offset=0"	
    response2 = urlopen(__create_request(address2))
    data2 = json.load(response2)
    # parse output basically we remove data from strain_info if the do not have required Sequence type
    for klucz in data2['straindata']:
        to_remove = 1
        for scheme in data2['straindata'][klucz]['sts']:
            if 'cgMLST_v2' in scheme.values() and scheme['st_id'] == my_ST:
                strain_info[klucz].append(scheme['st_id'])
                to_remove = 0
        if to_remove:
            try:
                del(strain_info[klucz])
            except:
                print(f'No such key {klucz}')
    start = end

# na tym etapie strain_info powinno zawierac tylko recody z naszym ST
with open('enterobase_historical_data.txt', 'w') as f:
    f.write(f'Country\\tYear\\tST\\n'
    for klucz,wartosc in strain_info.items():
        f.write(f'{wartosc[0]}\\t{wartosc[1]}\\t{wartosc[2]}\\n'
"""
}


//process build_model {
  // Ogolny proces do budowania modelu
  // Zostaje do implementacji przez Michala K.
//}


process run_plasmidfinder {
  container  = 'salmonella_illumina:2.0'
  tag "Predicting plasmids for sample $x"
  publishDir "pipeline_wyniki/${x}/plasmidfinder_results", mode: 'copy'
  cpus params.cpus
  input:
  tuple val(x), path(fasta)
  output:
  tuple val(x), path('plasmidfinder_results/*')

  script:
  """
  # -i to oczywiscie input na podstawie jego rozszerzenia program wybiera metode do analizy (kma dla fastq i blastn dla fasta)
  # -o to katalog z wynikami, musi istniec
  # -p to sciezka do katalog z bazami (tymi z repo plasmidfinder_db)
  # -l to minimalny procent sekwencji jaki musi alignowac sie na genom
  # -t to minimalne sequence identity miedzy query a subject
  # -x printuj dodatkowe dane w output (alignmenty)
  mkdir plasmidfinder_results # program wymaga tworzenia katalogu samodzielnie
  /opt/docker/plasmidfinder/plasmidfinder.py  -i $fasta -o plasmidfinder_results -p /opt/docker/plasmidfinder_db -l 0.6 -t 0.9 -x
  """
}

// FUNKCJE DODANE DLA ANALIZY NANOPORE //

process run_flye {
  // nano-raw to input dla wersji przed guppym 5+ z duza liczba bledow
  // -g to estymowana wielkosc ocekiwanego genomu
  // -t to liczba watkow CPU
  // -i to liczba powtorzec wygladzania genomy
  // --no-alt-contig to informacja aby nie podawac haplotypow gdyby byly, w koncu salmonella to monoploid ...

  container  = 'salmonella_illumina:2.0'
  tag "Predicting scaffold with flye for sample $x"
  publishDir "pipeline_wyniki/${x}", mode: 'copy', pattern: "*"
  cpus params.cpus
  input:
  tuple val(x), path(fastq_gz)
  output:
  tuple val(x), path('output/assembly.fasta') 
  
  script:
  """
  # /data/Flye to sciezka z Flye instalowanego z github, uwaga
  # w kontenerze tez jest flye instalowant przez etoki i ten jest w PATH
  /data/Flye/bin/flye --nano-raw ${fastq_gz} -g 6m -o output -t ${params.cpus} -i 3 --no-alt-contig

  """
 
}

process clean_fastq_SE {
  // Prosta funkcja zaimplementowa w etoki do czyszczenia plikow fastq, wersja dla Single End
  // wiec de facto pozostawiam tylko quality trim
  // Zmieniono quality trimming na 5
  container  = 'salmonella_illumina:2.0'
  tag "Fixing fastq dla sample $x"
  input:
  tuple val(x), path(read)
  output:
  tuple val(x), path('prep_out_L1_SE.fastq.gz')
  script:
  """
  python /opt/docker/EToKi/EToKi.py prepare --se ${read} -c ${params.cpus} -q 5 -p prep_out
  """
}



process run_minimap2 {
  // Proces do mapowania odczytow na scaffold
  tag "Remapping of reads to predicted scaffold for sample $x"
  container  = 'salmonella_illumina:2.0' 
  publishDir "pipeline_wyniki/${x}", mode: 'copy', pattern: "*"
  input:
  tuple val(x), path(fasta), path(reads) 
  output:
  tuple val(x), path('sorted.bam'), path(fasta)
  script:
  """
  # minmap jest zarowno w PATH z etoki/externals jak i w /data/Flye/bin

  minimap2 -a -x map-ont -t ${params.cpus} $fasta $reads | samtools view -bS -F 2052 - | samtools sort -@ ${params.cpus} -o sorted.bam -
  
  """
}

process run_minimap2_2nd {
  // Proces do mapowania odczytow na scaffold
  // W nanopre w jednym workflow uzywam go 2 razy wiec musze zrobic ta glupia kopie
  tag "Remapping of reads to predicted scaffold for sample $x"
  container  = 'salmonella_illumina:2.0'
  publishDir "pipeline_wyniki/${x}", mode: 'copy', pattern: "*"
  input:
  tuple val(x), path(fasta), path(reads)
  output:
  tuple val(x), path('sorted.bam'), path(fasta)
  script:
  """
  # minmap jest zarowno w PATH z etoki/externals jak i w /data/Flye/bin

  minimap2 -a -x map-ont -t ${params.cpus} $fasta $reads | samtools view -bS -F 2052 - | samtools sort -@ ${params.cpus} -o sorted.bam -

  """
}

process run_pilon_nanopore {
  // De facto slave run_pilon, ale akcpetujemy jeden bam
  // Do testow nie wiem jak program zachowa sie na danych nanopre ktore sa kiepskie
  // Moze wykladzac medaka ?

  container  = 'salmonella_illumina:2.0'
  // w kontenerze 2.0 dodalem pilona ale nie chce kasowac 1.0
  tag "Pilon for sample $x"
  publishDir "pipeline_wyniki/pilon/${x}", mode: 'copy', pattern: 'latest_pilon.*'
  input:
  tuple val(x), path(bam1), path(fasta)
  output:
  // tuple val(x), path('latest_pilon.fasta'), path('latest_pilon.changes'), emit: ALL
  tuple val(x), path('latest_pilon.fasta'), emit: ONLY_GENOME
  script:
  """
  # indeksacja bam-ow
  /opt/docker/EToKi/externals/samtools index  $bam1

  # wywolanie pilon-a
  java -jar /opt/docker/EToKi/externals/pilon.jar  --threads ${params.cpus} --defaultqual 5 --genome $fasta --unpaired $bam1 --output pilon_polish --vcf --changes --outdir pilon_out

  # przygotowanie outputu
  # cat pilon_bwa/pilon_polish.fasta | awk  -v ALA=${x} 'BEGIN{OFS=""}; {if(\$0 ~ />/) print \$0,"_", ALA; else print \$0}' >> last_pilon.fasta
  cat pilon_out/pilon_polish.fasta >> latest_pilon.fasta
  cp  pilon_out/pilon_polish.changes latest_pilon.changes
  """
}

process run_medaka {
  // wygladzanie genomu medaka, nanopolish nie dziala bo nie mamy pliko fast5

  container  = 'salmonella_illumina:2.0'
  tag "Pilon for sample $x"
  publishDir "pipeline_wyniki/pilon/${x}", mode: 'copy', pattern: 'latest_pilon.*'
  input:
  tuple val(x), path(bam1), path(fasta)
  output:
  // tuple val(x), path('postmedaka.fasta'), path('medaka_annotated_filtered.vcf.gz'), emit: ALL
  tuple val(x), path('postmedaka.fasta'), emit: ONLY_GENOME
  script:
  """
  # indeksacja bam-ow
  samtools index $bam1
 
  MODEL="r941_min_hac_g507"  
  
  medaka consensus --model \${MODEL} \
                     --threads ${params.cpus} \
                     $bam1 \
                     forvariants.hdf

  
  medaka variant $fasta forvariants.hdf medaka.vcf
  medaka tools annotate medaka.vcf $fasta $bam1 medaka_annotated.vcf
  bcftools sort medaka_annotated.vcf >> medaka_annotated_sorted.vcf
  bgzip medaka_annotated_sorted.vcf
  tabix medaka_annotated_sorted.vcf.gz
  
  qual=13
  min_cov=20
  
  bcftools filter -O z -o medaka_annotated_filtered.vcf.gz -i "GQ > \${qual} && DP >= \${min_cov}" medaka_annotated_sorted.vcf.gz
  tabix medaka_annotated_filtered.vcf.gz


  cat $fasta | bcftools consensus medaka_annotated_filtered.vcf.gz >> postmedaka.fasta
  
  """
}

process run_kraken2_nanopore {
  // kopia processu do illuminy, ale z uwzglednieniem ze jest tylko jeden plik z odczytami
  tag "kraken2:${x}"
  container  = 'salmonella_illumina:2.0'
  publishDir "pipeline_wyniki/${x}/kraken2", mode: 'copy', pattern: "report_kraken2_individualreads.txt"
  publishDir "pipeline_wyniki/${x}/kraken2", mode: 'copy', pattern: "report_kraken2.txt"
  publishDir "pipeline_wyniki/${x}/", mode: 'copy', pattern: "summary_kraken.txt"
  containerOptions "--volume ${params.kraken2_db_absolute_path_on_host}:/home/external_databases/kraken2"
  maxForks 5

  input:
  tuple val(x), path(reads)

  output:
  tuple val(x), path('report_kraken2.txt'), path('report_kraken2_individualreads.txt'), path('summary_kraken.txt')

  script:
  """
  kraken2 --db /home/external_databases/kraken2 \
          --report report_kraken2.txt \
          --threads ${params.cpus} \
          --gzip-compressed \
          --minimum-base-quality ${params.quality_initial} \
          --use-names ${reads} >> report_kraken2_individualreads.txt 2>&1
  # parse kraken extract two most abundant species
  SPEC1=`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}' | grep -w S | sort -rnk 1 | head -1 | tr -s " " | cut -f6 | tr -d "="`
  SPEC2=`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}' | grep -w S | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f6 | tr -d "="`
  ILE1==`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}'| grep -w S | sort -rnk 1 | head -1 | tr -s " " | cut -f1 | tr -d " "`
  ILE2==`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}'| grep -w S | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f1 | tr -d " "`

  echo -e "${x}\t\${SPEC1}\${ILE1}%\t\${SPEC2}\${ILE2}%" >> summary_kraken.txt
  """
}



// SUB WORKFLOWS NANOPORE //

workflow pilon_first_nanopore {
// pilona puscimy tylko raz bo zakladam ze flye po cos te rundy filtrowania robi 
// zamieniono pilona na medaka
take:
initial_scaffold_inner
processed_fastq_inner_SE
main:
// laczymy kanaly ze scaffoldem genomu i odczytmai
for_remaping_SE_inner = initial_scaffold_inner.join(processed_fastq_inner_SE, by : 0, remainder : true)
// mapujemy odczyty na genom
single_bams_and_genome = run_minimap2(for_remaping_SE_inner)

// laczymy bam-a z mapowania z genomem

// Zamieniamy pilona na medaka
// final_assembly = run_pilon_nanopore(single_bams_and_genome)

final_assembly = run_medaka(single_bams_and_genome)

// liczymy pokrycia i fitrujemy slabe contig

for_remaping_polished_assembly = final_assembly.join(processed_fastq_inner_SE, by : 0, remainder : true)
single_bams_and_polished_genome = run_minimap2_2nd(for_remaping_polished_assembly)

final_assembly_filtered = run_coverage(single_bams_and_polished_genome)

emit:
// sub pipeline zwraca identyfikator probki + nowy scaffold
// bamy sa zbedne bo w kolejenj iteracji musza byc remapowane na poprawiony genom
final_assembly_filtered.ONLY_GENOME
}


// SUB WORKFLOWS ILLUMINA //

workflow pilon_first {
take:
initial_scaffold_inner
processed_fastq_inner_PE
processed_fastq_inner_SE

// initial_scaffold_inner jest kanalem zawierajacym
// ientyfikator probki
// plik fasta z genomem

// processed_fastq_inner_PE i SE kazdy jest kanalem zawierajacym
// PE -> identyfikator probki i odczyty PE
// SE -> identyfikator probki i odczyty SE
main:
// przygotowanie do mapowan
for_remaping_PE_inner = initial_scaffold_inner.join(processed_fastq_inner_PE, by : 0, remainder : true)
for_remaping_SE_inner = initial_scaffold_inner.join(processed_fastq_inner_SE, by : 0, remainder : true)

// mapowania
single_bams_inner = bwa_single(for_remaping_SE_inner)
paired_bams_inner = bwa_paired(for_remaping_PE_inner)

// laczenie wynikow mapowan
merged_bams_inner = paired_bams_inner.join(single_bams_inner, by : 0, remainder : true)
merged_bams_and_scaffold_inner = merged_bams_inner.join(initial_scaffold_inner,  by : 0, remainder : true)
run_pilon(merged_bams_and_scaffold_inner)

emit:
// sub pipeline zwraca identyfikator probki + nowy scaffold
// bamy sa zbedne bo w kolejenj iteracji musza byc remapowane na poprawiony genom 
run_pilon.out.ONLY_GENOME
}


workflow pilon_second {
// kopia pilon_first bez komentarzy
take:
initial_scaffold_inner
processed_fastq_inner_PE
processed_fastq_inner_SE
main:
for_remaping_PE_inner = initial_scaffold_inner.join(processed_fastq_inner_PE, by : 0, remainder : true)
for_remaping_SE_inner = initial_scaffold_inner.join(processed_fastq_inner_SE, by : 0, remainder : true)

single_bams_inner = bwa_single(for_remaping_SE_inner)
paired_bams_inner = bwa_paired(for_remaping_PE_inner)

merged_bams_inner = paired_bams_inner.join(single_bams_inner, by : 0, remainder : true)
merged_bams_and_scaffold_inner = merged_bams_inner.join(initial_scaffold_inner,  by : 0, remainder : true)
run_pilon(merged_bams_and_scaffold_inner)

emit:
run_pilon.out.ONLY_GENOME
}

workflow pilon_third {
// kopia pilon_first bez komentarzy
take:
initial_scaffold_inner
processed_fastq_inner_PE
processed_fastq_inner_SE
main:
for_remaping_PE_inner = initial_scaffold_inner.join(processed_fastq_inner_PE, by : 0, remainder : true)
for_remaping_SE_inner = initial_scaffold_inner.join(processed_fastq_inner_SE, by : 0, remainder : true)

single_bams_inner = bwa_single(for_remaping_SE_inner)
paired_bams_inner = bwa_paired(for_remaping_PE_inner)

merged_bams_inner = paired_bams_inner.join(single_bams_inner, by : 0, remainder : true)
merged_bams_and_scaffold_inner = merged_bams_inner.join(initial_scaffold_inner,  by : 0, remainder : true)
run_pilon(merged_bams_and_scaffold_inner)

emit:
run_pilon.out.ONLY_GENOME
}


workflow calculate_coverage {
take:
initial_scaffold_inner
processed_fastq_inner_PE
processed_fastq_inner_SE
main:

for_remaping_PE_inner = initial_scaffold_inner.join(processed_fastq_inner_PE, by : 0, remainder : true)
for_remaping_SE_inner = initial_scaffold_inner.join(processed_fastq_inner_SE, by : 0, remainder : true)

single_bams_inner = bwa_single(for_remaping_SE_inner)
paired_bams_inner = bwa_paired(for_remaping_PE_inner)

merged_bams_inner = paired_bams_inner.join(single_bams_inner, by : 0, remainder : true)

// W odroznieniu od workflow z piloenm laczymy pliki bam PE i SE w jeden plik
merged_bams_into_onefile_inner = merge_bams(merged_bams_inner)

merged_bams_into_onefile_and_scaffold_inner = merged_bams_into_onefile_inner.join(initial_scaffold_inner,  by : 0, remainder : true)
pokrycie = run_coverage(merged_bams_into_onefile_and_scaffold_inner)

emit:
// chyba .out jest mozliwa tylko gdy funkcja zwraca dwa emity, w jednym emit dostajemy sie
// do niego bez posrednictwa .out
pokrycie.ONLY_GENOME
}



// MAIN WORKFLOW //

workflow {

if(params.machine == 'Illumina') {
Channel
  .fromFilePairs(params.reads)
  .set {initial_fastq}

// FASTQC
run_fastqc(initial_fastq)

// kraken2 analiza kontmainacji innymi organizmami
run_kraken2_illumina(initial_fastq)

// Czyszczenie fastq jak w etoki
processed_fastq = clean_fastq(initial_fastq)

// Pierwsyz scaffold
initial_scaffold = spades(processed_fastq.All_path)

// Puszczamy 3-krotne wygladzanie pilon-em
first_polish_run = pilon_first(initial_scaffold, processed_fastq.PE_path, processed_fastq.SE_path)
second_polish_run = pilon_second(first_polish_run, processed_fastq.PE_path, processed_fastq.SE_path)
third_polish_run = pilon_third(second_polish_run, processed_fastq.PE_path, processed_fastq.SE_path)

// Liczymy pokrycie dla ostatecznie policzonych contigow, 
// usuwamy contigi o pokryciu ponizej 0.2 sredniego globalnego pokrycia zgodnie z kodem z etoki
// params.min_coverage_ratio = 0.2
// Analize robimy jednak moim kodem

final_assembly = calculate_coverage(third_polish_run, processed_fastq.PE_path, processed_fastq.SE_path)

} 

else {

//log.info "Nanopore"
//log.info params.reads
// Dla spojnosci z illumina niech kanal inital_fastq tez niech bedzie tuplem z pierwszym elementem jako identyfikatorem
Channel
  .fromPath(params.reads)
  .map {it -> tuple(it.getName().split("\\.")[0], it)}
  .set {initial_fastq}

// FASTQC
run_fastqc(initial_fastq)
// Aanliza zanieczyszczen
run_kraken2_nanopore(initial_fastq)

// Czyszczenie fastq
// Teoretycznie nie jest konieczne wedlug autorow flye
processed_fastq = clean_fastq_SE(initial_fastq)

// scaffold flye + 3 rundy wygladzania tym anrzedziem
initial_scaffold = run_flye(processed_fastq)

// jedna runda pilona poki co
final_assembly = pilon_first_nanopore(initial_scaffold, processed_fastq)

}

// Post-analizy
// Sprawdzono dzialaja zarowno dla illuminy jak i nnaopore wiec sa poza if-em
extract_final_stats(final_assembly)
MLST_out = run_7MLST(final_assembly)
parse_7MLST(MLST_out)
run_Seqsero(final_assembly)
run_sistr(final_assembly)
run_pointfinder(final_assembly)
run_plasmidfinder(final_assembly)
cgMLST_out = run_cgMLST(final_assembly)
parse_cgMLST_out = parse_cgMLST(cgMLST_out)
run_pHierCC(parse_cgMLST_out)
// extract_historical_data(parse_cgMLST_out) // NIE WYKONYWAC POKI CO POKI ENTEROBASE NIE POTWIERDZI ZE TAK MOZNA
prokka_out = run_prokka(final_assembly)
run_VFDB(prokka_out)
run_spifinder(final_assembly)

} 

