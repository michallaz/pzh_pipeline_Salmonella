params.cpus = 25
params.min_coverage_ratio = 0.1 
// Etoki odrzuca contigi ktorych srednie pokrycie jest mniejsze niz 0.2 globalnego pokrycia
// my nie jestesmy tak ostrzy dajemy 0.1 co tlumaczymy uzyciem innego alignera
// params.min_coverage_ratio_nanopore = 0.1 // Ten parametr nie jest moze  posluzy do nanopre'a

params.species = 's.enterica' // inne mozliwosci to c.jejuni i e.coli
params.reads = '/mnt/sda1/michall/Salmonella/*_R{1,2}_001.fastq.gz'

// eqluowane w main worflow jesli jest cos innego niz illumina
// to wchozi sciezka nanoporowe'a

params.machine = 'Illumina'

params.quality_initial = 5 // Parametr stosowany aktualnie tylko przez krakena


// bazy specyficzne dla organizmu
if ( params.species  == 's.enterica' ) {
	params.Achtman7GeneMLST_db_absolute_path_on_host = "/mnt/sda1/michall/db/Salmonella/MLST_Achtman"  
	params.cgMLST_db_absolute_path_on_host = "/mnt/sda1/michall/db/Salmonella/cgMLST_v2" // sciezka do alleli z cgMLST
	params.phiercc_db_absolute_path_on_host = "/mnt/sda1/michall/db/Salmonella/pHierCC_local" // wyniki wlasnego klastrowania
	params.Enterobase_db_absolute_path_on_host = "/mnt/sda1/michall/db/Salmonella/Enterobase"
} else if ( params.species  == 'e.coli' ) {
	params.Achtman7GeneMLST_db_absolute_path_on_host = "/mnt/sda1/michall/db/Ecoli/MLST_Achtman" //ponownie na sztywno do poprawy docelowo
        params.cgMLST_db_absolute_path_on_host = "/mnt/sda1/michall/db/Ecoli/cgMLST_v1" // sciezka do alleli z cgMLST
        params.phiercc_db_absolute_path_on_host = "/mnt/sda1/michall/db/Ecoli/pHierCC_local"
        params.Enterobase_db_absolute_path_on_host = "/mnt/sda1/michall/db/Ecoli/Enterobase"
} else if ( params.species  == 'c.jejuni' ) {
	println("This option is not yet implemented")
} else {
	println("Incorrect species provided")
	System.exit(0)
}

// TOKEN
params.enterobase_api_token = "eyJhbGciOiJIUzI1NiIsImlhdCI6MTcyMDQzNjQxMSwiZXhwIjoxNzM2MjA0NDExfQ.eyJfIjoibUsyNFlZSHd4SyIsInVzZXJuYW1lIjoiTWljaGFsX0xhem5pZXdza2kiLCJpZCI6ODg4MCwiYWRtaW5pc3RyYXRvciI6bnVsbCwiZW1haWwiOiJtbGF6bmlld3NraUBwemguZ292LnBsIiwiYXBpX2FjY2Vzc19jbG9zdHJpZGl1bSI6IlRydWUiLCJhcGlfYWNjZXNzX2Vjb2xpIjoiVHJ1ZSIsImFwaV9hY2Nlc3Nfc2VudGVyaWNhIjoiVHJ1ZSJ9.VEsyVPv8sn1zG7d3uFqEjfk6XFS2qP8P5Y5mh9VPE9w" // klucz api nadawany przez ENTEROBASE po rejestracji na stronie + wystapieniu o klucz 

// bazy ogolne
params.AMRFINDER_db_absolute_path_on_host = "/mnt/sda1/michall/db/AMRfider_plus"
params.metaphlan_db_absolute_path_on_host = "/mnt/sda1/michall/db/Metaphlan"
params.kmerfinder_db_absolute_path_on_host = "/mnt/sda1/michall/db/Kmerfinder/kmerfinder_db"
params.kraken2_db_absolute_path_on_host = "/home/michall/kraken2/kraken2_db/kraken2_sdb/" 

// Kontenery uzywane w tym skrypcie 
// salmonella_illumina:2.0 - bazowy kontener z programami o kodem
// staphb/prokka:latest - kontener z prokka, w notatkach mam ze budowanie programu od 0 jest meczace bo to kod sprzed ponad 4 lat 



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
  maxForks 5
  input:
  tuple val(x), path(reads)
  output:
  tuple path("*forward_fastqc.txt"), path("*reverse_fastqc.txt"), path("*fastqc.html")

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

process run_fastqc_nanopore {
  // Lekko zmodyfkowany modul z sars-a
  tag "fastqc for sample ${x}"
  container  = 'salmonella_illumina:2.0'
  publishDir "pipeline_wyniki/${x}/QC", mode: 'copy'
  maxForks 5
  input:
  tuple val(x), path(reads)
  output:
  tuple path("*_fastqc.txt"), path("*fastqc.html")

  script:
  """
  QUALITY=2
  fastqc --format fastq \
         --threads ${params.cpus} \
         --memory 2024 \
         --extract \
         --delete \
         --outdir . \
         ${reads}
    r1=\$(basename ${reads} .fastq.gz)
    /data/parse_fastqc_output.py \${r1}_fastqc/fastqc_data.txt \${QUALITY} >> ${x}_fastqc.txt
  """
}

process clean_fastq {
  // Prosta funkcja zaimplementowa w etoki do czyszczenia plikow fastq, ma na celu rozdzielenie odczytow ktore sa sparowane
  // od tych ktore pary nie maja, trimmowanie odczytow na podstawie jakosci, zmiana nazwy odczytow .
  // Ta funkcja generalnie mimikuje dzialanie trimmomatic-a. Wiec teoretycznie mozna uzyc jego
  // Default na trimming to quality 6

  container  = 'salmonella_illumina:2.0'
  tag "Fixing fastq dla sample $x"
  maxForks 5
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
    # w niektorych przypadkach plik SE moze nic nie zawierac wtedy Etoki sie gubi i ie generuje poprawnie pliku SE
  # tworzymy samodzilenie plik z dummy sekwencja tak by reszta skryptu dzialala
  if [ -e prep_out.1.0.s.fastq.gz ]; then
      echo "@0_SE_0" >> prep_out_L1_SE.fastq
      echo "GTACTGACCAAACTGGTCGTGTAGCGTTTCATGCCACATCGTATTTTCGGCCATTGGCTGATACCTCCATTGTTAACACCCGTAAAAAAAGGGCGCAACATCATAGCTAACAATGACCGTGGATGCACGGTCATTATTTCAGCAATAGGAT" >> prep_out_L1_SE.fastq
      echo "+" >> prep_out_L1_SE.fastq
      echo "AFFFFFFFFFFFFAFFFFFFFFAFFFAFFF/FAFFAAAFF/=FFAAFFAFFFF/FFFFFF//FFFFFFFFA/FFFFFAFFAA=FFFFFFFF/FFFFFFFFF/FFF/FFFFFFFFAFFFAFFFFFAFAF/FF=FFFFAFFFFF/FF/FAFAF" >> prep_out_L1_SE.fastq
      gzip prep_out_L1_SE.fastq
  fi
  """
}


process spades {
  // Funkcja do odpalania spadesa
  // Powstaly plik fasta jest poprawiany bo 
  // Hapo-G nie akceptuje "." w plikach fasta
  cpus params.cpus
  container  = 'salmonella_illumina:2.0'
  tag "Spades dla sample $x"
  maxForks 5
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
  maxForks 5
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
  maxForks 5
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
  maxForks 5
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
  // Dokladne uzycie pilona jest tu https://github.com/broadinstitute/pilon/wiki/Requirements-&-Usage
  container  = 'salmonella_illumina:2.0'
  tag "Pilon for sample $x"
  maxForks 5
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
  maxForks 5
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
  /opt/docker/EToKi/EToKi.py MLSTdb -i /Achtman7GeneMLST_entero/all_allels.fasta -x 0.8 -m 0.5 -r /Achtman7GeneMLST_entero/MLST_Achtman_ref.fasta -d MLST_database.tab
  /opt/docker/EToKi/EToKi.py MLSType -i $fasta -r /Achtman7GeneMLST_entero/MLST_Achtman_ref.fasta -k ${x} -o MLSTout.txt -d MLST_database.tab
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

# zostawiamy wczytanie profili dla 7 genomwego bo jest szybkie i pozwala ladnie parsowac wynik
# nie jest jednak potrzebne dla samego szukania pasujacego ST dla zidentyfikowanego zestawu alleli w probce 

known_profiles, klucze = create_profile('/Achtman7GeneMLST_entero/profiles.list')
identified_profile = parse_MLST_fasta('MLSTout.txt')
matching_profile, min_value, _, _ = getST(MLSTout = identified_profile, \
                                       profile_file = '/Achtman7GeneMLST_entero/profiles.list')

  
if "${params.species}" == 's.enterica':
    klucze_sorted = ['aroC', 'dnaN', 'hemD', 'hisD', 'purE', 'sucA', 'thrA']
elif "${params.species}" == 'e.coli':
    klucze_sorted = ['adk', 'fumC', 'gyrB', 'icd', 'mdh', 'purA', 'recA']

elif "${params.species}" == 'campylobacter':
    pass

header="ST\\t"+"\\t".join(klucze_sorted) + "\\tDistance\\n"
formatted_string = ""
for klucz in klucze_sorted:
    try:
        formatted_string += f"{identified_profile[klucz]}\\t"
    except KeyError:
        formatted_string += f"-1\\t" 

with open('parsed_7MLST.txt', 'w') as f:
    f.write(header)
    f.write(f'{matching_profile}\\t{formatted_string}{min_value}\\n')
		

"""
}

process run_Seqsero {
  // SeqSero only works for Salmonella
  container  = 'salmonella_illumina:2.0'
  tag "Predicting OH for sample $x with Seqsero"
  publishDir "pipeline_wyniki/${x}", mode: 'copy'
  cpus params.cpus
  maxForks 5
  input:
  tuple val(x), path(fasta)
  output:
  tuple val(x), path('seqsero_out/SeqSero_result.txt')
  script:
  """
  # -m to rodzaj algorytmu -m to chyba opart o k-mery
  # -t 4 to informacja ze inputem sa contigi z genomem
  # -p to procki 
  if [ "${params.species}" == 's.enterica' ];then 
      python /opt/docker/SeqSero2/bin/SeqSero2_package.py -m k -t 4 -p ${task.cpus} -i $fasta -d seqsero_out
  else
      # Prepare empty output dor pipeline to work
      mkdir seqsero_out
      touch seqsero_out/SeqSero_result.txt
  fi
  """
}

process run_sistr {
  container  = 'salmonella_illumina:2.0'
  tag "Predicting OH for sample $x with Sistr"
  publishDir "pipeline_wyniki/${x}", mode: 'copy'
  cpus params.cpus
  maxForks 5
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
  if [ "${params.species}" == 's.enterica' ];then
    /usr/local/bin/sistr --qc -vv --alleles-output allele-results.json --novel-alleles novel-alleles.fasta --cgmlst-profiles cgmlst-profiles.csv -f tab -t ${task.cpus} -o sistr-output.tab $fasta
  else
    ## https://github.com/phac-nml/ecoli_serotyping.git
    touch sistr-output.tab
  fi

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
  # params.species jest dobrany tak by rozumial go resfinder
  # ale dla campylo sa 2 wersje c.jejuni i c.coli 
  
  python -m resfinder -o resfinder_out/ -s ${params.species}  -l 0.6 -t 0.8 --acquired --point -k /opt/docker/kma/kma -db_disinf /opt/docker/disinfinder_db/ -db_res /opt/docker/resfinder_db/ -db_point /opt/docker/pointfinder_db/ -ifa ${fasta}
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
  maxForks 5
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
  /data/run_blastn_ver11.sh $fasta ${task.cpus}
  cat log.log | cut -f1,2 > cgMLST.txt
  """
}

process parse_cgMLST {
  // Funkcja do parsowania wynikow run_cgMLST
  // Zwraca plik parsed_cgMLST.txt ktory zawiera w pierwszej kolumnie identyfikator ST (albo pojedyncza cyfra jesli dany ST jest juz znany bazie enterobase) albo 
  // local_cyfra (jesli ST jest nieznany bazie enterobase)
  // druga kolumna w tym pliku to albo 0 (dany ST zostal znaleziony w bazie) albo NOVEL jesli w bazie entero ani local nie ma takiego sample'a 
  container  = 'salmonella_illumina:2.0'
  containerOptions "--volume ${params.cgMLST_db_absolute_path_on_host}:/cgMLST2_entero"
  tag "Parsing cgMLST for sample $x"
  publishDir "pipeline_wyniki/${x}", mode: 'copy', pattern: 'parsed_cgMLST.txt'
  publishDir "pipeline_wyniki/${x}", mode: 'copy', pattern: 'matching_allel_list.txt'
  maxForks 1 // ustawiamy maxforks 1 dzieki temu mam nadzieje mamy pewnosc ze baza "local" bedzie poprawnie updatowana
  // choc na pewno zadziala to wolniej
  input:
  tuple val(x), path('cgMLST.txt')
  output:
  tuple val(x), path('parsed_cgMLST.txt'), path('matching_allel_list.txt')
  script:
"""
#!/usr/bin/python
import  sys
sys.path.append('/data')
from all_functions_salmonella import *

#known_profiles, klucze = create_profile('/cgMLST2_entero/profiles.list')
identified_profile = parse_MLST_blastn('cgMLST.txt')
matching_profile, min_value, allele_list_lowest_difference, sample_profile = getST(MLSTout = identified_profile, 
                                                                   profile_file = '/cgMLST2_entero/profiles.list')


if min_value == 0:
    # found a matching profile in "core" database
    with open('parsed_cgMLST.txt', 'w') as f:
        f.write(f'{matching_profile}\\t{min_value}\\tIN_ENTERO\\n')
else:
    # look for a profile in "local" database 
    matching_profile_local, min_value_local, allele_list_lowest_difference_local, _ = getST(MLSTout = identified_profile,
                                                                                   profile_file = '/cgMLST2_entero/local/profiles_local.list')
    if min_value_local == 0:
        # Znalazlem wpis w bazie local (nie musze updatowac lokalnych plikow)
        # ale zapisuje odleglosc zarowno do ST z bazy entero (potrzebuje do phierCC)
        with open('parsed_cgMLST.txt', 'w') as f:
            f.write(f'{matching_profile_local}\\t{min_value}\\tIN_LOCAL\\n')
    else:
        # nie znaleziono pasujacego profilu rowniez w bazie local
        # dodaje nowy wpis do tej bazy
        # Musze pobrac wpis local ST maja postac local_1, local_2 itd ...
        last_ST = 0
        with open('/cgMLST2_entero/local/profiles_local.list') as f:
            for line in f:
                line = line.rsplit()
                if 'local' in line[0]:
                    last_ST = line[0].split('_')[1]
        
        # dodajemy nowy ST do bazy local
        
        novel_profile_ST = f'local_{int(last_ST)+1}'
        to_save = "\t".join(map(str, sample_profile))
        write_novel_sample(f'{novel_profile_ST}\\t{to_save}\\n', '/cgMLST2_entero/local/profiles_local.list')
        
        # no na koniec zapisuje wynik do outputu, ALE UWAGA MIN_VALUE TO ZAWSZE ODLEGLOSC DO NAJBLIZSZEGO ST Z ENTERO !

        with open('parsed_cgMLST.txt', 'w') as f:
            f.write(f'{novel_profile_ST}\\t{min_value}\\tNOVEL_LOCAL\\n')

# ZAPISUJE NAJBLIZSZY ST I JEGO ALLELE Z BAZY ENETRO  potrzebujeto do liczenia hierCC
with open('matching_allel_list.txt', 'w') as f:
    f.write(allele_list_lowest_difference)
         


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
  maxForks 5
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
  maxForks 5
  input:
  // inputem jest output procesu run_prokka
  tuple val(x), path(gff), path(faa), path(ffn), path(tsv)
  output:
  tuple val(x), path('VFDB_summary.txt')
  script:
  """
  if [ ${params.species} == 's.enterica' ]; then
      SPEC="Salmonella" # wpradzie jest juz paramter params.species, ale baza VFDB go nie zrozumie.
  elif [ ${params.species} == 'e.coli' ]; then
      SPEC="Escherichia"
  elif [ ${params.species} == 'c.jejuni' ]; then
      SPEC="Campylobacter"
  fi 
  PIDENT=80 # minimalna identycznosc sekwencyjna aby stwierdzic ze jest hit 
  COV=80 # minimalne pokrycie query i hitu aby stwierdzic ze jest hit 
  EVAL=0.01  # maksymalne e-value
  /opt/docker/EToKi/externals/run_VFDB.sh $ffn ${task.cpus} \${SPEC} \${PIDENT} \${EVAL} \${COV}
  """

}

//process parse_VFDB {

//  Tutaj bedziemy parsowac output VFDB w celu okreslenia dla ECOLI co to za podtyp

//}

process run_spifinder {
  container  = 'salmonella_illumina:2.0'
  tag "Predicting virulence islands for sample $x"
  publishDir "pipeline_wyniki/${x}", mode: 'copy'
  cpus params.cpus
  maxForks 5
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
  if [ ${params.species} == 's.enterica' ]; then
      python /opt/docker/spifinder/spifinder.py -i $fasta -o spifinder_results -mp blastn -p /opt/docker/spifinder_db/ -l 0.6 -t 0.9 -x
  else
      # spifider jest tylko dla Salmonelli
      touch spifinder_results/dummy.txt
  fi
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

process run_metaphlan_illumina {
  // modul wziety z pipeline do SARS
  // poprawilem tylko pare rzeczy zwiazane z parametrami i kontenerami ktore u mnie sa inne
  // kraken2 instalowany jest przez ETOKI i jest w path kontenera do salmonelli
  tag "kraken2:${x}"
  container  = 'salmonella_illumina:2.0'
  publishDir "pipeline_wyniki/${x}/metaphlan", mode: 'copy', pattern: "report_metaphlan*"
  containerOptions "--volume ${params.metaphlan_db_absolute_path_on_host}:/bowtie_db"
  maxForks 5
  cpus params.cpus
  input:
  tuple val(x), path(reads)

  output:
  tuple val(x), path('report_metaphlan_SGB.txt'), path('report_metaphlan_species.txt'), path('report_metaphlan_genra.txt')

  script:
  """
  metaphlan ${reads[0]},${reads[1]} --bowtie2out metagenome.bowtie2.bz2 --nproc ${task.cpus} --input_type fastq -o profiled_metagenome.txt --bowtie2db /bowtie_db/ --unclassified_estimation
  # Parsujemy wyniki
  metaphlan metagenome.bowtie2.bz2 --input_type bowtie2out --bowtie2db /bowtie_db/  --nproc ${task.cpus} --tax_lev 't' -o report_metaphlan_SGB.txt
  metaphlan metagenome.bowtie2.bz2 --input_type bowtie2out --bowtie2db /bowtie_db/  --nproc ${task.cpus} --tax_lev 's' -o report_metaphlan_species.txt
  metaphlan metagenome.bowtie2.bz2 --input_type bowtie2out --bowtie2db /bowtie_db/  --nproc ${task.cpus} --tax_lev 'g' -o report_metaphlan_genra.txt
  """
}

process run_kmerfinder {
  tag "kmerfinder:${x}"
  container  = 'salmonella_illumina:2.0'
  publishDir "pipeline_wyniki/${x}/kmerfinder", mode: 'copy', pattern: "results*"
  containerOptions "--volume ${params.kmerfinder_db_absolute_path_on_host}:/kmerfinder_db"
  maxForks 5
  cpus params.cpus
  input:
  tuple val(x), path(reads)

  output:
  tuple val(x),path('results.spa'), path('results.txt')

  script:
  """
  /opt/docker/kmerfinder/kmerfinder.py -i ${reads[0]} ${reads[1]} -o ./kmerfider_out -db /kmerfinder_db/bacteria/bacteria.ATG -tax /kmerfinder_db/bacteria/bacteria.tax -x -kp /opt/docker/kma/
  cp kmerfider_out/results.spa .
  cp kmerfider_out/results.txt .
  """
}

process run_pHierCC {
  
  // Funkcja odpytuje API Enterobase w celu wyciagniecia profuili z tej bazu
  // Ponadto Funkcja odpytuje dwa pliki przygotowane przeze mnie 
  // Jeden zawierajacy klastrowanie SINGLE linkage zbudowane na 430k profili z enterobase
  // Drugi zawierajacy klastrowanie COMPLETE linkage zbudowane na 430k profili z enterobae
  // Jesli ST jest nowy NIE obecny w zadnej z baz program nic nie zwraca ?
  maxRetries 3
  errorStrategy 'retry' 
  container  = 'salmonella_illumina:2.0'
  // containerOptions "--volume ${params.cgMLST_db_absolute_path_on_host}:/cgMLST2_entero"
  containerOptions "--volume ${params.phiercc_db_absolute_path_on_host}:/pHierCC_local"
  tag "Predicting hierCC levels for sample $x"
  publishDir "pipeline_wyniki/${x}/phiercc", mode: 'copy'
  input:
  tuple val(x), path('parsed_cgMLST.txt'), path('matching_allel_list.txt')
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
import re
import numpy as np
import time

# rzadkich przypadkach gdyby bylo na raz za duzo requestow
# Enterobase zwraca blad wiec zrandomizujmy moment wejscia
time.sleep(np.random.randint(2,20))

API_TOKEN = "${params.enterobase_api_token}"

def __create_request(request_str):
    base64string = base64.b64encode('{0}: '.format(API_TOKEN).encode('utf-8'))
    headers = {"Authorization": "Basic {0}".format(base64string.decode())}
    request = urllib.request.Request(request_str, None, headers)
    return request

def getST(my_file):
    # prosta funkcja zwraca numer ST sample'a oraz odleglosc do najblzszego znanego ST z bazy enterobase
    with open(my_file) as f:
        for line in f:
            line = line.rsplit()
    return line[0], int(line[1])

def get_matching_ST(my_file):
    # funkcja bierze plik matching_allel_list.txt' i wyciaga z niego pierwszy element (najblizy ST z bazye entero)
    # i elementy 1: (profil alleli tego ST)
    with open(my_file) as f:
        for line in f:
            elementy = line.rsplit()
    return elementy[0], elementy[1:]

if "${params.species}" == 's.enterica':
    DATABASE="senterica"
    scheme_name="cgMLST_v2"
    phiercc_header='ST\\tHC0\\tHC2\\tHC5\\tHC10\\tHC20\\tHC50\\tHC100\\tHC200\\tHC400\\tHC900(ceBG)\\tHC2000\\tHC2600\\tHC2850(subsp.)\\n'
    lista_kluczy = ['d0', 'd2', 'd5', 'd10', 'd20', 'd50', 'd100', 'd200' , 'd400', 'd900', 'd2000', 'd2600', 'd2850']
elif "${params.species}" == 'e.coli':
    DATABASE="ecoli"
    scheme_name="cgMLST"
    phiercc_header='ST\\tHC0\\tHC2\\tHC5\\tHC10\\tHC20\\tHC50\\tHC100\\tHC200\\tHC400\\tHC1100(cgST Cplx)\\tHC1500\\tHC2000\\tHC2350(subsp.)\\n'
    lista_kluczy = ['d0', 'd2', 'd5', 'd10', 'd20', 'd50', 'd100', 'd200' , 'd400', 'd1100', 'd1500', 'd2000', 'd2350']
elif "${params.species}" == 'c.jejuni':
    # No campylobacter in Enterobase 
    # need to write separate module for that bacteria
    pass


my_ST, my_dist = getST('parsed_cgMLST.txt')
matching_ST, matching_ST_allels =  get_matching_ST('matching_allel_list.txt')


# mamy nastepujace mozliwosci  my_dist wynosi 0 lub nie-0  (czyli znalzlem cos niwego badz nie)
# jesli wynosi 0 nie ma 
# jesli jest cos nowego do ide sciezka (cos starego i z bazy entero) ale modyfikuje linijki w zaleznosci od dystansu i dodaje je do bazy local


# Tworzenie lokalnych plikow wynikowych
with open('parsed_phiercc_enterobase.txt', 'w') as f:
    f.write(phiercc_header)

with open('parsed_phiercc_minimum_spanning_tree.txt', 'w') as f:
    f.write(phiercc_header)

with open('parsed_phiercc_maximum_spanning_tree.txt', 'w') as f:
    f.write(phiercc_header)

if "${params.species}" == 'c.jejuni':
    # Output files exist, module will not crash
    exit(1)
 
# 1. Szukanie w bazie enterobase
with open('parsed_phiercc_enterobase.txt', 'a') as f:
    address = f"https://enterobase.warwick.ac.uk/api/v2.0/{DATABASE}/{scheme_name}/sts?st_id={matching_ST}&scheme={scheme_name}&limit=5"
    #try:
    response = urlopen(__create_request(address))
    data = json.load(response)
    
    lista_poziomow = [data['STs'][0]['info']['hierCC'][x] for x in lista_kluczy] # lista z uporzadkowanymi poziomami
            
    # modyfikacja wartosci w lista_poziomow na podstawie 
    # wartsoci z hierCC z odelgoscia mniejsza niz dystans do najblzszego znangeo przedstawiciela z enterobase sa podmieniane na my_ST 
    try:
        last_index = np.where(list(map(lambda x: int(re.findall('\\d+', x)[0]) < my_dist, lista_kluczy)))[0][-1]
        lista_poziomow[:(last_index+1)] = [my_ST] * (last_index + 1)
    except IndexError:
        pass
    formatted_string = "\t".join(list(map(str, lista_poziomow)))
    f.write(f'{my_ST}\\t{formatted_string}\\n')
    
    # Wywalm exccepta na rzecz maxRetries 
    #except HTTPError as Response_error:
    #    print(f"{Response_error.code} {Response_error.reason}. URL: {Response_error.geturl()}\\n Reason: {Response_error.read()}")

# 2. Szukanie w wynikach mojego klastrowania z uzyciem single linkage
with open('parsed_phiercc_minimum_spanning_tree.txt', 'a') as f, gzip.open('/pHierCC_local/profile_single_linkage.HierCC.gz') as f2, open('/pHierCC_local/profile_single_linkage.HierCC.index') as f3:

    # wstepne szukanie w f3, wynikow szukamy w f2, a zapisujemy do f
    pointer = 0
    for line in f3:
        line = line.rsplit()
        try:
            if int(matching_ST) < int(line[0]):
                break
            else:
                pointer = int(line[1])

        except ValueError:
            pass
            # pierwszy wiersz w indekszie to naglowek
    # ustaw kursow blizej lokalziacji przed szukanym ST
    f2.seek(pointer)
	
    for line in f2:
        line = list(map(lambda x: x.decode('utf-8', errors='replace'), line.split()))
        if line[0] == matching_ST:
            # Znalzlem linijke z najblizszym ST poprawiam ja aby uwzglednic nie idealny hit
            
            if "${params.species}" == 's.enterica':
                lista_poziomow = [line[1], line[3], line[6], line[11], line[21], line[51], line[101],  line[201], line[401], line[901], line[2001], line[2601], line[2851]]
            elif "${params.species}" == 'e.coli':
                lista_poziomow = [line[1], line[3], line[6], line[11], line[21], line[51], line[101],  line[201], line[401], line[1101], line[1501], line[2001], line[2351]]
            try:
                last_index = np.where(list(map(lambda x: int(re.findall('\\d+', x)[0]) < my_dist, lista_kluczy)))[0][-1]
                lista_poziomow[:(last_index + 1)] = [my_ST] * (last_index +1)
            except IndexError:
                pass
            formatted_string = "\t".join(list(map(str, lista_poziomow)))
            f.write(f'{my_ST}\\t{formatted_string}\\n')
            # nie ma potrzeby dalszego ogladania pliku
            break

# 3. Szukanie pHierCC w wynikach 
with open('parsed_phiercc_maximum_spanning_tree.txt', 'a') as f, gzip.open('/pHierCC_local/profile_complete_linkage.HierCC.gz') as f2, open('/pHierCC_local/profile_complete_linkage.HierCC.index') as f3:

    # wstepne szukanie w f3, wynikow szukamy w f2, a zapisujemy do f
    pointer = 0
    for line in f3:
        line = line.rsplit()
        try:
            if int(matching_ST) < int(line[0]):
                break
            else:
                pointer = int(line[1])

        except ValueError:
            pass
            # pierwszy wiersz w indekszie to naglowek
    # ustaw kursow blizej lokalziacji przed szukanym ST
    f2.seek(pointer)

    for line in f2:
        line = list(map(lambda x: x.decode('utf-8', errors='replace'), line.split()))
        if line[0] == matching_ST:
            # Znalzlem linijke z najblizszym ST poprawiam ja aby uwzglednic nie idealny hiti
            if "${params.species}" == 's.enterica':
                lista_poziomow = [line[1], line[3], line[6], line[11], line[21], line[51], line[101],  line[201], line[401], line[901], line[2001], line[2601], line[2851]]
            elif "${params.species}" == 'e.coli':
                lista_poziomow = [line[1], line[3], line[6], line[11], line[21], line[51], line[101],  line[201], line[401], line[1101], line[1501], line[2001], line[2351]]
            try:
                last_index = np.where(list(map(lambda x: int(re.findall('\\d+', x)[0]) < my_dist, lista_kluczy)))[0][-1]
                lista_poziomow[:(last_index+1)] = [my_ST] * (last_index+1)
            except:
                formatted_string = "\t".join(list(map(str, lista_poziomow)))
                f.write(f'{my_ST}\\t{formatted_string}\\n')
            # nie ma potrzeby dalszego ogladania pliku
            break            
    
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
  containerOptions "--volume ${params.Enterobase_db_absolute_path_on_host}:/Enterobase"
  publishDir "pipeline_wyniki/${x}/", mode: 'copy'
  input:
  tuple val(x), path('parsed_phiercc_enterobase.txt'), path('parsed_phiercc_minimum_spanning_tree.txt'), path('parsed_phiercc_maximum_spanning_tree.txt')
  output:
  tuple val(x), path('enterobase_historical_data.txt')
  script:
"""
#!/usr/bin/python
import numpy as np
def get_hiercc_level(my_file):
    with open(my_file) as f:
        i = 0
        for line in f:
            if i == 0:
                klucze = line.rsplit()
            elif i == 1:
                wartosci = line.rsplit()
            i+=1
    slownik = {x:y for x,y in zip(klucze,wartosci)}
    return slownik
# Tutaj mamy informacje o poziomach mojej probki
slownik_hiercc = get_hiercc_level('parsed_phiercc_enterobase.txt')

# Szukamy strainow o takim wspolnym poziomie, to powinno byc user defined
phiercc_level_userdefined = '5'


common_STs = slownik_hiercc[f'HC{phiercc_level_userdefined}']

# teraz lecimy po tablicy STs i szukamy ST z pozioem jak common_STs na poziomie phiercc_level_userdefined

expected_ST = {}
STs = np.load('/Enterobase/sts_table.npy', allow_pickle=True)
STs = STs.item()

for klucz,wartosc in STs.items(): 
    if wartosc[f'd{phiercc_level_userdefined}'] == common_STs:
        expected_ST[klucz] = ''
    		

# na tym etapie znamy STs ktore na poziomie phiercc_level_userdefined maja wartosc common_STs
# te STs sa znane w tablicy straindata

straindata = np.load('/Enterobase/straindata_table.npy',  allow_pickle=True)
straindata = straindata.item()

# niestety lecimy po wszystkich elementach w straindata
dane_historyczne = {} # kluczem jest nazwa szczepu , wartoscia 2 elementowa lista z krajem i rokiem
for strain, wartosc in straindata.items():
    for scheme in wartosc['sts']:
        if 'cgMLST_v2' in scheme.values() and str(scheme['st_id']) in expected_ST.keys():
            dane_historyczne[strain] = [wartosc['country'], wartosc['collection_year'], scheme['st_id']]  

# zapisujemy dane
with open('enterobase_historical_data.txt', 'w') as f:
    f.write(f'Strain_id\\tCountry\\tYear\\tST\\n')
    for klucz, wartosc in dane_historyczne.items(): 
        f.write(f'{klucz}\\t{wartosc[0]}\\t{wartosc[1]}\\t{wartosc[2]}\\n')
"""
}

process run_amrfinder {
  container  = 'salmonella_illumina:2.0'
  tag "Predicting microbial resistance with AMRfinder for sample $x"
  publishDir "pipeline_wyniki/${x}/AMRplus_fider", mode: 'copy', pattern: "AMRfinder*"
  containerOptions "--volume ${params.AMRFINDER_db_absolute_path_on_host}:/AMRfider"
  input:
  tuple val(x), path(fasta)
  output:
  tuple val(x), path('AMRfinder_resistance.txt'), path('AMRfinder_virulence.txt')
  // -n input, plik z sekwencja nukleotydowa
  // -d input, sciezka do bazy AMRfindera 
  // -i input, seq identity miedzy targetem a query
  // -c input, query coverage z blasta
  // -O input, nazwa organizmu 
  // -o outpu, nazwa pliku z outputem
  // --plus input, Add the plus genes to the report
  // The 'plus' subset include a less-selective set of genes of interest including genes involved in virulence, biocide, heat, metal, and acid resistance
  // --blast_bin input sciezka do binarek blast-a, podaje wxplicite po w kontenerze sa 2 binarki blasta te z etoki i instalowane recznie
  // Te z etoki sa za stare 
  script:
  """
  if [ ${params.species} == 's.enterica' ]; then
      SPEC="Salmonella" # wpradzie jest juz paramter params.species, ale baza VFDB go nie zrozumie.
  elif [ ${params.species} == 'e.coli' ]; then
      SPEC="Escherichia"
  elif [ ${params.species} == 'c.jejuni' ]; then
      SPEC="Campylobacter"
  fi
  amrfinder --blast_bin /blast/bin -n $fasta -d /AMRfider  -i 0.9 -c 0.5 -o initial_output.txt -O \${SPEC} --plus
  cat initial_output.txt | grep -w AMR >> AMRfinder_resistance.txt
  cat initial_output.txt | grep -w VIRULENCE >> AMRfinder_virulence.txt

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
  // --deterministic -  perform disjointig assembly single-threaded program zwraca dla tych samych danych rozne wyniki
  // wedlug issue z githuba https://github.com/mikolmogorov/Flye/issues/640
  // ten problem dalej istnieje 

  container  = 'salmonella_illumina:2.0'
  tag "Predicting scaffold with flye for sample $x"
  publishDir "pipeline_wyniki/${x}", mode: 'copy', pattern: "*"
  cpus params.cpus
  maxForks 5
  input:
  tuple val(x), path(fastq_gz)
  output:
  tuple val(x), path('output/assembly.fasta') 
  
  script:
  """
  # /data/Flye to sciezka z Flye instalowanego z github, uwaga
  # w kontenerze tez jest flye instalowant przez etoki i ten jest w PATH
  /data/Flye/bin/flye --nano-raw ${fastq_gz} -g 6m -o output -t ${params.cpus} -i 3 --no-alt-contig --deterministic

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
  maxForks 5
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
  maxForks 5
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
  publishDir "pipeline_wyniki/${x}/pilon", mode: 'copy', pattern: 'latest_pilon.*'
  maxForks 5
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
  publishDir "pipeline_wyniki/${x}/medaka", mode: 'copy', pattern: 'latest_pilon.*'
  maxForks 5
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
  LEVEL="G" # G to chyba rodzzaj (genus ?)  S to pewnie gatunek (SPECIES)

  SPEC1=`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -1 | tr -s " " | cut -f6 | tr -d "="`
  SPEC2=`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f6 | tr -d "="`
  ILE1==`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}'| grep -w \${LEVEL} | sort -rnk 1 | head -1 | tr -s " " | cut -f1 | tr -d " "`
  ILE2==`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}'| grep -w \${LEVEL} | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f1 | tr -d " "`

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

// Metaphlan
run_metaphlan_illumina(initial_fastq)

//Kmerfinder
run_kmerfinder(initial_fastq)

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

else if (params.machine == 'Nanopore') {

//log.info "Nanopore"
//log.info params.reads
// Dla spojnosci z illumina niech kanal inital_fastq tez niech bedzie tuplem z pierwszym elementem jako identyfikatorem
Channel
  .fromPath(params.reads)
  .map {it -> tuple(it.getName().split("\\.")[0], it)}
  .set {initial_fastq}

// FASTQC
run_fastqc_nanopore(initial_fastq)
// Aanliza zanieczyszczen
run_kraken2_nanopore(initial_fastq)

// Czyszczenie fastq
// Teoretycznie nie jest konieczne wedlug autorow flye
processed_fastq = clean_fastq_SE(initial_fastq)

// scaffold flye + 3 rundy wygladzania tym anrzedziem
initial_scaffold = run_flye(processed_fastq)

// jedna runda pilona poki co
final_assembly = pilon_first_nanopore(initial_scaffold, processed_fastq)

} else {
  println("Incorrect option provided")
  System.exit(0)
}

// Post-analizy
// Sprawdzono dzialaja zarowno dla illuminy jak i nnaopore wiec sa poza if-em
extract_final_stats(final_assembly)
MLST_out = run_7MLST(final_assembly)
parse_7MLST(MLST_out)
run_Seqsero(final_assembly)
run_sistr(final_assembly)
run_pointfinder(final_assembly)
run_amrfinder(final_assembly)
run_plasmidfinder(final_assembly)
cgMLST_out = run_cgMLST(final_assembly)
parse_cgMLST_out = parse_cgMLST(cgMLST_out)
run_pHierCC_out = run_pHierCC(parse_cgMLST_out)
extract_historical_data(run_pHierCC_out) 
prokka_out = run_prokka(final_assembly)
run_VFDB(prokka_out)
run_spifinder(final_assembly)

} 

