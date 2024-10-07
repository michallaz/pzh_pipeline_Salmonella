// Input Parameters
params.genus = '' // Can be set up by a User.
                  // This parameter is used only to print a hint to a user that the pipeline will
                  // not produce any valid output for genuses other than Salmonella, Escherichia or Campylobacter

params.reads =''  // Must be set up by User, path to reads i.e. '/mnt/sda1/michall/Salmonella/*_R{1,2}_001.fastq.gz'. Required
params.machine = '' // Can be set to either 'Illumina' or 'Nanopore'. Required


// Minimal quality of a base for kraken2, fastqc 
if ( params.machine  == 'Illumina' ) {
  params.quality = 6
  // Here we will put parameters to determine quliaty of a sample at fidderent stages
  params.min_number_of_reads = 50000 // Like this, if  fastqfile has less than that many reads the QC_status is switched to nie 
                                     // and downstream modules will not produce any valid output
  params.min_median_quality = 10
} else if (params.machine  == 'Nanopore') {
  params.quality = 2
  params.min_number_of_reads = 10000
  params.min_median_quality = 5
} else {
  println("Incorrect sequnecing platform, avalable options are : Illumina and Nanopore")
  System.exit(0)
}

params.cpus = 25 
params.min_coverage_ratio = 0.1 // Minumum coverage of a contig with respect to average coverage for all contigs. Etoki uses 0.2

// Paths to PREDEFINED structure of /db directory
// // Databases for ALL available species
params.AMRFINDER_db_absolute_path_on_host = "/mnt/sda1/michall/db/AMRfider_plus"
params.metaphlan_db_absolute_path_on_host = "/mnt/sda1/michall/db/Metaphlan"
params.kmerfinder_db_absolute_path_on_host = "/mnt/sda1/michall/db/Kmerfinder/kmerfinder_db"
params.kraken2_db_absolute_path_on_host = "/home/michall/kraken2/kraken2_db/kraken2_sdb/"

// Inform user that program is intended for only these 3 genra
if ( params.genus  == 'Salmonella' || params.genus  == 'Escherichia' || params.genus == 'Campylobacter') {
        params.db_absolute_path_on_host="/mnt/sda1/michall/db/"
} else {
        params.db_absolute_path_on_host="/mnt/sda1/michall/db/"
	println("This program is intended to work with following genera: Salmonella, Escherichia, and Campylobacter")
	println("It will continue but unless one of this genera is identified, sub-programs will not execute")
	// System.exit(0)
}


// API TOKEN for ENTEROBASE API 2.0
params.enterobase_api_token = "eyJhbGciOiJIUzI1NiIsImlhdCI6MTcyMDQzNjQxMSwiZXhwIjoxNzM2MjA0NDExfQ.eyJfIjoibUsyNFlZSHd4SyIsInVzZXJuYW1lIjoiTWljaGFsX0xhem5pZXdza2kiLCJpZCI6ODg4MCwiYWRtaW5pc3RyYXRvciI6bnVsbCwiZW1haWwiOiJtbGF6bmlld3NraUBwemguZ292LnBsIiwiYXBpX2FjY2Vzc19jbG9zdHJpZGl1bSI6IlRydWUiLCJhcGlfYWNjZXNzX2Vjb2xpIjoiVHJ1ZSIsImFwaV9hY2Nlc3Nfc2VudGVyaWNhIjoiVHJ1ZSJ9.VEsyVPv8sn1zG7d3uFqEjfk6XFS2qP8P5Y5mh9VPE9w" 


// Kontenery uzywane w tym skrypcie ustawione NA SZTYWNO !!! 
// salmonella_illumina:2.0 - bazowy kontener z programami o kodem
// staphb/prokka:latest - publiczny kontener z prokka, w notatkach mam ze budowanie programu od 0 jest meczace bo to kod sprzed ponad 4 lat 


process run_fastqc_illumina {
  tag "fastqc for sample ${x}"
  container  = 'salmonella_illumina:2.0'
  publishDir "pipeline_wyniki/${x}/QC", mode: 'copy'
  maxForks 5
  input:
  tuple val(x), path(reads), val(QC_STATUS)
  output:
  tuple val(x), path("*csv"), emit: publishdir // Wykresy, kopiujemy bo json wskazuje do publishdir
  tuple val(x), path("*json"), emit: json // Sam json do kopiowania w celu zlozenia "ostatecznego jsona"
  tuple val(x), env(QC_STATUS_EXIT), emit: qcstatus // Sam QC status
  tuple val(x), env(QC_STATUS_EXIT), env(TOTAL_BASES), emit: qcstatus_and_values // QC status + informccje o calkowitej liczbie zasad i ilosci odczytow
  script:
  if (QC_STATUS == null) { QC_STATUS="tak" } // Domyslna wartosc w przypadku gdy user nie poda QC status
  """
  # Set up QC_STATUS to "tak" if a given module does not provide this value via input
  

  # Pipeline bakteryjny ma tylko liczenie QC na etapie 
  DANE_FORWARD=(`python /opt/docker/EToKi/externals/run_fastqc_and_generate_json.py -i ${reads[0]} -m 4048 -c ${params.cpus} -s tak -e pre-filtering -p "pipeline_wyniki/${x}/QC" -o forward.json`)
  STATUS_FORWARD="\${DANE_FORWARD[0]}"
  ILOSC_ODCZYTOW_FORWARD="\${DANE_FORWARD[1]}"
  MEDIANA_JAKOSCI_FORWARD="\${DANE_FORWARD[2]}"
  BASES_FORWARD="\${DANE_FORWARD[3]}"

  if [[ \${STATUS_FORWARD} != "tak" || \${ILOSC_ODCZYTOW_FORWARD} -lt ${params.min_number_of_reads} || \${MEDIANA_JAKOSCI_FORWARD} -lt ${params.min_median_quality} ]]; then
    STATUS_FORWARD_ALL="nie"
  else
    STATUS_FORWARD_ALL="tak"
  fi
 
  DANE_REVERSE=(`python /opt/docker/EToKi/externals/run_fastqc_and_generate_json.py -i ${reads[1]} -m 4048 -c ${params.cpus} -s tak -e pre-filtering -p "pipeline_wyniki/${x}/QC" -o reverse.json`)
  STATUS_REVERSE="\${DANE_REVERSE[0]}"
  ILOSC_ODCZYTOW_REVERSE="\${DANE_REVERSE[1]}"
  MEDIANA_JAKOSCI_REVERSE="\${DANE_REVERSE[2]}"
  BASES_REVERSE="\${DANE_REVERSE[3]}"
 
  TOTAL_BASES=`echo "\${BASES_FORWARD} + \${BASES_REVERSE}" | bc -l`

  if [[ \${STATUS_REVERSE} != "tak" || \${ILOSC_ODCZYTOW_REVERSE} -lt ${params.min_number_of_reads} || \${MEDIANA_JAKOSCI_REVERSE} -lt ${params.min_median_quality} ]]; then
    STATUS_REVERSE_ALL="nie"
  else
    STATUS_REVERSE_ALL="tak"
  fi

  if [[ \${STATUS_FORWARD_ALL} == "nie"  || \${STATUS_REVERSE_ALL} == "nie" ]]; then
    QC_STATUS_EXIT="nie" # moduly "nizej" dostaja niei, bo blad jest na tym etapie
  else
    QC_STATUS_EXIT="tak"
  fi 
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
         --memory 4048 \
         --extract \
         --delete \
         --outdir . \
         ${reads}
    r1=\$(basename ${reads} .fastq.gz)
    /data/parse_fastqc_output.py \${r1}_fastqc/fastqc_data.txt \${QUALITY} >> ${x}_fastqc.txt
  """
}

process clean_fastq_illumina {
  // Prosta funkcja zaimplementowa w etoki do czyszczenia plikow fastq, ma na celu rozdzielenie odczytow ktore sa sparowane
  // od tych ktore pary nie maja, trimmowanie odczytow na podstawie jakosci, zmiana nazwy odczytow .
  // Ta funkcja generalnie mimikuje dzialanie trimmomatic-a. Wiec teoretycznie mozna uzyc jego
  // Default na trimming to quality 6 (-q to zmienia)

  container  = 'salmonella_illumina:2.0'
  tag "Fixing fastq dla sample $x"
  maxForks 5
  input:
  tuple val(x), path(reads), val(QC_status)
  output:
  tuple val(x), path('prep_out_L1_R1.fastq.gz'), path('prep_out_L1_R2.fastq.gz'), emit: PE_path
  tuple val(x), path('prep_out_L1_SE.fastq.gz'), emit: SE_path
  tuple val(x), path('prep_out_L1_R1.fastq.gz'), path('prep_out_L1_R2.fastq.gz'), path('prep_out_L1_SE.fastq.gz'), val(QC_status), emit: All_path
  script:
  read_1 = reads[0]
  read_2 = reads[1]
  """
  if [ ${QC_status} == "nie" ]; then
     touch prep_out_L1_R1.fastq.gz
     touch prep_out_L1_R2.fastq.gz
     touch prep_out_L1_SE.fastq.gz
  else
    python /opt/docker/EToKi/EToKi.py prepare --pe ${read_1},${read_2} -p prep_out -c ${params.cpus} -q ${params.quality}
    # w niektorych przypadkach plik SE moze byc pusty (bo pewnie ktos go wczesniej czyscil) 
    # Etoki sie gubi i nie generuje poprawnie pliku SE
    # tworzymy samodzilenie plik z dummy sekwencja tak by reszta skryptu dzialala
    # NIE JEST TO POWOD DO ZMIANY PARAMETRU QC BO PLIKI pair-end SA POPRAWNE
    if [ -e prep_out.1.0.s.fastq.gz ]; then
      echo "@0_SE_0" >> prep_out_L1_SE.fastq
      echo "GTACTGACCAAACTGGTCGTGTAGCGTTTCATGCCACATCGTATTTTCGGCCATTGGCTGATACCTCCATTGTTAACACCCGTAAAAAAAGGGCGCAACATCATAGCTAACAATGACCGTGGATGCACGGTCATTATTTCAGCAATAGGAT" >> prep_out_L1_SE.fastq
      echo "+" >> prep_out_L1_SE.fastq
      echo "AFFFFFFFFFFFFAFFFFFFFFAFFFAFFF/FAFFAAAFF/=FFAAFFAFFFF/FFFFFF//FFFFFFFFA/FFFFFAFFAA=FFFFFFFF/FFFFFFFFF/FFF/FFFFFFFFAFFFAFFFFFAFAF/FF=FFFFAFFFFF/FF/FAFAF" >> prep_out_L1_SE.fastq
      gzip prep_out_L1_SE.fastq
    fi
  fi
  """
}


process run_initial_mlst_illumina {
  // Pierwszy check MLSDT 7 genomwego dla celow QC
  // z wykorzystaniem soft cge
  // Proces wymaga zarwono golych odczytow, jak i informacji o gatunku, wiec bedzie 'wpiety' po species prediciton
  // Jego QC jest wazne i jego wynik bedzie wpiety jako input do skladania genomu
  // w tym celu wstawimy zmienna QC_status, ktora bedzie 'podrozowac' razem z samplem
  // jesli qc_status bedzie fail to procesy potomne sie nie uruchomia ? zwroca dummy wartosci ?
  container  = 'salmonella_illumina:2.0'
  tag "Initial MLST for sample $x"
  maxForks 5
  input:
  tuple val(x), path(reads), val(SPECIES), val(GENUS), val(QC_status)
  output:
  tuple val(x), path(reads), env(QC_status_exit)
  script:
  read_1 = reads[0]
  read_2 = reads[1]
  """
  mkdir tmp
  if [ ${QC_status} == "nie" ]; then
    QC_status_exit="nie"
  else
    if [[ "${GENUS}" == *"Salmo"* ]]; then
    python /opt/docker/mlst/mlst.py -i ${read_1} ${read_2} -s senterica -p /opt/docker/mlst_db/ -mp kma -t tmp/
    elif [[ "${GENUS}" == *"Escher"* ]]; then
    python /opt/docker/mlst/mlst.py -i ${read_1} ${read_2} -s ecoli -p /opt/docker/mlst_db/ -mp kma -t tmp/
    elif [ ${GENUS} == "Campylobacter" ]; then
    # w tej bazie podgatunki campylo okreslane sa typowo z cjejuni, clari itd .. 
    python /opt/docker/mlst/mlst.py -i ${read_1} ${read_2} -s c${SPECIES} -p /opt/docker/mlst_db/ -mp kma -t tmp/
    fi

    # Parsowanie wyniku

    # Warunki do QC 
    QC_status_exit="tak"
  fi
  """
}



process run_initial_mlst_nanopore {
  // Proces dla nanopore
  container  = 'salmonella_illumina:2.0'
  tag "Initial MLST for sample $x"
  maxForks 5
  input:
  tuple val(x), path(reads), val(SPECIES), val(GENUS)
  output:
  env(QC_status)
  script:
  """
  mkdir tmp
  if [[ "${GENUS}" == *"Salmo"* ]]; then
  python /opt/docker/mlst/mlst.py -i ${reads} -s senterica -p /opt/docker/mlst_db/ -mp kma -t tmp/
  elif [[ "${GENUS}" == *"Escher"* ]]; then
  python /opt/docker/mlst/mlst.py -i ${reads} -s ecoli -p /opt/docker/mlst_db/ -mp kma -t tmp/
  elif [ ${GENUS} == "Campylobacter" ]; then
  # w tej bazie podgatunki campylo okreslane sa typowo z cjejuni, clari itd ..
  python /opt/docker/mlst/mlst.py -i ${reads} -s c${SPECIES} -p /opt/docker/mlst_db/ -mp kma -t tmp/
  fi

  # Parsowanie wyniku

  # Warunki do QC
  QC_status="pass"
  """
}


process spades {
  // Funkcja do odpalania spadesa
  // Powstaly plik fasta jest poprawiany bo Hapo-G nie akceptuje "." w nazwach sekwencji w pliku fasta
  // Modul definiuje QC_status az do modulu extract_final_stats, gdzie status moze ulec zmianie
  cpus params.cpus
  container  = 'salmonella_illumina:2.0'
  tag "Spades dla sample $x"
  maxForks 5
  input:
  tuple val(x), path('R1.fastq.gz'), path('R2.fastq.gz'), path('SE.fastq.gz'), val(QC_status)
  output:
  tuple val(x),  path('scaffolds_fix.fasta'), env(QC_status)
  script:
  """
  #  for 150bp reads SPAdes uses k-mer sizes 21, 33, 55, 77
  #  We strongly recommend not to change -k parameter unless you are clearly aware about the effect.
  # --isolate This flag is highly recommended for high-coverage isolate and multi-cell Illumina data;
  # --careful Tries to reduce the number of mismatches and short indels. ale nie dziala z isolate
  #  wywolanie spades w ramach etoki to tylko definicja inputu + liczby procesorow
  # /opt/docker/EToKi/externals/spades.py -t 8 --pe-1 1 /Salomenlla/test_etoki/id_151_novel_etoki_with_comments/prep_out_L1_R1.fastq.gz --pe-2 1 /Salomenlla/test_etoki/id_151_novel_etoki_with_comments/prep_out_L1_R2.fastq.gz --pe-s 2 /Salomenlla/test_etoki/id_151_novel_etoki_with_comments/prep_out_L1_SE.fastq.gz -o spades

  if [ ${QC_status} == "nie" ]; then
    # We create dummy files so that pipeline can continue
    echo ">dummy_contig" >> scaffolds_fix.fasta
    echo "AAAAAAAAAAAAA" >> scaffolds_fix.fasta
    # zwracanie json z informacja ze na tym etapie pipeline zakoczyl dzialanie z powodu zbyt malej liczby odczytow
    # wszystkie podrogramy przecwytuja komentarz i zwracaja odpowidni komunikat 
  else
    python /opt/docker/EToKi/externals/spades.py --isolate -t ${task.cpus} --pe-1 1 R1.fastq.gz --pe-2 1 R2.fastq.gz --pe-s 2 SE.fastq.gz  -o spades_manual
    cat spades_manual/scaffolds.fasta  | awk '{if(\$0 ~ />/) {split(\$0, ala, "_"); print ">NODE"ala[2]} else print \$0}' >> scaffolds_fix.fasta
  
  fi
  """
}


process bwa_paired {
  // Funkcja do mapowania odczytow Pair-end
  cpus params.cpus
  container  = 'salmonella_illumina:2.0'
  tag "RE-mapowanie PE dla sample $x"
  maxForks 5
  input:
  tuple val(x), path('genomic_fasta.fasta'), val(QC_status), path(read_1),  path(read_2)
  output:
  tuple val(x), path('mapowanie_bwa_PE.bam')
  script:
  """

  if [ ${QC_status} == "nie" ]; then
      touch mapowanie_bwa_PE.bam
  else
    /opt/docker/EToKi/externals/bwa index genomic_fasta.fasta
    /opt/docker/EToKi/externals/bwa mem -t ${task.cpus} -T 30 genomic_fasta.fasta ${read_1} ${read_2} |  /opt/docker/EToKi/externals/samtools fixmate -m -@ ${task.cpus} - - | /opt/docker/EToKi/externals/samtools sort -@ ${task.cpus} -  |  /opt/docker/EToKi/externals/samtools markdup -r -@ ${task.cpus}  -O BAM - mapowanie_bwa_PE.bam
  fi
  """
}


process bwa_single {
  // Funkcja do mapowania odczytow Single-end na genom
  cpus params.cpus
  container  = 'salmonella_illumina:2.0'
  tag "RE-mapowanie SE dla sample $x"
  maxForks 5
  input:
  tuple val(x), path('genomic_fasta.fasta'), val(QC_status), path(reads)
  output:
  tuple val(x), path('mapowanie_bwa_SE.bam')
  // bwa_single przekazuje QC_status do merge_bams, nie ma potrzeby aby robily to oba moduly do bwa
  script:
  """
  if [ ${QC_status} == "nie" ]; then
     touch mapowanie_bwa_SE.bam
  else
    /opt/docker/EToKi/externals/bwa index genomic_fasta.fasta
    /opt/docker/EToKi/externals/bwa mem -t ${task.cpus} -T 30 genomic_fasta.fasta ${reads} |  /opt/docker/EToKi/externals/samtools sort -@ ${task.cpus} -  |  /opt/docker/EToKi/externals/samtools markdup -r -@ ${task.cpus} -O BAM - mapowanie_bwa_SE.bam
  fi
  """
}

process merge_bams {
  // process do mergowania bam-ow
  // Uzywany w workflow calculate_coverage
  container  = 'salmonella_illumina:2.0'
  tag "Merging bam files for sample $x"
  input:
  tuple val(x), path(bam1), path(bam2)
  output:
  tuple val(x), path('merged_bams.bam')
  script:
  """
  PAIRED_BAM_SIZE=`wc -l $bam1 | cut -d " " -f1`
  if [ \${PAIRED_BAM_SIZE} -eq 0 ]; then
     touch merged_bams.bam
  else
    /opt/docker/EToKi/externals/samtools merge -f merged_bams.bam $bam1 $bam2
  fi
  """ 
}

process run_pilon {
  // Dokladne uzycie pilona jest tu https://github.com/broadinstitute/pilon/wiki/Requirements-&-Usage
  // Pilon tylko przekazuje QC_status bez mozliwosci zmiany
  container  = 'salmonella_illumina:2.0'
  tag "Pilon for sample $x"
  maxForks 5
  // publishDir "pilon_wyniki/${x}", mode: 'copy'
  input:
  tuple val(x), path(bam1), path(bam2), path('genomic_fasta.fasta'), val(QC_status)
  output:
  tuple val(x), path('latest_pilon.fasta'), path('latest_pilon.changes'), emit: ALL
  tuple val(x), path('latest_pilon.fasta'), val(QC_status), emit: ONLY_GENOME
  script:
  """
  if [ ${QC_status} == "nie" ]; then
    echo ">dummy_contig" >> latest_pilon.fasta
    echo "AAAAAAAAAAAAA" >> latest_pilon.fasta
    touch latest_pilon.changes  
  else 
    # indeksacja bam-ow
    /opt/docker/EToKi/externals/samtools index  $bam1
    /opt/docker/EToKi/externals/samtools index  $bam2

    # wywolanie pilon-a
    java -jar /opt/docker/EToKi/externals/pilon.jar --threads ${params.cpus} --genome genomic_fasta.fasta --frags $bam1 --unpaired $bam2 --output pilon_polish --vcf --changes --outdir pilon_bwa
 
    # przygotowanie outputu
    # cat pilon_bwa/pilon_polish.fasta | awk  -v ALA=${x} 'BEGIN{OFS=""}; {if(\$0 ~ />/) print \$0,"_", ALA; else print \$0}' >> last_pilon.fasta
    cat pilon_bwa/pilon_polish.fasta >> latest_pilon.fasta
    cp  pilon_bwa/pilon_polish.changes latest_pilon.changes
  fi
  """
}


process extract_final_contigs {
  // Liczenie pokrycia dla kazdego contiga polaczona z filtorwanie odczytow
  // Przy uzyciu NASZEGO skryptu
  // Modukl przekazuje bez mozliwosci zmiany QC_stauts
  container  = 'salmonella_illumina:2.0'
  tag "Coverage-based filtering for sample $x"
  publishDir "pipeline_wyniki/${x}", mode: 'copy'
  maxForks 5
  input:
  tuple val(x), path(bam1), path('genomic_fasta.fasta'), val(QC_status)
  output:
  tuple val(x), path('final_scaffold_filtered.fa'), val(QC_status), emit: ONLY_GENOME
  tuple val(x), path('final_scaffold_filtered.fa'), path('Rejected_contigs.fa'), val(QC_status), emit: ONLY_GENOME_AND_REJECT
  // final_scaffold_filtered.fa to nazwa ustawiona NA SZTYWNO w skrypcie coverage_filter.py
  script:
  """
  if [ ${QC_status} == "nie" ]; then
    touch final_scaffold_filtered.fa
    touch Rejected_contigs.fa   
  else
    /opt/docker/EToKi/externals/samtools index  $bam1
    python  /opt/docker/EToKi/externals/coverage_filter.py genomic_fasta.fasta $bam1 ${params.min_coverage_ratio}
  fi
  """
}

process extract_final_stats {
  container  = 'salmonella_illumina:2.0'
  tag "Calculating basic statistics for sample $x"
  publishDir "pipeline_wyniki/${x}", mode: 'copy'
  input:
  tuple val(x), path(fasta), path(fasta_reject), val(QC_status)
  output:
  tuple val(x), path('Summary_statistics.txt'), path('Summary_statistics_with_reject.txt'), env(QC_status), emit: STATISTICS
  tuple val(x), path(fasta), env(QC_status), emit: GENOME
  script:
  """
  if [ ${QC_status} == "nie" ]; then
     touch Summary_statistics.txt
     touch Summary_statistics_with_reject.txt
     QC_status="nie" # jako ze output to env, musimy w srodowisku ustalic ta zmienne inaczej program nie "przechwyci" do output QC_status
  else
    cat $fasta $fasta_reject >> all_contigs.fasta
    python  /opt/docker/EToKi/externals/calculate_stats.py $fasta all_contigs.fasta
    # Tutaj parsujemy output Summary_statistics.txt aby ocenic czy genom ma odpowiednia dlugosc, ilosc contigow nie jest absurdalna itd
    # Wiec tu jest potencjalny switch z pass do fail
    QC_status="tak"
  fi
  """
}

process run_7MLST {
  container  = 'salmonella_illumina:2.0'
  containerOptions "--volume ${params.db_absolute_path_on_host}:/db"
  tag "Predicting MLST for sample $x"
  // publishDir "pipeline_wyniki/${x}", mode: 'copy'
  input:
  tuple val(x), path(fasta), val(QC_status), val(SPECIES), val(GENUS), val(QC_status_contaminations)
  output:
  tuple val(x), path('MLSTout.txt'), val(SPECIES), val(GENUS), val(QC_status), val(QC_status_contaminations)
  // when:
  // GENUS == 'Salmonella' || GENUS == 'Escherichia' || GENUS == 'Campylobacter'
  script:
  """
  if [[ ${QC_status} == "nie"  || ${QC_status_contaminations} == "nie" ]]; then
    touch MLSTout.txt
    # json na zle QC
  else
    CAMPYLO_SPECIES='concisus fetus helveticus hyointestinalis insulaenigrae jejuni lanienae lari sputorum upsaliensis'
    if [[ "${SPECIES}" == *"Salmo"* || "${SPECIES}" == *"Escher"* ]]; then
      /opt/docker/EToKi/EToKi.py MLSTdb -i /db/${GENUS}/MLST_Achtman/all_allels.fasta -x 0.8 -m 0.5 -r MLST_Achtman_ref.fasta -d MLST_database.tab
      /opt/docker/EToKi/EToKi.py MLSType -i $fasta -r MLST_Achtman_ref.fasta -k ${x} -o MLSTout.txt -d MLST_database.tab
    elif [[ \${CAMPYLO_SPECIES[@]} =~ "${SPECIES}" ]]; then
      # sciezka dla Campylobacter
      /opt/docker/EToKi/EToKi.py MLSTdb -i /db/${GENUS}/${SPECIES}/MLST/all_allels.fasta -x 0.8 -m 0.5 -r MLST_Achtman_ref.fasta -d MLST_database.tab
     /opt/docker/EToKi/EToKi.py MLSType -i $fasta -r MLST_Achtman_ref.fasta -k ${x} -o MLSTout.txt -d MLST_database.tab
    else
       echo "Provided species ${SPECIES} is not part of any MLST databases" >> MLSTout.txt
       # json na zly gatunek
    fi # koniec if-a na zly gatunek
  fi # koniec if-a na zle QC
  """
}


process parse_7MLST {
  // Process that extract a ST given a set of loci calculated with run_7MLST process
  // Process returns to 3 files 
  
  // // One MLST_parsed_output.txt with four columns (1st is ST identified for a sample, 2nd is the closest ST fromamong known ST in database. 
  // // 3nd is the distance between a sample allelic profile and allelic profile of closes matching ST, can be between
  // //  0 to 7, locus with no allelic varaiant are omitted); 4th column is a Comment
  // // If distance is 0 1st and 2nd column should be identical
  
  // // Second file MLST_sample_full_list_of_allels.txt contains eight columns (1st is the ST identified for a sample, followed by allel versions for each locus in a given scheme)
  
  // // Third file is MLST_closest_ST_full_list_of_allels.txt contains eight columns (1st is the ST identified for a sample, followed by allel versions for each locus of that ST)
  // // If the distance between between identified and expexted profiles is 0 (As indicated in MLST_parsed_output.txt) this file should be identical to MLST_sample_full_list_of_allels.txt

  container  = 'salmonella_illumina:2.0'
  containerOptions "--volume ${params.db_absolute_path_on_host}:/db"
  tag "Parsing MLST for sample $x"
  publishDir "pipeline_wyniki/${x}/MLST", mode: 'copy', pattern: 'MLST*txt'
  maxForks 1
  input:
  tuple val(x), path('MLSTout.txt'), val(SPECIES), val(GENUS), val(QC_status), val(QC_status_contaminations) 
  output:
  tuple val(x), path('MLST_parsed_output.txt'), path('MLST_sample_full_list_of_allels.txt'), path('MLST_closest_ST_full_list_of_allels.txt')
  // when:
  // GENUS == 'Salmonella' || GENUS == 'Escherichia' || GENUS == 'Campylobacter'
  script:
"""
#!/usr/bin/python
import sys
import re

sys.path.append('/data')
from all_functions_salmonella import *

species="$SPECIES"
genus="$GENUS"
qc_status="$QC_status"
qc_status_contaminations="$QC_status_contaminations"

### Determine correct paths to /db
### Prepare dummy output if provided SPECIES is not handled

if qc_status == "nie" or qc_status_contaminations == "nie":
    with open('MLST_parsed_output.txt', 'w') as f1, open('MLST_sample_full_list_of_allels.txt', 'w') as f2, open('MLST_closest_ST_full_list_of_allels.txt', 'w') as f3:
        f1.write(f'ST\\tComment\\n')
        f1.write(f'unk\\tUnknown species: {species}\\n')

        f2.write(f'ST\\tComment\\n')
        f2.write(f'unk\\tUnknownn species: {species}')

        f3.write(f'ST\\tComment\\n')
        f3.write(f'unk\\tUnknown species: {species}')
    # json na zle QC
    sys.exit(0)

if re.findall('Salmo', genus) or re.findall('Esch', genus):
    sciezka=f'/db/{genus}/MLST_Achtman'
elif species in ['concisus','fetus','helveticus','hyointestinalis','insulaenigrae','jejuni','lanienae','lari','sputorum','upsaliensis']:
    sciezka=f'/db/{genus}/{species}/MLST'
else:
    with open('MLST_parsed_output.txt', 'w') as f1, open('MLST_sample_full_list_of_allels.txt', 'w') as f2, open('MLST_closest_ST_full_list_of_allels.txt', 'w') as f3:
        f1.write(f'ST\\tComment\\n')
        f1.write(f'unk\\tUnknown species: {species}\\n')
       
        f2.write(f'ST\\tComment\\n')
        f2.write(f'unk\\tUnknownn species: {species}')

        f3.write(f'ST\\tComment\\n')
        f3.write(f'unk\\tUnknown species: {species}')
    # json na zly gatunek
    sys.exit(0)
### Deterine ST of out Sample
known_profiles, klucze_sorted = create_profile(f'{sciezka}/profiles.list')
identified_profile = parse_MLST_fasta('MLSTout.txt')
matching_ST, min_value, matching_ST_profile, sample_profile, all_loci_names = getST(MLSTout = identified_profile, \
                                                                              profile_file = f'{sciezka}/profiles.list')


### matching_ST is a string
### matching_ST_profile, sample_profile, all_loci_names are all tab-separated strings
### min_valu is int

### Prep output

with open('MLST_closest_ST_full_list_of_allels.txt', 'w') as f3:
    f3.write(f'ST\\t{all_loci_names}\\n')
    f3.write(f'{matching_ST}\\t{matching_ST_profile}\\n')

with open('MLST_parsed_output.txt', 'w') as f1, open('MLST_sample_full_list_of_allels.txt', 'w') as f2:
    if min_value == 0:
        f1.write(f'ST_sample\\tST_database\\tDistance\\tComment\\n')
        f1.write(f'{matching_ST}\\t{matching_ST}\\t{min_value}\\tIn external database\\n')

        f2.write(f'ST\\t{all_loci_names}\\n')
        f2.write(f'{matching_ST}\\t{sample_profile}\\n')

    else:
        # Unknown profile
        # Check if this profile is in a "local" database and if not create a new entry
        known_profiles_local, klucze_sorted_local = create_profile(f'{sciezka}/local/profiles_local.list')
        matching_ST_local, min_value_local, matching_ST_profile_local, _, _ =  getST(MLSTout = identified_profile, \
                                                                               profile_file = f'{sciezka}/local/profiles_local.list')
    
        if min_value_local == 0:
            # found profile in local database
            # we still print info what is the distance to known ST in external database
            f1.write(f'ST_sample\\tST_database\\tDistance\\tComment\\n')
            f1.write(f'{matching_ST_local}\\t{matching_ST}\\t{min_value}\\tIn local database\\n')

            f2.write(f'ST\\t{all_loci_names}\\n')
            f2.write(f'{matching_ST_local}\\t{sample_profile}\\n')

        else:
            # new profile appending it to local database with new number  
            last_ST = 0
            with open(f'{sciezka}/local/profiles_local.list') as f:
                for line in f:
                    line = line.rsplit()
                    if 'local' in line[0]:
                        last_ST = line[0].split('_')[1]
            novel_profile_ST = f'local_{int(last_ST)+1}'

            write_novel_sample(f'{novel_profile_ST}\\t{sample_profile}\\n', f'{sciezka}/local/profiles_local.list')
            
            f1.write(f'ST_sample\\tST_database\\tDistance\\tComment\\n')
            f1.write(f'{novel_profile_ST}\\t{matching_ST}\\t{min_value}\\tNovel in local database\\n')

            f2.write(f'ST\\t{all_loci_names}\\n')
            f2.write(f'{novel_profile_ST}\\t{sample_profile}\\n')

"""
}

process run_Seqsero {
  // SeqSero only works for Salmonella
  container  = 'salmonella_illumina:2.0'
  tag "Predicting OH for sample $x with Seqsero"
  publishDir "pipeline_wyniki/${x}", mode: 'copy'
  cpus params.cpus
  maxForks 15
  input:
  tuple val(x), path(fasta), val(QC_status), val(SPECIES), val(GENUS), val(QC_status_contaminations)
  output:
  tuple val(x), path('seqsero/SeqSero_result.txt')
  // when:
  // GENUS == 'Salmonella'
  script:
  """
  # -m to rodzaj algorytmu -m to chyba opart o k-mery
  # -t 4 to informacja ze inputem sa contigi z genomem
  # -p to procki, proces jest szybki wiec ustawie 4 + maxforks 15

  if [[ ${QC_status} == "nie"  || ${QC_status_contaminations} == "nie" ]]; then
    mkdir seqsero
    touch seqsero/SeqSero_result.txt
    #json na zle QC
  else
    if [ ${GENUS} == "Salmonella" ]; then
      python /opt/docker/SeqSero2/bin/SeqSero2_package.py -m k -t 4 -p 4 -i $fasta -d seqsero
    else
      mkdir seqsero
      touch seqsero/SeqSero_result.txt
      #json na zly gatunek
    fi # koniec if-a na zly gatunek
  fi # koniec if-a na zle QC
  """
}

process run_sistr {
  // Sistr works only for Salmonella
  container  = 'salmonella_illumina:2.0'
  tag "Predicting OH for sample $x with Sistr"
  publishDir "pipeline_wyniki/${x}/sistr", mode: 'copy'
  cpus params.cpus
  maxForks 15
  input:
  tuple val(x), path(fasta), val(QC_status), val(SPECIES), val(GENUS), val(QC_status_contaminations)
  output:
  tuple val(x), path('sistr-output.tab')
  // when:
  // GENUS == 'Salmonella'
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
 if [[ ${QC_status} == "nie"  || ${QC_status_contaminations} == "nie" ]]; then
    # json na zle QC
    touch sistr-output.tab
  else
    if [ ${GENUS} == "Salmonella" ]; then
      /usr/local/bin/sistr --qc -vv --alleles-output allele-results.json --novel-alleles novel-alleles.fasta --cgmlst-profiles cgmlst-profiles.csv -f tab -t 4 -o sistr-output.tab $fasta
    else
      # json na zly gatunek
      touch sistr-output.tab
    fi # koniec if-a na zly gatunek
  fi # koniec if-a na zle QC
  """
}

process run_ectyper {
  // This process works only for Escherichia
  container  = 'salmonella_illumina:2.0'
  tag "Predicting OH for sample $x with ectyper"
  publishDir "pipeline_wyniki/${x}/ECTyper", mode: 'copy'
  cpus params.cpus
  maxForks 15
  input:
  tuple val(x), path(fasta), val(QC_status), val(SPECIES), val(GENUS), val(QC_status_contaminations)
  output:
  tuple val(x), path('output.tsv')
  // when:
  // GENUS == 'Escherichia'
  script:
  """
  # opcje:
  # -i to input, 
  # -c to liczba core'ow
  # -hpid to minimalny seq identity dla antygenu H ustawiam na 90 zamiast default (95) po analizie Strain-u 5 z EQA 2023
  # -o to katalog z output
  mkdir ectyper_out

  if [[ ${QC_status} == "nie"  || ${QC_status_contaminations} == "nie" ]]; then
    # json na zle QC
    touch output.tsv
  else
    if [ ${GENUS} == "Escherichia" ]; then
      ectyper -i $fasta -c 4 -hpid 90 -o ectyper_out
      cp ectyper_out/output.tsv .
    else
       #json na zly gatunek
       touch output.tsv
    fi # koniec if-a na zly gatunek
  fi # koniec if-a na zle QC
  """
}


process run_resfinder {
  container  = 'salmonella_illumina:2.0'
  tag "Predicting microbial resistance for sample $x"
  publishDir "pipeline_wyniki/${x}/", mode: 'copy', pattern: "resfinder/ResFinder_results_table.txt"
  publishDir "pipeline_wyniki/${x}/", mode: 'copy', pattern: "resfinder/pheno_table_*.txt"
  publishDir "pipeline_wyniki/${x}/", mode: 'copy', pattern: "resfinder/PointFinder_results.txt"
  input:
  tuple val(x), path(fasta), val(QC_status), val(SPECIES), val(GENUS), val(QC_status_contaminations)
  output:
  tuple val(x), path('resfinder/pheno_table*.txt'), path('resfinder/ResFinder_results_table.txt'), path('resfinder/PointFinder_results.txt')
  // when:
  // GENUS == 'Salmonella' || GENUS == 'Escherichia' || GENUS == 'Campylobacter'
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
  // pozostale opcje to sciezki do baz /instalowanych wraz z tworzeniem obrazu/, zastanowic sie nad "wypchnieciem ich na zewnatrz" 
  """
  
  # resfinder-owi mozna tez podac pliki fastq (-ifq)
  # resfinder operated on species-level, but for now we asume here that all Salmonella are s.enterica, or all Campylobaster are c.jejuni
  # even if different species is actually analyzed. 
  
  if [[ ${QC_status} == "nie"  || ${QC_status_contaminations} == "nie" ]]; then
    # Tworzenie json i output 
    mkdir resfinder
    cd resfinder
    touch pheno_table_1.txt
    touch ResFinder_results_table.txt
    touch PointFinder_results.txt
  else

    if [[ "${GENUS}" == *"Salmo"* ]]; then
        python -m resfinder -o resfinder/ -s 'senterica'  -l 0.6 -t 0.8 --acquired --point -k /opt/docker/kma/kma -db_disinf /opt/docker/disinfinder_db/ -db_res /opt/docker/resfinder_db/ -db_point /opt/docker/pointfinder_db/ -ifa ${fasta}
    elif [[ "${GENUS}" == *"Escher"* ]]; then
        python -m resfinder -o resfinder/ -s 'ecoli'  -l 0.6 -t 0.8 --acquired --point -k /opt/docker/kma/kma -db_disinf /opt/docker/disinfinder_db/ -db_res /opt/docker/resfinder_db/ -db_point /opt/docker/pointfinder_db/ -ifa ${fasta}
    elif [ ${GENUS} == "Campylobacter" ]; then
         python -m resfinder -o resfinder/ -s 'cjejuni'  -l 0.6 -t 0.8 --acquired --point -k /opt/docker/kma/kma -db_disinf /opt/docker/disinfinder_db/ -db_res /opt/docker/resfinder_db/ -db_point /opt/docker/pointfinder_db/ -ifa ${fasta}
    else
      # Tworzenie json i dummy output
      mkdir resfinder
      cd resfinder
      touch pheno_table_1.txt
      touch ResFinder_results_table.txt
      touch PointFinder_results.txt
    fi # koniec if-a na gatunek przy zdanym QC 
 
 fi # koniec if-a na status QC przekazany modulowi
 
  """
}

process run_cgMLST {
  container  = 'salmonella_illumina:2.0'
  containerOptions "--volume ${params.db_absolute_path_on_host}:/db"
  tag "Predicting cgMLST for sample $x"
  // publishDir "pipeline_wyniki/${x}", mode: 'copy'
  cpus params.cpus
  maxForks 5
  input:
  tuple val(x), path(fasta), val(QC_status), val(SPECIES), val(GENUS), val(QC_status_contaminations)
  output: 
  tuple val(x), path('cgMLST.txt'), path('cgMLST_all_identical_allels.txt'), val(SPECIES), val(GENUS), val(QC_status), val(QC_status_contaminations)
  // when:
  // GENUS == 'Salmonella' || GENUS == 'Escherichia' || SPECIES == 'jejuni'
  // among Campylobacter only c.jejuni has cgMLST scheme
  script:
  """
  if [[ ${QC_status} == "nie"  || ${QC_status_contaminations} == "nie" ]]; then
    touch cgMLST.txt; touch cgMLST_all_identical_allels.txt
    # json dla zlego QC
  else
    if [[ "${GENUS}" == *"Salmo"* ]]; then 
         /data/run_blastn_ver11.sh $fasta ${task.cpus} /db/${GENUS}/cgMLST_v2
    elif [[ "${GENUS}" == *"Esche"* ]]; then
         /data/run_blastn_ver11.sh $fasta ${task.cpus} /db/${GENUS}/cgMLST_v1
    elif [ "${SPECIES}" == "jejuni" ]; then
         /data/run_blastn_ver11.sh $fasta ${task.cpus} /db/${GENUS}/jejuni/cgMLST_v2
    else
         # This should never happen
         echo "Provided species $SPECIES is not part of any cgMLST databases" >> log.log
    fi # koniec if-a na zly gatunek
    cat log.log | cut -f1,2 > cgMLST.txt
  fi # koniec if-a na zle QC
  """
}

process parse_cgMLST {
  // Extracting ST given cgMLST profile, the output is identical to that from parse_7MLST
  container  = 'salmonella_illumina:2.0'
  containerOptions "--volume ${params.db_absolute_path_on_host}:/db"
  tag "Parsing cgMLST for sample $x"
  publishDir "pipeline_wyniki/${x}/cgMLST", mode: 'copy', pattern: 'cgMLST*txt'
  maxForks 1 // set to "1" thus we ensure that when multiple sequencing are analyzed we correctly assign cgST for each of the sample (if they are all new 
  // and not part of the Enterobase
  input:
  tuple val(x), path('cgMLST.txt'), path('cgMLST_all_identical_allels.txt'), val(SPECIES), val(GENUS), val(QC_status), val(QC_status_contaminations)
  output:
  tuple val(x), path('cgMLST_parsed_output.txt'), path('cgMLST_sample_full_list_of_allels.txt'), path('cgMLST_closest_ST_full_list_of_allels.txt'), val(SPECIES), val(GENUS), val(QC_status), val(QC_status_contaminations), emit: standard
  path('cgMLST_all_identical_allels.txt'), emit: output
  // when:
  // GENUS == 'Salmonella' || GENUS == 'Escherichia' || SPECIES == 'jejuni'
  script:
"""
#!/usr/bin/python
import  sys
import re
sys.path.append('/data')
from all_functions_salmonella import *

species="$SPECIES"
genus="$GENUS"
qc_status="$QC_status"
qc_status_contaminations="$QC_status_contaminations"

if qc_status == "nie" or qc_status_contaminations == "nie":
    with open('cgMLST_parsed_output.txt', 'w') as f1, open('cgMLST_sample_full_list_of_allels.txt', 'w') as f2, open('cgMLST_closest_ST_full_list_of_allels.txt', 'w') as f3:
        f1.write(f'ST\\tComment\\n')
        f1.write(f'unk\\tUnknown species: {species}\\n')

        f2.write(f'ST\\tComment\\n')
        f2.write(f'unk\\tUnknownn species: {species}')

        f3.write(f'ST\\tComment\\n')
        f3.write(f'unk\\tUnknown species: {species}')
    # json na zle QC
    sys.exit(0)


if genus == 'Salmonella':
    sciezka=f'/db/{genus}/cgMLST_v2'
elif genus == 'Escherichia':
    sciezka=f'/db/{genus}/cgMLST_v1'
elif species == 'jejuni':
    sciezka=f'/db/{genus}/jejuni/cgMLST_v2'
else:
    with open('cgMLST_parsed_output.txt', 'w') as f1, open('cgMLST_sample_full_list_of_allels.txt', 'w') as f2, open('cgMLST_closest_ST_full_list_of_allels.txt', 'w') as f3:
        f1.write(f'ST\\tComment\\n')
        f1.write(f'unk\\tUnknown species: {species}\\n')

        f2.write(f'ST\\tComment\\n')
        f2.write(f'unk\\tUnknownn species: {species}')

        f3.write(f'ST\\tComment\\n')
        f3.write(f'unk\\tUnknown species: {species}')
    # json na zly gatunek
    sys.exit(0) 

identified_profile = parse_MLST_blastn('cgMLST.txt')

matching_ST, min_value,  matching_ST_profile, sample_profile, all_loci_names = getST(MLSTout = identified_profile, \
                                                                                     profile_file = f'{sciezka}/profiles.list')

# Write info regarding closest ST from EXTERNAL database
with open('cgMLST_closest_ST_full_list_of_allels.txt', 'w') as f3:
    f3.write(f'ST\\t{all_loci_names}\\n')
    f3.write(f'{matching_ST}\\t{matching_ST_profile}\\n')

with open('cgMLST_parsed_output.txt', 'w') as f1, open('cgMLST_sample_full_list_of_allels.txt', 'w') as f2:

    # for f1 and f2 the only difference is ST assigned to the sample
        
    if min_value == 0:
        f1.write(f'ST_sample\\tST_database\\tDistance\\tComment\\n')
        f1.write(f'{matching_ST}\\t{matching_ST}\\t{min_value}\\tIn external database\\n')

        f2.write(f'ST\\t{all_loci_names}\\n')
        f2.write(f'{matching_ST}\\t{sample_profile}\\n')
   
    else:
        # look for a profile in "local" database 
        matching_ST_local, min_value_local, matching_ST_profile_local, _, _ = getST(MLSTout = identified_profile,
                                                                                    profile_file = f'{sciezka}/local/profiles_local.list')
        if min_value_local == 0:
            f1.write(f'ST_sample\\tST_database\\tDistance\\tComment\\n')
            f1.write(f'{matching_ST_local}\\t{matching_ST}\\t{min_value}\\tIn local database\\n')

            f2.write(f'ST\\t{all_loci_names}\\n')
            f2.write(f'{matching_ST_local}\\t{sample_profile}\\n')

        else:
            # new allelic profile not present in either external or local databases
            last_ST = 0
            with open(f'{sciezka}/local/profiles_local.list') as f:
                for line in f:
                    line = line.rsplit()
                    if 'local' in line[0]:
                        last_ST = line[0].split('_')[1]
        
        
            novel_profile_ST = f'local_{int(last_ST)+1}'
            
            write_novel_sample(f'{novel_profile_ST}\\t{sample_profile}\\n', f'{sciezka}/local/profiles_local.list')
        
            f1.write(f'ST_sample\\tST_database\\tDistance\\tComment\\n')
            f1.write(f'{novel_profile_ST}\\t{matching_ST}\\t{min_value}\\tNovel in local database\\n')

            f2.write(f'ST\\t{all_loci_names}\\n')
            f2.write(f'{novel_profile_ST}\\t{sample_profile}\\n')

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
 
  // I assume there is no point checking here what is the organism

  container  = 'staphb/prokka:latest'
  tag "Predicting genes for sample $x"
  // Do we really need prokka output in results dir ? 
  publishDir "pipeline_wyniki/${x}", mode: 'copy'
  cpus params.cpus
  maxForks 5
  input:
  tuple val(x), path(fasta), val(QC_status), val(SPECIES), val(GENUS), val(QC_status_contaminations)
  output:
  tuple val(x), path('prokka_out/prokka_out*gff'), path('prokka_out/*faa'), path('prokka_out/*.ffn'), path('prokka_out/*.tsv'), val(SPECIES), val(GENUS), val(QC_status), val(QC_status_contaminations)
  //when:
  //GENUS == 'Salmonella' || GENUS == 'Escherichia' || GENUS == 'Campylobacter'
  script:
  """
  if [[ ${QC_status} == "nie"  || ${QC_status_contaminations} == "nie" ]]; then
    mkdir prokka_out; touch prokka_out/prokka_out_dummy.gff; touch prokka_out/prokka_out_dummy.faa; touch prokka_out/prokka_out_dummy.ffn; touch prokka_out/prokka_out_dummy.tsv
    # json z informacja o bledzie jakosci
  else
    if [[ ${GENUS} == "Salmonella" || ${GENUS} == "Escherichia" || ${GENUS} == "Campylobacter" ]]; then
      prokka --metagenome --cpus ${params.cpus} --outdir prokka_out --prefix prokka_out --compliant --kingdom Bacteria $fasta 
    else
      mkdir prokka_out; touch prokka_out/prokka_out_dummy.gff; prokka_out/prokka_out_dummy.ffa; prokka_out/prokka_out_dummy.ffn; prokka_out/prokka_out_dummy.tsv
      # json z informacja o zlym gatunku

    fi
  fi
  """
}

process run_VFDB {
  // Baza z czynnikami wirulencji w roznych bakteriach w tym salmonelli
  // Przygotowana baza z instrukcja jak ja przygotowac jest w /mnt/sda1/michall/db/VFDB/README
  container  = 'salmonella_illumina:2.0'
  containerOptions '--volume /mnt/sda1/michall/db/VFDB:/db'
  // Ponownie montujemy na sztyno do /db bo taka lokalizacje na sztywno ma wpisany moj skrypt
  tag "Predicting VirulenceFactors for sample $x"
  publishDir "pipeline_wyniki/${x}/VFDB/", mode: 'copy'
  cpus params.cpus
  maxForks 5
  input:
  // inputem jest output procesu run_prokka
  tuple val(x), path(gff), path(faa), path(ffn), path(tsv), val(SPECIES), val(GENUS), val(QC_status), val(QC_status_contaminations)
  output:
  tuple val(x), path('VFDB_summary*txt'), val(SPECIES), val(GENUS), emit: non_ecoli
  tuple val(x), path('VFDB_summary_Escherichia.txt'), path('VFDB_summary_Shigella.txt'), val(SPECIES), val(GENUS), val(QC_status), val(QC_status_contaminations), optional: true, emit: ecoli
  //when:
  //GENUS == 'Salmonella' || GENUS == 'Escherichia' || GENUS == 'Campylobacter'
  script:
  """
  
  if [[ ${QC_status} == "nie"  || ${QC_status_contaminations} == "nie" ]]; then
    touch VFDB_summary_dummy.txt; touch VFDB_summary_Escherichia.txt  ; touch VFDB_summary_Shigella.txt
    # json na blad QC
  else
    if [[ ${GENUS} == "Salmonella" || ${GENUS} == "Escherichia" || ${GENUS} == "Campylobacter" ]]; then
      SPEC2=""

      if [ "${GENUS}" == "Escherichia" ]; then
        # for Escherichia we must also check Shigella
        SPEC2="Shigella"  
      fi

      PIDENT=80 # minimalna identycznosc sekwencyjna aby stwierdzic ze jest hit 
      COV=80 # minimalne pokrycie query i hitu aby stwierdzic ze jest hit 
      EVAL=0.01  # maksymalne e-value
      /opt/docker/EToKi/externals/run_VFDB.sh $ffn ${task.cpus} ${GENUS} \${PIDENT} \${EVAL} \${COV}

      mv VFDB_summary.txt VFDB_summary_${GENUS}.txt

      if [ \${SPEC2} == "Shigella" ]; then
        /opt/docker/EToKi/externals/run_VFDB.sh $ffn ${task.cpus} \${SPEC2} \${PIDENT} \${EVAL} \${COV}
        mv VFDB_summary.txt VFDB_summary_\${SPEC2}.txt
      fi
    else
      touch VFDB_summary_dummy.txt; touch VFDB_summary_Escherichia.txt  ; touch VFDB_summary_Shigella.txt
      # json na bledny rodzaj
    fi # koniec if-a na rodzja
  
  fi # koniec if-a na przejscie QC
  """

}

process parse_VFDB_ecoli {
  // Parser wynikow dla E.coli w celu okreslenia czy jest to STEC/VTEC itd ...
  tag "Predicting phenotype for sample $x"
  publishDir "pipeline_wyniki/${x}/", mode: 'copy'
  cpus params.cpus
  input:
  tuple val(x), path('VFDB_summary_Escherichia.txt'), path('VFDB_summary_Shigella.txt'), val(SPECIES), val(GENUS), val(QC_status), val(QC_status_contaminations)
  output:
  tuple val(x), path('VFDB_phenotype.txt')
  //when:
  //GENUS == 'Escherichia'
  script:
  """
    if [[ ${QC_status} == "nie"  || ${QC_status_contaminations} == "nie" ]]; then
       touch VFDB_phenotype.txt
       # json na zle QC
    else
      if [ ${GENUS} == 'Escherichia' ]; then 
        touch VFDB_phenotype.txt
        ### STEC ###
        ### Geny STX1 lub 2 
    
        STEC=0
        STEC=`cat VFDB_summary_Escherichia.txt | grep "stx1\\|stx2" | grep -v BRAK | wc -l`
        GENY_STEC=`cat VFDB_summary_Escherichia.txt | grep "stx1\\|stx2" | grep -v BRAK | cut -f4 | tr "\\n" " "`
        if [ \${STEC} -gt 0 ]; then
          echo -e "STEC\\t\${GENY_STEC}" >> VFDB_phenotype.txt
        fi
        ### Konice 

        ### EPEC ###
        ### EPEC definiujemy obecnosci intyminy
        ### ktora wystepuje w wynikach 2 razy jako czesc roznych sciezek (stad head -1 nizej)

        EPEC=0
        # eaeH to nie intymina
        EPEC=`cat VFDB_summary_Escherichia.txt | grep -w "eae"  | grep -v BRAK | wc -l`
        GENY_EPEC=`cat VFDB_summary_Escherichia.txt | grep -w "eae"  | grep -v BRAK | cut -f4 | tr "\\n" " "`
        if [ \${EPEC} -gt 0 ]; then
          echo -e "EPEC\\t\${GENY_EPEC}" >> VFDB_phenotype.txt
        fi
 
        ### Koniec

        ### EAEC ###
        ### Definicja Tomka i VFDB do uzgodnienia
        ### geny adherence-  aafA; aafB; aafC; aafD 
        ### geny wirulencji east1 oraz set1a i set1b
        ### dyspersyna - gen aap

        ### Geny AggR (kontroler kilku genow w tym dyspersyny), i gen aa1c czesc secration system
    
        EAST=0
        SET=0  
        AAF=0 
        AAP=0 

        AGGR=0 
        AAIC=0
 
        EAST=`cat VFDB_summary_Escherichia.txt | grep -w east1 | grep -v BRAK | wc -l`
        SET=`cat VFDB_summary_Escherichia.txt | grep "set1A\\|set1B" | grep -v BRAK | wc -l`
        AAF=`cat VFDB_summary_Escherichia.txt | grep "aafA\\|aafB\\|aafC\\|aafD" | grep -v BRAK | wc -l`
        AAP=`cat VFDB_summary_Escherichia.txt | grep -w aap | grep -v BRAK | wc -l`

        AGGR=`cat VFDB_summary_Escherichia.txt | grep -w aggR | grep -v BRAK | wc -l`
        AAIC=`cat VFDB_summary_Escherichia.txt | grep -w aaic-hcp |  grep -v BRAK | wc -l`
 
        if [ \${AGGR} -gt 0 ] && [ \${AAIC} -gt 0 ] ; then 
          echo -e "EAEC\\taggR\\taaiC" >> VFDB_phenotype.txt
        elif [ \${AGGR} -gt 0 ]; then
          echo -e "EAEC\\taggR" >> VFDB_phenotype.txt
        elif [ \${AAIC} -gt 0 ]; then
          echo -e "EAEC\\taaiC" >> VFDB_phenotype.txt
        fi

        ### Konice

        ### EIEC ###
        ### Tu jest prosta definicaja ale gen jest nie w bazie Ecoli a w bazie Shigella
        IPAH=0
  
        IPAH=`cat VFDB_summary_Shigella.txt | grep ipaH | grep -v BRAK | wc -l`
        IPAH_GENES=`cat VFDB_summary_Shigella.txt | grep ipaH | grep -v BRAK |  cut -f4 | tr "\\n" " "`
        if [ \${IPAH} -gt 0 ]; then
            echo -e "EIEC\\t\${IPAH_GENES}" >> VFDB_phenotype.txt
        fi
   
        ### Koniec 

        ### ETEC ###
        ELT=0 # Heat liable # wiele izoform
        EST=0 # Heat Stable
   
        ELT=`cat VFDB_summary_Escherichia.txt | grep elt | grep -v BRAK | wc -l`
        ELT_GENES=`cat VFDB_summary_Escherichia.txt | grep elt | grep -v BRAK | cut -f4 | tr "\\n" " "`
        EST=`cat VFDB_summary_Escherichia.txt | grep estIa | grep -v BRAK | wc -l`  
   
        if [ \${ELT} -gt 0 ] && [ \${EST} -gt 0 ] ; then
          echo -e "ETEC\\t\${ELT_GENES}\\testIa" >> VFDB_phenotype.txt
        elif [ \${ELT} -gt 0 ]; then
          echo -e "ETEC\\t\${ELT_GENES}" >> VFDB_phenotype.txt
        elif [ \${EST} -gt 0 ]; then
          echo -e "ETEC\testIa" >> VFDB_phenotype.txt
        fi

        ### Koniec  
    else
      touch VFDB_phenotype.txt
      # json na zly gatunek
    fi # koniec if-a na zly gatunek
  fi # Koniec if-a na Quality
  """
}

process run_spifinder {
  // works only for salmonella
  container  = 'salmonella_illumina:2.0'
  tag "Predicting virulence islands for sample $x"
  publishDir "pipeline_wyniki/${x}", mode: 'copy'
  // cpus params.cpus
  // maxForks 5
  input:
  tuple val(x), path(fasta), val(QC_status), val(SPECIES), val(GENUS), val(QC_status_contaminations)
  output:
  tuple val(x), path('spifinder_results/*')
  // when:
  // GENUS == 'Salmonella'  
  script:
  """
  mkdir spifinder_results # program wymaga tworzenia katalogu samodzielnie
  # -mp opcja jak szukamy wysp kma jest gdy inputem sa dane surowe, blastn gdy podajemy zlozony genom/contigi
  # -p sciezka do bazy, sciagana w trakcie budowy kontenera
  # -l i -t to parametrty na alignment coverage i seq id
  # -x to rozszerzony output

  if [[ ${QC_status} == "nie"  || ${QC_status_contaminations} == "nie" ]]; then  
    touch spifinder_results/dummy_file.txt
    # json na zle QC
  else
    if [ ${GENUS} == "Salmonella" ]; then
      python /opt/docker/spifinder/spifinder.py -i $fasta -o spifinder_results -mp blastn -p /opt/docker/spifinder_db/ -l 0.6 -t 0.9 -x
    else
      touch spifinder_results/dummy_file.txt
      # json na zly gatunek
    fi # koniec if-a na zly gatunek
  fi # koniec if-a na zle QC
  """
}



process run_kraken2_illumina {
  // kraken2 instalowany jest przez ETOKI i jest w path kontenera do salmonelli
  tag "Run kraken2 for sample:${x}"
  container  = 'salmonella_illumina:2.0'
  publishDir "pipeline_wyniki/${x}/kraken2", mode: 'copy', pattern: "report_kraken2_individualreads.txt"
  publishDir "pipeline_wyniki/${x}/kraken2", mode: 'copy', pattern: "report_kraken2.txt"
  publishDir "pipeline_wyniki/${x}/kraken2", mode: 'copy', pattern: "Summary_kraken_*.txt"
  containerOptions "--volume ${params.kraken2_db_absolute_path_on_host}:/home/external_databases/kraken2"
  maxForks 6

  input:
  tuple val(x), path(reads), val(QC_STATUS), val(TOTAL_BASES)

  output:
  tuple val(x), path('report_kraken2.txt'), path('report_kraken2_individualreads.txt'), path('Summary_kraken_genera.txt'), path('Summary_kraken_species.txt')

  script:
  """
  if [ ${QC_STATUS} == "nie" ]; then
    # upstream module failed we produce empty files do the pipeline can execute
    touch report_kraken2.txt
    touch report_kraken2_individualreads.txt
    touch Summary_kraken_genera.txt
    touch Summary_kraken_species.txt
  else
    kraken2 --db /home/external_databases/kraken2 \
            --report report_kraken2.txt \
            --threads ${params.cpus} \
            --gzip-compressed \
            --minimum-base-quality ${params.quality} \
            --use-names ${reads[0]} ${reads[1]} >> report_kraken2_individualreads.txt 2>&1
    # parse kraken extract two most abundant FAMILIES
    LEVEL="G" # G - genus, S - species
    SPEC1=`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -1 | tr -s " " | cut -f6 | tr -d "="`
    SPEC2=`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f6 | tr -d "="`
    ILE1==`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}'| grep -w \${LEVEL} | sort -rnk 1 | head -1 | tr -s " " | cut -f1 | tr -d " "`
    ILE2==`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}'| grep -w \${LEVEL} | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f1 | tr -d " "`

    echo -e "${x}\t\${SPEC1}\${ILE1}%\t\${SPEC2}\${ILE2}%" >> Summary_kraken_genera.txt

    LEVEL="S"
    SPEC1=`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -1 | tr -s " " | cut -f6 | tr -d "="`
    SPEC2=`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f6 | tr -d "="`
    ILE1==`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}'| grep -w \${LEVEL} | sort -rnk 1 | head -1 | tr -s " " | cut -f1 | tr -d " "`
    ILE2==`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}'| grep -w \${LEVEL} | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f1 | tr -d " "`

    echo -e "${x}\t\${SPEC1}\${ILE1}%\t\${SPEC2}\${ILE2}%" >> Summary_kraken_species.txt
  fi
  """
}

process run_metaphlan_illumina {
  tag "Run Metaphlan for sample:${x}"
  container  = 'salmonella_illumina:2.0'
  publishDir "pipeline_wyniki/${x}/metaphlan", mode: 'copy', pattern: "report_metaphlan*"
  containerOptions "--volume ${params.metaphlan_db_absolute_path_on_host}:/bowtie_db"
  maxForks 6
  cpus params.cpus
  input:
  tuple val(x), path(reads), val(QC_STATUS), val(TOTAL_BASES)

  output:
  tuple val(x), path('report_metaphlan_SGB.txt'), path('report_metaphlan_species.txt'), path('report_metaphlan_genera.txt')

  script:
  """
  if [ ${QC_STATUS} == "nie" ]; then
    # upstream module failed we produce empty files do the pipeline can execute
    touch report_metaphlan_SGB.txt
    touch report_metaphlan_species.txt
    touch report_metaphlan_genera.txt
  else
    metaphlan ${reads[0]},${reads[1]} --bowtie2out metagenome.bowtie2.bz2 --nproc ${task.cpus} --input_type fastq -o profiled_metagenome.txt --bowtie2db /bowtie_db/ --unclassified_estimation
    # Parsujemy wyniki
    metaphlan metagenome.bowtie2.bz2 --input_type bowtie2out --bowtie2db /bowtie_db/  --nproc ${task.cpus} --tax_lev 't' -o report_metaphlan_SGB.txt
    metaphlan metagenome.bowtie2.bz2 --input_type bowtie2out --bowtie2db /bowtie_db/  --nproc ${task.cpus} --tax_lev 's' -o report_metaphlan_species.txt
    metaphlan metagenome.bowtie2.bz2 --input_type bowtie2out --bowtie2db /bowtie_db/  --nproc ${task.cpus} --tax_lev 'g' -o report_metaphlan_genera.txt
  fi
  """
}

process run_kmerfinder_illumina {
  // This module unlike krakern and metaphlan will pass QC_STATUS and TOTAL_BASES to get_species_illumina
  tag "kmerfinder:${x}"
  container  = 'salmonella_illumina:2.0'
  publishDir "pipeline_wyniki/${x}/kmerfinder", mode: 'copy', pattern: "results.spa"
  containerOptions "--volume ${params.kmerfinder_db_absolute_path_on_host}:/kmerfinder_db"
  maxForks 5
  cpus params.cpus
  input:
  tuple val(x), path(reads), val(QC_STATUS), val(TOTAL_BASES)

  output:
  tuple val(x),path('results.spa'), path('results.txt'), env(SPECIES), env(GENUS), val(QC_STATUS), val(TOTAL_BASES)

  script:
  """
  if [ ${QC_STATUS} == "nie" ]; then
    # upstream module failed we produce empty files do the pipeline can execute
    touch results.spa
    touch results.txt
    SPECIES=""
    GENUS=""
  else
    /opt/docker/kmerfinder/kmerfinder.py -i ${reads[0]} ${reads[1]} -o ./kmerfider_out -db /kmerfinder_db/bacteria/bacteria.ATG -tax /kmerfinder_db/bacteria/bacteria.tax -x -kp /opt/docker/kma/
    cp kmerfider_out/results.spa .
    cp kmerfider_out/results.txt .
  
    SPECIES=`python /data/parse_kmerfinder.py kmerfider_out/data.json species`
    GENUS=`python /data/parse_kmerfinder.py kmerfider_out/data.json genus`
  fi
  """
}




process run_pHierCC_enterobase {
  // Funkcja odpytuje API Enterobase w celu wyciagniecia profuili z tej bazu
  // It work only for Salmonella and Escherichia as Campylo is not present in Enterobase
  maxRetries 3
  errorStrategy 'retry' // in case there is aroblem with internet connection
  container  = 'salmonella_illumina:2.0'
  tag "Predicting hierCC from enterobase for sample $x"
  publishDir "pipeline_wyniki/${x}/pHierCC", mode: 'copy'
  input:
  tuple val(x), path('cgMLST_parsed_output.txt'), path('cgMLST_sample_full_list_of_allels.txt'), path('cgMLST_closest_ST_full_list_of_allels.txt'), val(SPECIES), val(GENUS), val(QC_status), val(QC_status_contaminations)
  output:
  tuple val(x), path('parsed_phiercc_enterobase.txt'), val(SPECIES), val(GENUS), val(QC_status), val(QC_status_contaminations)
  // when:
  // GENUS == 'Salmonella' || GENUS == 'Escherichia'
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

species="$SPECIES"
genus="$GENUS"
qc_status="$QC_status"
qc_status_contaminations="$QC_status_contaminations"

if qc_status == "nie" or qc_status_contaminations == "nie":
    with open('parsed_phiercc_enterobase.txt', 'w') as f1:
        f1.write(f'ST\\tComment\\n')
        f1.write(f'unk\\tUnknown species: {species}\\n')


    # json na zle QC
    sys.exit(0)

time.sleep(np.random.randint(2,20))

API_TOKEN = "${params.enterobase_api_token}"

def __create_request(request_str):
    base64string = base64.b64encode('{0}: '.format(API_TOKEN).encode('utf-8'))
    headers = {"Authorization": "Basic {0}".format(base64string.decode())}
    request = urllib.request.Request(request_str, None, headers)
    return request

def getST(my_file):
    with open(my_file) as f:
        for line in f:
            line = line.rsplit()
            if line[0] == "ST_sample":
                # pass the first line
                continue
            else:
                # the second line in a file is what we want
                return line[0], line[1], int(line[2])


ST_sample, ST_matching, my_dist = getST('cgMLST_parsed_output.txt')

if genus == 'Salmonella':
    DATABASE="senterica" 
    scheme_name="cgMLST_v2"
    phiercc_header='ST\\tHC0\\tHC2\\tHC5\\tHC10\\tHC20\\tHC50\\tHC100\\tHC200\\tHC400\\tHC900(ceBG)\\tHC2000\\tHC2600\\tHC2850(subsp.)\\n'
    lista_kluczy = ['d0', 'd2', 'd5', 'd10', 'd20', 'd50', 'd100', 'd200' , 'd400', 'd900', 'd2000', 'd2600', 'd2850']
elif genus == 'Escherichia':
    DATABASE="ecoli"
    scheme_name="cgMLST"
    phiercc_header='ST\\tHC0\\tHC2\\tHC5\\tHC10\\tHC20\\tHC50\\tHC100\\tHC200\\tHC400\\tHC1100(cgST Cplx)\\tHC1500\\tHC2000\\tHC2350(subsp.)\\n'
    lista_kluczy = ['d0', 'd2', 'd5', 'd10', 'd20', 'd50', 'd100', 'd200' , 'd400', 'd1100', 'd1500', 'd2000', 'd2350']
else:
    with open('parsed_phiercc_enterobase.txt', 'w') as f:
        f.write('Provided genus: {genus} is not part of the Enterobase')
    sys.exit(0)
    # json for wrong species

with open('parsed_phiercc_enterobase.txt', 'w') as f:
    f.write(phiercc_header)
    address = f"https://enterobase.warwick.ac.uk/api/v2.0/{DATABASE}/{scheme_name}/sts?st_id={ST_matching}&scheme={scheme_name}&limit=5"
    try:
        response = urlopen(__create_request(address))
        data = json.load(response)
        lista_poziomow = [data['STs'][0]['info']['hierCC'][x] for x in lista_kluczy] # lista z uporzadkowanymi poziomami
        # modify outpuy to include distance > 0 that mean we have "local" STs
        try:
            last_index = np.where(list(map(lambda x: int(re.findall('\\d+', x)[0]) < my_dist, lista_kluczy)))[0][-1]
            lista_poziomow[:(last_index+1)] = [ST_sample] * (last_index + 1)
        except IndexError:
            pass
        formatted_string = "\t".join(list(map(str, lista_poziomow)))
        f.write(f'{ST_sample}\\t{formatted_string}\\n')
    except HTTPError as Response_error:
        print(f"{Response_error.code} {Response_error.reason}. URL: {Response_error.geturl()}\\n Reason: {Response_error.read()}")
        sys.exit(1)

"""
}

process run_pHierCC_pubmlst {
  container  = 'salmonella_illumina:2.0'
  containerOptions "--volume ${params.db_absolute_path_on_host}:/db"
  tag "Predicting hierCC with local database for sample $x"
  publishDir "pipeline_wyniki/${x}/pHierCC", mode: 'copy'
  input:
  tuple val(x), path('cgMLST_parsed_output.txt'), path('cgMLST_sample_full_list_of_allels.txt'), path('cgMLST_closest_ST_full_list_of_allels.txt'), val(SPECIES), val(GENUS),  val(QC_status), val(QC_status_contaminations)
  output:
  tuple val(x), path('parsed_phiercc_pubmlst.txt'), val(SPECIES), val(GENUS),  val(QC_status), val(QC_status_contaminations)
  //when:
  //SPECIES == 'jejuni'
  script:
"""
#!/usr/bin/python
import  sys
import gzip
import re
import numpy as np

species="$SPECIES"
genus="$GENUS"
qc_status="$QC_status"
qc_status_contaminations="$QC_status_contaminations"

if qc_status == "nie" or qc_status_contaminations == "nie":
    with open('parsed_phiercc_pubmlst.txt', 'w') as f1:
        f1.write(f'ST\\tComment\\n')
        f1.write(f'unk\\tUnknown species: {species}\\n')


    # json na zle QC
    sys.exit(0)


def getST(my_file):
    with open(my_file) as f:
        for line in f:
            line = line.rsplit()
            if line[0] == "ST_sample":
                # pass the first line
                continue
            else:
                # the second line in a file is what we want
                return line[0], line[1], int(line[2])


if species != "jejuni":
    with open('parsed_phiercc_pubmlst.txt', 'w') as f1:
        f1.write(f'ST\\tComment\\n')
        f1.write(f'unk\\tUnknown species: {species}\\n')


    # json na zly gatunek
    sys.exit(0)

ST_sample, ST_matching, my_dist = getST('cgMLST_parsed_output.txt')

phiercc_header='ST\\tHC5\\tHC10\\tHC25\\tHC50\\tHC100\\tHC200\\n'
lista_kluczy = ['d5', 'd10', 'd25', 'd50','d100', 'd200']
directory=f'/db/{genus}/jejuni/pubmlst'


STs_data = np.load(f'{directory}/sts_table.npy', allow_pickle=True).item() 

# Matching ST must be in the data

lista_poziomow = [STs_data[ST_matching][level] for level in lista_kluczy]
try:
    last_index = np.where(list(map(lambda x: int(re.findall('\\d+', x)[0]) < my_dist, lista_kluczy)))[0][-1]
    lista_poziomow[:(last_index + 1)] = [ST_sample] * (last_index +1)
except IndexError:
    pass

with open('parsed_phiercc_pubmlst.txt', 'w') as f:
    f.write(phiercc_header)
    formatted_string = "\t".join(list(map(str, lista_poziomow)))
    f.write(f'{ST_sample}\\t{formatted_string}\\n')

"""
}

process run_pHierCC_local {
  
  // Ponadto Funkcja odpytuje dwa pliki przygotowane przeze mnie 
  // Jeden zawierajacy klastrowanie SINGLE linkage zbudowane na 430k profili z enterobase
  // Drugi zawierajacy klastrowanie COMPLETE linkage zbudowane na 430k profili z enterobae
  container  = 'salmonella_illumina:2.0'
  containerOptions "--volume ${params.db_absolute_path_on_host}:/db"
  tag "Predicting hierCC with local database for sample $x"
  publishDir "pipeline_wyniki/${x}/pHierCC", mode: 'copy'
  input:
  tuple val(x), path('cgMLST_parsed_output.txt'), path('cgMLST_sample_full_list_of_allels.txt'), path('cgMLST_closest_ST_full_list_of_allels.txt'), val(SPECIES), val(GENUS),  val(QC_status), val(QC_status_contaminations)
  output:
  tuple val(x), path('parsed_phiercc_minimum_spanning_tree.txt'), path('parsed_phiercc_maximum_spanning_tree.txt')
  // when:
  // GENUS == 'Salmonella' || GENUS == 'Escherichia' || SPECIES == 'jejuni'
  script:
"""
#!/usr/bin/python
import  sys
import gzip
import re
import numpy as np

species="$SPECIES"
genus="$GENUS"
qc_status="$QC_status"
qc_status_contaminations="$QC_status_contaminations"

if qc_status == "nie" or qc_status_contaminations == "nie":
    with open('parsed_phiercc_minimum_spanning_tree.txt', 'w') as f1, open('parsed_phiercc_maximum_spanning_tree.txt', 'w') as f2:
        f1.write(f'ST\\tComment\\n')
        f1.write(f'unk\\tUnknown species: {species}\\n')

        f2.write(f'ST\\tComment\\n')
        f2.write(f'unk\\tUnknownn species: {species}')

    # json na zle QC
    sys.exit(0)


def getST(my_file):
    with open(my_file) as f:
        for line in f:
            line = line.rsplit()
            if line[0] == "ST_sample":
                # pass the first line
                continue
            else:
                # the second line in a file is what we want
                return line[0], line[1], int(line[2])

ST_sample, ST_matching, my_dist = getST('cgMLST_parsed_output.txt')

if genus ==  'Salmonella':
    directory=f'/db/{genus}/pHierCC_local'
    phiercc_header='ST\\tHC0\\tHC2\\tHC5\\tHC10\\tHC20\\tHC50\\tHC100\\tHC200\\tHC400\\tHC900(ceBG)\\tHC2000\\tHC2600\\tHC2850(subsp.)\\n'
    lista_kluczy = ['d0', 'd2', 'd5', 'd10', 'd20', 'd50', 'd100', 'd200' , 'd400', 'd900', 'd2000', 'd2600', 'd2850']
elif genus == 'Escherichia':
    directory=f'/db/{genus}/pHierCC_local'
    phiercc_header='ST\\tHC0\\tHC2\\tHC5\\tHC10\\tHC20\\tHC50\\tHC100\\tHC200\\tHC400\\tHC1100(cgST Cplx)\\tHC1500\\tHC2000\\tHC2350(subsp.)\\n'
    lista_kluczy = ['d0', 'd2', 'd5', 'd10', 'd20', 'd50', 'd100', 'd200' , 'd400', 'd1100', 'd1500', 'd2000', 'd2350']
elif species == 'jejuni':
    # Provisal selecetion of jejuni
    phiercc_header='ST\\tHC5\\tHC10\\tHC25\\tHC50\\tHC100\\tHC200\\n' 
    lista_kluczy = ['d5', 'd10', 'd25', 'd50', 'd100', 'd200']
    directory=f'/db/{genus}/jejuni/pHierCC_local'

else:
    with open('parsed_phiercc_minimum_spanning_tree.txt', 'w') as f1, open('parsed_phiercc_maximum_spanning_tree.txt', 'w') as f2:
        f1.write('Provided species: {species} is not part of any cgMLST scheme')
        f2.write('Provided species: {species} is not part of any cgMLST scheme')
        # json na zly gatunek
        sys.exit(0)

# 1. Szukanie w wynikach mojego klastrowania z uzyciem single linkage
with open('parsed_phiercc_minimum_spanning_tree.txt', 'w') as f, gzip.open(f'{directory}/profile_single_linkage.HierCC.gz') as f2, open(f'{directory}/profile_single_linkage.HierCC.index') as f3:
    f.write(phiercc_header)
    pointer = 0
    for line in f3:
        line = line.rsplit()
        try:
            if int(ST_matching) < int(line[0]):
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
        if line[0] == ST_matching:
            if re.findall('Salmo', species):    
                lista_poziomow = [line[1], line[3], line[6], line[11], line[21], line[51], line[101], line[201], line[401], line[901], line[2001], line[2601], line[2851]]
            elif re.findall('Escher', species):
                lista_poziomow = [line[1], line[3], line[6], line[11], line[21], line[51], line[101], line[201], line[401], line[1101], line[1501], line[2001], line[2351]]
            elif re.findall('jejun', species):
                lista_poziomow = [line[6], line[11], line[26], line[51], line[101], line[201]]
            else:
                pass

            try:
                last_index = np.where(list(map(lambda x: int(re.findall('\\d+', x)[0]) < my_dist, lista_kluczy)))[0][-1]
                lista_poziomow[:(last_index + 1)] = [ST_sample] * (last_index +1)
            except IndexError:
                pass
            formatted_string = "\t".join(list(map(str, lista_poziomow)))
            f.write(f'{ST_sample}\\t{formatted_string}\\n')
            # nie ma potrzeby dalszego ogladania pliku
            break

# 2. Szukanie pHierCC w wynikach  maximum spanning tree
with open('parsed_phiercc_maximum_spanning_tree.txt', 'w') as f, gzip.open(f'{directory}/profile_complete_linkage.HierCC.gz') as f2, open(f'{directory}/profile_complete_linkage.HierCC.index') as f3:
    f.write(phiercc_header)
    pointer = 0
    for line in f3:
        line = line.rsplit()
        try:
            if int(ST_matching) < int(line[0]):
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
        if line[0] == ST_matching:
            if re.findall('Salmo', species):
                lista_poziomow = [line[1], line[3], line[6], line[11], line[21], line[51], line[101],  line[201], line[401], line[901], line[2001], line[2601], line[2851]]
            elif re.findall('Escher', species):
                lista_poziomow = [line[1], line[3], line[6], line[11], line[21], line[51], line[101],  line[201], line[401], line[1101], line[1501], line[2001], line[2351]]
            elif re.findall('jejun', species):
                lista_poziomow = [line[6], line[11], line[26], line[51], line[101], line[201]]
            else:
                pass
            try:
                last_index = np.where(list(map(lambda x: int(re.findall('\\d+', x)[0]) < my_dist, lista_kluczy)))[0][-1] 
                lista_poziomow[:(last_index + 1)] = [ST_sample] * (last_index +1)
            except IndexError:
                pass 
            formatted_string = "\t".join(list(map(str, lista_poziomow)))
            f.write(f'{ST_sample}\\t{formatted_string}\\n')
            # nie ma potrzeby dalszego ogladania pliku
            break
"""
}


process extract_historical_data_enterobase {
  container  = 'salmonella_illumina:2.0'
  tag "Extracting historical data for sample $x"
  containerOptions "--volume ${params.db_absolute_path_on_host}:/db"
  publishDir "pipeline_wyniki/${x}/", mode: 'copy'
  input:
  tuple val(x), path('parsed_phiercc_enterobase.txt'), val(SPECIES), val(GENUS),  val(QC_status), val(QC_status_contaminations)
  output:
  tuple val(x), path('enterobase_historical_data.txt'), val(SPECIES), val(GENUS), val(QC_status), val(QC_status_contaminations)
  // when:
  // GENUS == 'Salmonella' || GENUS == 'Escherichia'
  script:
"""
#!/usr/bin/python
import numpy as np
import sys
import re 
species="${SPECIES}"
genus="$GENUS"

qc_status="$QC_status"
qc_status_contaminations="$QC_status_contaminations"

if qc_status == "nie" or qc_status_contaminations == "nie":
    with open('enterobase_historical_data.txt', 'w') as f1:
        f1.write(f'ST\\tComment\\n')
        f1.write(f'unk\\tUnknown species: {species}\\n')

    # json na zle QC
    sys.exit(0)

if genus != 'Salmonella' and genus != 'Escherichia':
    with open('enterobase_historical_data.txt', 'w') as f1:
        f1.write(f'ST\\tComment\\n')
        f1.write(f'unk\\tUnknown species: {species}\\n')

    # json na zle QC 
    sys.exit(0)

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

# we keep the data as a dict with "level" i.e. 'HC2' as a key and cluster id i.e. 40 as value
slownik_hiercc = get_hiercc_level('parsed_phiercc_enterobase.txt')

# As a default we look for all the strains that belong to the same cluster at THIS level
# BE AWARE THAT SOMETIMES KEYS HAVE STRANGE NOTATION LIKE "HC1100(cgST Cplx)", THIS WANT WORK IN THAT CASE
try: 
    if re.findall('Salmo', species): 
        phiercc_level_userdefined = '5'
        scheme_name="cgMLST_v2"
    elif re.findall('Escher', species):
        phiercc_level_userdefined = '20' # Escherichia seems to be much more diverse
        scheme_name="cgMLST"
    common_STs = slownik_hiercc[f'HC{phiercc_level_userdefined}']
except KeyError:
    print('Provided key does not exist')
    sys.exit(0)

directory=f'/db/{genus}/Enterobase'
# Now we look for ALL STs belonging to the same cluster as our cluster given phiercc_level
# for that we query our local enterobase instance
expected_ST = {}
STs = np.load(f'{directory}/sts_table.npy', allow_pickle=True)
STs = STs.item()

for klucz,wartosc in STs.items(): 
    if wartosc[f'd{phiercc_level_userdefined}'] == common_STs:
        expected_ST[klucz] = ''
    		

straindata = np.load(f'{directory}/straindata_table.npy',  allow_pickle=True)
straindata = straindata.item()

dane_historyczne = {} # kluczem jest nazwa szczepu , wartoscia 2 elementowa lista z krajem i rokiem
for strain, wartosc in straindata.items():
    for scheme in wartosc['sts']:
        if scheme_name in scheme.values() and str(scheme['st_id']) in expected_ST.keys():
            dane_historyczne[strain] = [wartosc['country'], wartosc['collection_year'], scheme['st_id']]  

# zapisujemy dane
with open('enterobase_historical_data.txt', 'w') as f:
    f.write(f'Strain_id\\tCountry\\tYear\\tST\\n')
    for klucz, wartosc in dane_historyczne.items():
        if wartosc[0] == 'United States' or wartosc[0] == 'USA':
            f.write(f'{klucz}\\tUnited States of America\\t{wartosc[1]}\\t{wartosc[2]}\\n')
        else:
            f.write(f'{klucz}\\t{wartosc[0]}\\t{wartosc[1]}\\t{wartosc[2]}\\n')
"""
}

process plot_historical_data_enterobase {
  container  = 'salmonella_illumina:2.0'
  tag "Plot historical data for sample $x"
  publishDir "pipeline_wyniki/${x}/", mode: 'copy'
  input:
  tuple val(x), path('enterobase_historical_data.txt'), val(SPECIES), val(GENUS), val(QC_status), val(QC_status_contaminations)
  output:
  tuple val(x), path('enterobase_historical_data.html')
  // when:
  // GENUS == 'Salmonella' || GENUS == 'Escherichia'
  script:
// The script requires a geojeson file that is a part of our container 
// The 3 parameters are input file, output prefix, year fromwhich plot the data, data before that year are ignored (to save html size)
// Same options are used for pubmlst version of that script
"""
if [[ ${QC_status} == "nie"  || ${QC_status_contaminations} == "nie" ]]; then
    # Tworzenie json i output
    touch enterobase_historical_data.html
    
else
  if [[ ${GENUS} == "Salmonella" || ${GENUS} == "Escherichia" ]]; then
    python /data/plot_historical_data_plotly.py enterobase_historical_data.txt enterobase_historical_data 2009
  else
    touch enterobase_historical_data.html
    # json na zly gatunek
  fi
fi
"""

}


process extract_historical_data_pubmlst {
  // The script is nearly identical to extract_historical_data_enterobase
  // However it is eqecuted on a different "branch" of a pipeline
  // and must be duplicated 
  container  = 'salmonella_illumina:2.0'
  tag "Extracting historical data for sample $x"
  containerOptions "--volume ${params.db_absolute_path_on_host}:/db"
  publishDir "pipeline_wyniki/${x}/", mode: 'copy'
  input:
  tuple val(x), path('parsed_phiercc_pubmlst.txt'), val(SPECIES), val(GENUS), val(QC_status), val(QC_status_contaminations)
  output:
  tuple val(x), path('pubmlst_historical_data.txt'), val(SPECIES), val(GENUS), val(QC_status), val(QC_status_contaminations)
  // when:
  // SPECIES == 'jejuni'

  script:
"""
#!/usr/bin/python
import numpy as np
import sys
import re
species="${SPECIES}"

qc_status="$QC_status"
qc_status_contaminations="$QC_status_contaminations"

if qc_status == "nie" or qc_status_contaminations == "nie":
    with open('pubmlst_historical_data.txt', 'w') as f1:
        f1.write(f'ST\\tComment\\n')
        f1.write(f'unk\\tUnknown species: {species}\\n')

    # json na zle QC
    sys.exit(0)

if species != 'jejuni':
    with open('pubmlst_historical_data.txt', 'w') as f1:
        f1.write(f'ST\\tComment\\n')
        f1.write(f'unk\\tUnknown species: {species}\\n')

    # json na zle QC
    sys.exit(0)



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



directory='/db/Campylobacter/jejuni/pubmlst' # where are the data for this species 
scheme_name='C. jejuni / C. coli cgMLST v2' # name of the scheme in that database 

slownik_hiercc = get_hiercc_level('parsed_phiercc_pubmlst.txt') # data for our sample 

phiercc_level_userdefined = '25' # pubmlst keeps only values of 5, 10, 25, 50, 100 and 200 for Campylo

common_STs = slownik_hiercc[f'HC{phiercc_level_userdefined}'] 


straindata = np.load(f'{directory}/straindata_table.npy',  allow_pickle=True)
straindata = straindata.item()

dane_historyczne = {} # strain id is a key, as a value we keep country and date
for strain, wartosc in straindata.items():
    if str(wartosc['hiercc'][f'd{phiercc_level_userdefined}']) == str(common_STs):
        for scheme in  wartosc['sts']:
            if scheme['scheme_name'] == scheme_name:
                dane_historyczne[strain] = [wartosc['country'], wartosc['year'], scheme['st_id']]

# zapisujemy dane
with open('pubmlst_historical_data.txt', 'w') as f:
    f.write(f'Strain_id\\tCountry\\tYear\\tST\\n')
    for klucz, wartosc in dane_historyczne.items():
        if wartosc[0] == 'United States' or wartosc[0] == 'USA':
            f.write(f'{klucz}\\tUnited States of America\\t{wartosc[1]}\\t{wartosc[2]}\\n')
        elif 'UK' in wartosc[0]:
            f.write(f'{klucz}\\tUnited Kingdom\\t{wartosc[1]}\\t{wartosc[2]}\\n')
        else:
            f.write(f'{klucz}\\t{wartosc[0]}\\t{wartosc[1]}\\t{wartosc[2]}\\n')
"""
}


process plot_historical_data_pubmlst {
  container  = 'salmonella_illumina:2.0'
  tag "Plot historical data for sample $x"
  publishDir "pipeline_wyniki/${x}/", mode: 'copy'
  input:
  tuple val(x), path('pubmlst_historical_data.txt'), val(SPECIES), val(GENUS), val(QC_status), val(QC_status_contaminations)
  output:
  tuple val(x), path('*html')
  // when:
  // SPECIES == 'jejuni'
  script:
// The script requires a geojeson file that is a part of our container
"""

if [[ ${QC_status} == "nie"  || ${QC_status_contaminations} == "nie" ]]; then
    # Tworzenie json i output
    touch empty.html

else
  if [ ${SPECIES} == "jejuni" ]; then
    python /data/plot_historical_data_plotly.py pubmlst_historical_data.txt pubmlst_historical_data 2009
  else
    touch empty.html
    # json na zly gatunek
  fi
fi
"""

}

process run_amrfinder {
  container  = 'salmonella_illumina:2.0'
  tag "Predicting microbial resistance with AMRfinder for sample $x"
  publishDir "pipeline_wyniki/${x}/AMRplus_fider", mode: 'copy', pattern: "AMRfinder*"
  containerOptions "--volume ${params.AMRFINDER_db_absolute_path_on_host}:/AMRfider"
  input:
  tuple val(x), path(fasta), val(QC_status), val(SPECIES), val(GENUS), val(QC_status_contaminations)
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
  // when:
  // GENUS == 'Salmonella' || GENUS == 'Escherichia' || GENUS == 'Campylobacter'
  script:
  """
  if [[ ${QC_status} == "nie"  || ${QC_status_contaminations} == "nie" ]]; then
    touch AMRfinder_resistance.txt
    touch AMRfinder_virulence.txt
    # json z wynikami
    # Komentarz NIE uzywac exit 0 wewnatrz script
  else
    if [[ ${GENUS} == "Salmonella" || ${GENUS} == "Escherichia" || ${GENUS} == "Campylobacter" ]]; then
      amrfinder --blast_bin /blast/bin -n $fasta -d /AMRfider  -i 0.9 -c 0.5 -o initial_output.txt -O ${GENUS} --plus 
 
      cat initial_output.txt  | awk 'BEGIN{FS="\\t"}; {if(\$9 == "AMR" || \$1 == "Protein identifier") print \$0}' > AMRfinder_resistance.txt
      cat initial_output.txt  | awk 'BEGIN{FS="\\t"}; {if(\$9 == "VIRULENCE" || \$1 == "Protein identifier") print \$0}' > AMRfinder_virulence.txt
    else
      touch AMRfinder_resistance.txt
      touch AMRfinder_virulence.txt
      # json z wynikami
    fi
  fi
  """ 
}

//process build_model {
  // Ogolny proces do budowania modelu
  // Zostaje do implementacji przez Michala K.
//}


process run_plasmidfinder {
  container  = 'salmonella_illumina:2.0'
  tag "Predicting plasmids for sample $x"
  publishDir "pipeline_wyniki/${x}/", mode: 'copy'
  input:
  tuple val(x), path(fasta), val(QC_status), val(SPECIES), val(GENUS), val(QC_status_contaminations)
  output:
  tuple val(x), path('plasmidfinder/results_tab.tsv')
  // when:
  // GENUS == 'Salmonella' || GENUS == 'Escherichia' || GENUS == 'Campylobacter'
  script:
  """
  # -i to oczywiscie input na podstawie jego rozszerzenia program wybiera metode do analizy (kma dla fastq i blastn dla fasta)
  # -o to katalog z wynikami, musi istniec
  # -p to sciezka do katalog z bazami (tymi z repo plasmidfinder_db)
  # -l to minimalny procent sekwencji jaki musi alignowac sie na genom
  # -t to minimalne sequence identity miedzy query a subject
  # -x printuj dodatkowe dane w output (alignmenty)
  
  if [[ ${QC_status} == "nie"  || ${QC_status_contaminations} == "nie" ]]; then
    mkdir plasmidfinder; touch plasmidfinder/results_tab.tsv
    # json zle QC
  else  
    if [[ ${GENUS} == "Salmonella" || ${GENUS} == "Escherichia" || ${GENUS} == "Campylobacter" ]]; then
      mkdir plasmidfinder # program wymaga tworzenia katalogu samodzielnie
      /opt/docker/plasmidfinder/plasmidfinder.py  -i $fasta -o plasmidfinder -p /opt/docker/plasmidfinder_db -l 0.6 -t 0.9 -x
    else
      mkdir plasmidfinder; touch plasmidfinder/results_tab.tsv
     # json zly gatunek
    fi
  fi
  """
}

process run_virulencefinder {
  container  = 'salmonella_illumina:2.0'
  tag "Predicting plasmids for sample $x"
  publishDir "pipeline_wyniki/${x}/virulencefinder", mode: 'copy', pattern: "results_tab.tsv"
  input:
  tuple val(x), path(fasta),  val(QC_status), val(SPECIES), val(GENUS), val(QC_status_contaminations)
  output:
  tuple val(x), path('results_tab.tsv')
  // when:
  // GENUS == 'Escherichia'
  script:
  """
  # -i to oczywiscie input na podstawie jego rozszerzenia program wybiera metode do analizy (kma dla fastq i blastn dla fasta)
  # -o to katalog z wynikami, musi istniec
  # -p to sciezka do katalog z bazami (katalog z repo virulencefinder_db)
  # -l to minimalny procent sekwencji jaki musi alignowac sie na genom
  # -t to minimalne sequence identity miedzy query a subject
  # -d to nazwa bazy/organizmu
  # -x printuj dodatkowe dane w output (alignmenty)
  # nie wiem czy to kwestia niezgodnosci mojego virulencefindera z bazami
  # ale program printuje mase output (choc tworzy poprawne pliki w koncu to parser blasta)
  mkdir virulencefinder # program wymaga tworzenia katalogu samodzielnie

  if [[ ${QC_status} == "nie"  || ${QC_status_contaminations} == "nie" ]]; then
    touch results_tab.tsv
    # json na zle QC
  else
    if [ ${GENUS} == "Escherichia" ]; then
      /opt/docker/virulencefinder/virulencefinder.py  -i $fasta -o virulencefinder -p /opt/docker/virulencefinder_db  -d virulence_ecoli -l 0.6 -t 0.9 -x >> log 2>&1
      cp virulencefinder/results_tab.tsv .
    else
      touch results_tab.tsv
      #json na zly gatunek
    fi # koniec if-a na zly gatunek
  fi # koniec if-a na zle QC
  """
}

// process parse_virulencefinder_ecoli {
// Po ustaleniu listy genow
// }




process get_species_illumina {
// Process laczy ouputy predykcji z krakena2, metaphlan i kmerfindera
container  = 'salmonella_illumina:2.0'
tag "Predicting species for ${x}"
publishDir "pipeline_wyniki/${x}", mode: 'copy', pattern: "predicted_genus_and_species.txt"
publishDir "pipeline_wyniki/${x}/json_output", mode: 'copy', pattern: "contaminations.json"
input:
tuple val(x), path('report_kraken2.txt'), path('report_kraken2_individualreads.txt'), path('Summary_kraken_genera.txt'), path('Summary_kraken_species.txt'), path('report_metaphlan_SGB.txt'), path('report_metaphlan_species.txt'), path('report_metaphlan_genera.txt'),  path('results.spa'), path('results.txt'), val(KMERFINDER_SPECIES), val(KMERFINDER_GENUS), val(QC_STATUS), val(TOTAL_BASES)

output:
tuple val(x), env(FINALE_SPECIES), env(FINAL_GENUS), env(QC_status_contaminations), emit: species_and_qcstatus
path('predicted_genus_and_species.txt'), emit: to_pubdir 
tuple val(x), path('contaminations.json'), emit: json
tuple val(x), env(QC_status_contaminations), emit: qcstatus_only
script:
"""
if [ ${QC_STATUS} == "nie" ]; then
  QC_status_contaminations="nie"
  FINALE_SPECIES="unknown"
  FINAL_GENUS="unknown"
  echo "This module was eneterd with failed QC and poduced no valid output" >> predicted_genus_and_species.txt
  # contaminations json
  python /opt/docker/EToKi/externals/json_output_contaminations.py bacterial_illumina contaminations.json nie
else 
  QC_status_contaminations="tak"
  PRE_FINALE_SPECIES=""
  cat report_kraken2.txt | grep -w "S" | sort -rnk1 | head -1 | awk '{print \$6,\$7}' >> intermediate.txt
  cat report_metaphlan_species.txt  | grep -v "#" | sort -rnk3 | head -1 | awk '{print \$1}' | sed s'/s__//'g | sed s'/_/ /'g >> intermediate.txt
  echo ${KMERFINDER_SPECIES} | cut -d ' ' -f1-2 >> intermediate.txt

  KRAKEN_GENUS_LEVEL=`cat report_kraken2.txt | grep -w "S" | sort -rnk1 | head -1 | awk '{print int(\$1)}'`
  METAPHLAN_GENUS_LEVEL=`cat report_metaphlan_species.txt  | grep -v "#" | sort -rnk3 | head -1 | awk '{print int(\$3)}'`
  KMERFINDER_COVERAGE=`cat results.txt | head -2 |  tail -1 | cut -f9 | awk '{print int (\$1)}'`

  if [[ \${KRAKEN_GENUS_LEVEL} -lt 50 && \${METAPHLAN_GENUS_LEVEL} -lt 50 && \${KMERFINDER_COVERAGE} -lt 30 ]]; then
    # kraken2 i metaphlan zwracaja ponziej 50% odczytow nalezacych do glownego gatunku
    # kmerfinder zwraca pokrycie pierwszego gatunku ponad 30 
    QC_status_contaminations="nie"
    FINALE_SPECIES="unknown"
    FINAL_GENUS="unknown"
    echo -e "The sample is contaminated or lacks sufficient number of reads" >> predicted_genus_and_species.txt
    python /opt/docker/EToKi/externals/json_output_contaminations.py bacterial_illumina contaminations.json blad
  else
    PRE_FINALE_SPECIES=`cat intermediate.txt | sort | uniq -c | tr -s " " | sort -rnk1 | head -1 | cut -d " " -f3,4`

    if [[ "\${PRE_FINALE_SPECIES}" == *"Salmonel"* ]]; then
      FINALE_SPECIES="\${PRE_FINALE_SPECIES}"
      GENOME_SIZE=5400000
    elif [[ "\${PRE_FINALE_SPECIES}" == *"Escher"* || "\${PRE_FINALE_SPECIES}" == *"Shigella"* ]]; then
      FINALE_SPECIES="Escherichia coli"
      GENOME_SIZE=4800000
    elif [[ "\${PRE_FINALE_SPECIES}" == "Campylobacter coli" ]]; then
      FINALE_SPECIES="jejuni"
      GENOME_SIZE=1800000
    elif [[ "\${PRE_FINALE_SPECIES}" == "Campylobacter jejuni" ]]; then
      FINALE_SPECIES="jejuni"
      GENOME_SIZE=1800000
    elif [[ "\${PRE_FINALE_SPECIES}" == "Campylobacter concisus" ]]; then
      FINALE_SPECIES="concisus"
      GENOME_SIZE=1800000
    elif [[ "\${PRE_FINALE_SPECIES}" == "Campylobacter curvus" ]]; then
      FINALE_SPECIES="concisus"
      GENOME_SIZE=1800000
    elif [[ "\${PRE_FINALE_SPECIES}" == "Campylobacter fetus" ]]; then
      FINALE_SPECIES="fetus"
      GENOME_SIZE=1800000
    elif [[ "\${PRE_FINALE_SPECIES}" == "Campylobacter helveticus" ]]; then
      FINALE_SPECIES="helveticus"
      GENOME_SIZE=1800000
    elif [[ "\${PRE_FINALE_SPECIES}" == "Campylobacter hyointestinalis" ]]; then
      FINALE_SPECIES="hyointestinalis"
      GENOME_SIZE=1800000
    elif [[ "\${PRE_FINALE_SPECIES}" == "Campylobacter insulaenigrae" ]]; then
      FINALE_SPECIES="insulaenigrae"
      GENOME_SIZE=1800000
    elif [[ "\${PRE_FINALE_SPECIES}" == "Campylobacter lanienae" ]]; then
      FINALE_SPECIES="lanienae"
      GENOME_SIZE=1800000
    elif [[ "\${PRE_FINALE_SPECIES}" == "Campylobacter lari" ]]; then
      FINALE_SPECIES="lari"
      GENOME_SIZE=1800000
    elif [[ "\${PRE_FINALE_SPECIES}" == "Campylobacter sputorum" ]]; then
      FINALE_SPECIES="sputorum"
      GENOME_SIZE=1800000
    elif [[ "\${PRE_FINALE_SPECIES}" == "Campylobacter upsaliensis" ]]; then
      FINALE_SPECIES="upsaliensis"
      GENOME_SIZE=1800000
    else
      FINALE_SPECIES="\${PRE_FINALE_SPECIES}"
      GENOME_SIZE=6000000
    fi

    cat report_kraken2.txt | grep -w "G" | sort -rnk1 | head -1 | awk '{print \$6,\$7}'  >> intermediate_genus.txt
    cat report_metaphlan_genera.txt  | grep -v "#" | sort -rnk3 | head -1 | awk '{print \$1" "}' | sed s'/g__//'g | sed s'/_/ /'g >> intermediate_genus.txt
    echo -e "${KMERFINDER_GENUS} " >> intermediate_genus.txt

    #Ecoli and Schigella are the same thing
    FINAL_GENUS=`cat intermediate_genus.txt | sed s'/Shigella/Escherichia/'g | sort | uniq -c | tr -s " " | sort -rnk1 | head -1 | cut -d " " -f3`
    echo -e "Final genus:\t\${FINAL_GENUS}\nFinal species:\t\${FINALE_SPECIES}" >> predicted_genus_and_species.txt

    # ostani przelacnik liczba zasad to co najmniej 30x dlugosc "oczekiwanego" genomu
    TEORETICAL_COVERAGE=`awk -v g_size="\${GENOME_SIZE}" -v t_bases="${TOTAL_BASES}" 'BEGIN {print t_bases/g_size}'`
    if [ `awk -v tot_cov=\${TEORETICAL_COVERAGE} 'BEGIN {if(tot_cov > 20) {print 1} else {print 0}}'` -eq 1 ]; then
      # Tu jestesmy lagodniejsi, ostateczna wartosc sredniego pokrycia uzyjemy dopiero przy skladaniu genomu
      # aLe przy liczbie zasad mniejszej niz 20 x teoretycznego pokrycia nawet nie ma co probowac
      python /opt/docker/EToKi/externals/json_output_contaminations.py bacterial_illumina contaminations.json tak
      QC_status_contaminations="tak"
      echo " OK w \${TEORETICAL_COVERAGE}" >> tmp.log
   else
     echo " Blad w \${TEORETICAL_COVERAGE}" >> tmp.log
     python /opt/docker/EToKi/externals/json_output_contaminations.py bacterial_illumina contaminations.json blad
     QC_status_contaminations="nie"
   fi
 fi

fi
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
  // publishDir "pipeline_wyniki/${x}", mode: 'copy', pattern: "*"
  cpus params.cpus
  // in --deterministic is awfully slow because in uses 1 cpu for some task
  // hence I increased maxforks to 10
  maxForks 10 
  input:
  tuple val(x), path(fastq_gz), val(SPECIES), val(GENUS), val(QC_status_contaminations)
  output:
  tuple val(x), path('output/assembly.fasta'), env(QC_status)
  
  script:
  """
  # /data/Flye to sciezka z Flye instalowanego z github, uwaga
  # w kontenerze tez jest flye instalowant przez etoki i ten jest w PATH

  ### QC_status_contaminations na razie ignorujemy, zostanie ono przekazane do kanalu jeszcze raz po zaproponowaniu genomu
  ### i tam bedzie uzywane

  QC_status="pass"
  ### Tutaj umieszcamy warunki na QC_status zwiazny ze startowa sekwencja (np za malo odczytow)
  #### Do ustalenia, poki co 1000 czyli skrajnie malo
  NO_READS=`zcat $fastq_gz | grep "^+" | wc -l`
  
  if [ \${NO_READS} -lt 10000 ]; then
    QC_status="nie"
    mkdir output
    echo ">dummy_contig" >> output/assembly.fasta
    echo "AAAAAAAAAAAAA" >> output/assembly.fasta
  else
    if [[ "${GENUS}" == *"Salmo"* ]]; then
        GENOME_SIZE="5.4m"
    elif [[ "${GENUS}" == *"Escher"* ]]; then
        GENOME_SIZE="4.6m"
    elif [ ${GENUS} == "Campylobacter" ]; then
        GENOME_SIZE="1.8m"
    else
        GENOME_SIZE="5m" # trafilem na zly organizm wiec wpisuje 5m ten genom i tak nie bedzie wykorzystany
    fi
    /opt/docker/Flye/bin/flye --nano-raw ${fastq_gz} -g \${GENOME_SIZE} -o output -t ${task.cpus} -i 3 --no-alt-contig --deterministic
  fi


  

  """
 
}

process clean_fastq_nanopore {
  container  = 'salmonella_illumina:2.0'
  tag "Fixing fastq dla sample $x"
  input:
  tuple val(x), path(read)
  output:
  tuple val(x), path('prep_out_L1_SE.fastq.gz')
  script:
  """
  python /opt/docker/EToKi/EToKi.py prepare --se ${read} -c ${params.cpus} -q ${params.quality} -p prep_out 
  """
}



process run_minimap2 {
  // Proces do mapowania odczytow na scaffold
  tag "Remapping of reads to predicted scaffold for sample $x"
  container  = 'salmonella_illumina:2.0' 
  // publishDir "pipeline_wyniki/${x}", mode: 'copy', pattern: "*"
  maxForks 5
  input:
  tuple val(x), path(fasta), val(QC_status), path(reads) 
  output:
  tuple val(x), path('sorted.bam'), path(fasta), val(QC_status)
  script:
  """
  # minmap jest zarowno w PATH z etoki/externals jak i w /data/Flye/bin
  if [ ${QC_status} == "nie" ]; then
    touch sorted.bam
    # dummy output aby skypt poszedl dalej
  else
    minimap2 -a -x map-ont -t ${params.cpus} $fasta $reads | samtools view -bS -F 2052 - | samtools sort -@ ${params.cpus} -o sorted.bam -
  fi
  """
}

process run_minimap2_2nd {
  // Proces do mapowania odczytow na scaffold
  // W nanopre w jednym workflow uzywam go 2 razy wiec musze zrobic ta glupia kopie
  tag "Remapping of reads to predicted scaffold for sample $x"
  container  = 'salmonella_illumina:2.0'
  // publishDir "pipeline_wyniki/${x}", mode: 'copy', pattern: "*"
  maxForks 5
  input:
  tuple val(x), path(fasta), val(QC_status), path(reads)
  output:
  tuple val(x), path('sorted.bam'), path(fasta), val(QC_status)
  script:
  """
  if [ ${QC_status} == "nie" ]; then
    touch sorted.bam
    # dummy output aby skypt poszedl dalej
  else
    # minmap jest zarowno w PATH z etoki/externals jak i w /data/Flye/bin

    minimap2 -a -x map-ont -t ${params.cpus} $fasta $reads | samtools view -bS -F 2052 - | samtools sort -@ ${params.cpus} -o sorted.bam -
  fi
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
  tag "Medaka for sample $x"
  // publishDir "pipeline_wyniki/${x}/medaka", mode: 'copy', pattern: 'latest_pilon.*'
  maxForks 5
  input:
  tuple val(x), path(bam1), path(fasta), val(QC_status)
  output:
  // tuple val(x), path('postmedaka.fasta'), path('medaka_annotated_filtered.vcf.gz'), emit: ALL
  tuple val(x), path('postmedaka.fasta'), val(QC_status), emit: ONLY_GENOME
  script:
  """
  if [ ${QC_status} == "nie" ]; then
    echo ">dummy_contig" >> postmedaka.fasta
    echo "AAAAAAAAAAAAA" >> postmedaka.fasta
  else
    # indeksacja bam-ow
    samtools index $bam1
 
    MODEL="r941_min_hac_g507"  
  
    medaka inference --model \${MODEL} \
                     --threads ${params.cpus} \
                     $bam1 \
                     forvariants.hdf

  
    medaka vcf forvariants.hdf  $fasta medaka.vcf
    medaka tools annotate medaka.vcf $fasta $bam1 medaka_annotated.vcf
    bcftools sort medaka_annotated.vcf >> medaka_annotated_sorted.vcf
    bgzip medaka_annotated_sorted.vcf
    tabix medaka_annotated_sorted.vcf.gz
  
    qual=13
    min_cov=20
  
    bcftools filter -O z -o medaka_annotated_filtered.vcf.gz -i "GQ > \${qual} && DP >= \${min_cov}" medaka_annotated_sorted.vcf.gz
    tabix medaka_annotated_filtered.vcf.gz


    cat $fasta | bcftools consensus medaka_annotated_filtered.vcf.gz >> postmedaka.fasta
  fi
  """
}

process run_kraken2_nanopore {
  // kopia processu do illuminy, ale z uwzglednieniem ze jest tylko jeden plik z odczytami
  tag "kraken2:${x}"
  container  = 'salmonella_illumina:2.0'
  publishDir "pipeline_wyniki/${x}/kraken2", mode: 'copy', pattern: "report_kraken2_individualreads.txt"
  publishDir "pipeline_wyniki/${x}/kraken2", mode: 'copy', pattern: "report_kraken2.txt"
  publishDir "pipeline_wyniki/${x}/kraken2", mode: 'copy', pattern: "Summary_kraken*.txt"
  containerOptions "--volume ${params.kraken2_db_absolute_path_on_host}:/home/external_databases/kraken2"
  maxForks 6

  input:
  tuple val(x), path(reads)

  output:
  tuple val(x), path('report_kraken2.txt'), path('report_kraken2_individualreads.txt'), path('Summary_kraken_genera.txt'), path('Summary_kraken_species.txt')

  script:
  """
  kraken2 --db /home/external_databases/kraken2 \
          --report report_kraken2.txt \
          --threads ${params.cpus} \
          --gzip-compressed \
          --minimum-base-quality ${params.quality} \
          --use-names ${reads} >> report_kraken2_individualreads.txt 2>&1
  
  LEVEL="G"

  SPEC1=`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -1 | tr -s " " | cut -f6 | tr -d "="`
  SPEC2=`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f6 | tr -d "="`
  ILE1==`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}'| grep -w \${LEVEL} | sort -rnk 1 | head -1 | tr -s " " | cut -f1 | tr -d " "`
  ILE2==`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}'| grep -w \${LEVEL} | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f1 | tr -d " "`

  echo -e "${x}\t\${SPEC1}\${ILE1}%\t\${SPEC2}\${ILE2}%" >> Summary_kraken_genera.txt

  LEVEL="S"
  SPEC1=`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -1 | tr -s " " | cut -f6 | tr -d "="`
  SPEC2=`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}' | grep -w \${LEVEL} | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f6 | tr -d "="`
  ILE1==`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}'| grep -w \${LEVEL} | sort -rnk 1 | head -1 | tr -s " " | cut -f1 | tr -d " "`
  ILE2==`cat report_kraken2.txt  | awk '{if(\$1 != 0.00) print \$0}'| grep -w \${LEVEL} | sort -rnk 1 | head -2 | tail -1 | tr -s " " | cut -f1 | tr -d " "`

  echo -e "${x}\t\${SPEC1}\${ILE1}%\t\${SPEC2}\${ILE2}%" >> Summary_kraken_species.txt
  """
}

process run_kmerfinder_nanopore {
  tag "kmerfinder:${x}"
  container  = 'salmonella_illumina:2.0'
  publishDir "pipeline_wyniki/${x}/kmerfinder", mode: 'copy', pattern: "results.spa"
  containerOptions "--volume ${params.kmerfinder_db_absolute_path_on_host}:/kmerfinder_db"
  maxForks 5
  cpus params.cpus
  input:
  tuple val(x), path(reads)

  output:
  tuple val(x),path('results.spa'), path('results.txt'), env(SPECIES), env(GENUS)

  script:
  """
  /opt/docker/kmerfinder/kmerfinder.py -i $reads -o ./kmerfider_out -db /kmerfinder_db/bacteria/bacteria.ATG -tax /kmerfinder_db/bacteria/bacteria.tax -x -kp /opt/docker/kma/
  cp kmerfider_out/results.spa .
  cp kmerfider_out/results.txt .
  
  SPECIES=`python /data/parse_kmerfinder.py kmerfider_out/data.json species`
  GENUS=`python /data/parse_kmerfinder.py kmerfider_out/data.json genus`
  """
}

process get_species_nanopore {
// Process laczy ouputy predykcji z krakena2, metaphlan i kmerfindera
container  = 'salmonella_illumina:2.0'
tag "Predicting species for ${x}"
publishDir "pipeline_wyniki/${x}", mode: 'copy', pattern: "predicted_genus_and_species.txt"
input:
tuple val(x), path('report_kraken2.txt'), path('report_kraken2_individualreads.txt'), path('Summary_kraken_genera.txt'), path('Summary_kraken_species.txt'),  path('results.spa'), path('results.txt'), val(KMERFINDER_SPECIES), val(KMERFINDER_GENUS)

output:
tuple val(x), env(FINALE_SPECIES), env(FINAL_GENUS), env(QC_status_contaminations), emit: species
path('predicted_genus_and_species.txt'), emit: output
script:
"""
QC_status_contaminations="pass"
PRE_FINALE_SPECIES=""
cat report_kraken2.txt | grep -w "S" | sort -rnk1 | head -1 | awk '{print \$6,\$7}' >> intermediate.txt
echo ${KMERFINDER_SPECIES} | cut -d ' ' -f1-2  >> intermediate.txt
PRE_FINALE_SPECIES=`cat intermediate.txt | sort | uniq -c | tr -s " " | sort -rnk1 | head -1 | cut -d " " -f3,4`

KRAKEN_GENUS_LEVEL=`cat report_kraken2.txt | grep -w "S" | sort -rnk1 | head -1 | awk '{print int(\$1)}'`
KMERFINDER_COVERAGE=`cat results.txt | head -2 |  tail -1 | cut -f9 | awk '{print int (\$1)}'`

if [[ \${KRAKEN_GENUS_LEVEL} -lt 50 && \${KMERFINDER_COVERAGE} -lt 30 ]]; then
# kraken2 zwraca mniej nizd 50% odczytow nalezacych do glownego gatunku
# kmerfinder zwraca pokrycie pierwszego gatunku ponizej 30
QC_status_contaminations="nie"
fi


if [[ "\${PRE_FINALE_SPECIES}" == *"Salmonel"* ]]; then
    FINALE_SPECIES="\${PRE_FINALE_SPECIES}"
elif [[ "\${PRE_FINALE_SPECIES}" == *"Escher"* || "\${PRE_FINALE_SPECIES}" == *"Shigella"* ]]; then
    FINALE_SPECIES="Escherichia coli"
elif [[ "\${PRE_FINALE_SPECIES}" == "Campylobacter coli" ]]; then
    FINALE_SPECIES="jejuni"
elif [[ "\${PRE_FINALE_SPECIES}" == "Campylobacter jejuni" ]]; then
    FINALE_SPECIES="jejuni"
elif [[ "\${PRE_FINALE_SPECIES}" == "Campylobacter concisus" ]]; then
    FINALE_SPECIES="concisus"
elif [[ "\${PRE_FINALE_SPECIES}" == "Campylobacter curvus" ]]; then
    FINALE_SPECIES="concisus"
elif [[ "\${PRE_FINALE_SPECIES}" == "Campylobacter fetus" ]]; then
    FINALE_SPECIES="fetus"
elif [[ "\${PRE_FINALE_SPECIES}" == "Campylobacter helveticus" ]]; then
    FINALE_SPECIES="helveticus"
elif [[ "\${PRE_FINALE_SPECIES}" == "Campylobacter hyointestinalis" ]]; then
    FINALE_SPECIES="hyointestinalis"
elif [[ "\${PRE_FINALE_SPECIES}" == "Campylobacter insulaenigrae" ]]; then
    FINALE_SPECIES="insulaenigrae"
elif [[ "\${PRE_FINALE_SPECIES}" == "Campylobacter lanienae" ]]; then
    FINALE_SPECIES="lanienae"
elif [[ "\${PRE_FINALE_SPECIES}" == "Campylobacter lari" ]]; then
    FINALE_SPECIES="lari"
elif [[ "\${PRE_FINALE_SPECIES}" == "Campylobacter sputorum" ]]; then
    FINALE_SPECIES="sputorum"
elif [[ "\${PRE_FINALE_SPECIES}" == "Campylobacter upsaliensis" ]]; then
    FINALE_SPECIES="upsaliensis"
else
    FINALE_SPECIES="\${PRE_FINALE_SPECIES}"
fi


cat report_kraken2.txt | grep -w "G" | sort -rnk1 | head -1 | awk '{print \$6,\$7}' >> intermediate_genus.txt
echo -e "${KMERFINDER_GENUS} " >> intermediate_genus.txt

FINAL_GENUS=`cat intermediate_genus.txt | sort | uniq -c | tr -s " " | sort -rnk1 | head -1 | cut -d " " -f3`

echo -e "Final genus:\t\${FINAL_GENUS}\nFinal species:\t\${FINALE_SPECIES}" >> predicted_genus_and_species.txt

"""
}

// SUB WORKFLOWS NANOPORE //

workflow polishing_with_medaka {
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

final_assembly_filtered = extract_final_contigs(single_bams_and_polished_genome)

emit:
extract_final_contigs.out[0] // fasta with contigs that passed coverage filter
extract_final_contigs.out[1] // fasta with contigs that passed coverage filter and second fasta with refejted contigs

}


// SUB WORKFLOWS ILLUMINA //

workflow pilon_first {
// First round of polishing initial scaffold with pilon
take:
initial_scaffold_inner
processed_fastq_inner_PE
processed_fastq_inner_SE

// initial_scaffold_inner jest kanalem zawierajacym
// // ientyfikator probki
// // plik fasta z genomem

// processed_fastq_inner_PE i SE kazdy jest kanalem zawierajacym
// // PE -> identyfikator probki i odczyty PE
// // SE -> identyfikator probki i odczyty SE
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
// Polishing scaffold obtained after first pilon run
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
// Polishing scaffold obtained after second pilon run
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
// This workflow takes scafold polished with 3 rounds of pilon
// And save only contigs that coverage is above params.min_coverage_ratio of average coverage
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
extract_final_contigs_out = extract_final_contigs(merged_bams_into_onefile_and_scaffold_inner)

emit:
extract_final_contigs.out[0] // fasta with contigs that passed coverage filter
extract_final_contigs.out[1] // fasta with contigs that passed coverage filter and second fasta with refejted contigs
}


workflow predict_species_illumina {
// Workflow puszcza 3 programy (kraken2, metaphlan i kmerfinder)
// jego celem jest zwrocenie informacji o gatunku (tak naprawde wazne tylko dla Campylo)
// tak bym nizej mogl okreslic poprawnie baze
take:
initial_fastq
inital_qc_status // Tu bedzie tylko jeden z emitow kanalu z FASTQC

main:
// join FASTQ with QC_status from FASTQC
initial_fastq_and_status = initial_fastq.join(inital_qc_status, by : 0, remainder : true)
// kraken
kraken2_out = run_kraken2_illumina(initial_fastq_and_status)

// Metaphlan
metaphlan_out = run_metaphlan_illumina(initial_fastq_and_status)

// Kmerfinder
kmerfinder_out = run_kmerfinder_illumina(initial_fastq_and_status)

merge1 = metaphlan_out.join(kmerfinder_out, by : 0, remainder : true)
programs_out = kraken2_out.join(merge1,  by : 0, remainder : true)
final_species = get_species_illumina(programs_out)
emit:
final_species.species_and_qcstatus
final_species.qcstatus_only
}

workflow predict_species_nanopore {
// Workflow puszcza 2 programy (kraken2,  kmerfinder)
// Metaphlan puszcza bowtie2 ktoru jest dla krtokich odczytow
take:
initial_fastq
main:
// kraken
kraken2_out = run_kraken2_nanopore(initial_fastq)

// Kmerfinder
kmerfinder_out = run_kmerfinder_nanopore(initial_fastq)

programs_out = kraken2_out.join(kmerfinder_out,  by : 0, remainder : true)
final_species = get_species_nanopore(programs_out)
emit:
final_species.species

}
// MAIN WORKFLOW //

workflow {

// Get data
if(params.machine == 'Illumina') {

Channel
  .fromFilePairs(params.reads)
  .set {initial_fastq}

// FASTQC
run_fastqc_illumina_out = run_fastqc_illumina(initial_fastq)

// Species prediction
(predict_species_out, predict_species_qc) = predict_species_illumina(initial_fastq, run_fastqc_illumina_out.qcstatus_and_values)


//Inilat MLST

initial_mlst_out = run_initial_mlst_illumina(initial_fastq.join(predict_species_out, by : 0)) // initial_mlst_out  przekazuje zarowno fastq jak i qc_status
// FASTQ trimming
processed_fastq = clean_fastq_illumina(initial_mlst_out)

// Initial scaffold
initial_scaffold = spades(processed_fastq.All_path)

// Polishing scaffold with pilon ( 3 times )
first_polish_run = pilon_first(initial_scaffold, processed_fastq.PE_path, processed_fastq.SE_path)
second_polish_run = pilon_second(first_polish_run, processed_fastq.PE_path, processed_fastq.SE_path)
third_polish_run = pilon_third(second_polish_run, processed_fastq.PE_path, processed_fastq.SE_path)

// Remove contigs with coverage less than 0.1 avarage coverage

(final_assembly, final_assembly_with_reject) = calculate_coverage(third_polish_run, processed_fastq.PE_path, processed_fastq.SE_path)

} else if (params.machine == 'Nanopore') {

// Get Data
Channel
  .fromPath(params.reads)
  .map {it -> tuple(it.getName().split("\\.")[0], it)}
  .set {initial_fastq}

// FASTQC
run_fastqc_nanopore(initial_fastq)

// Contaminations/subspecies prediction
predict_species_out = predict_species_nanopore(initial_fastq)

// FASTQ trimming
processed_fastq = clean_fastq_nanopore(initial_fastq)

// initial scaffold with Flye (with 3 internalrounds of polishing)
initial_scaffold = run_flye(processed_fastq.join(predict_species_out))

// one round of assembly polishing with medaka
(final_assembly, final_assembly_with_reject) = polishing_with_medaka(initial_scaffold, processed_fastq)

}


// Assembly quality
(extract_final_stats_statistics, extract_final_stats_genome) = extract_final_stats(final_assembly_with_reject)
 

// Species prediction with Achtman and core genome shemes
final_assembly_with_species = extract_final_stats_genome.join(predict_species_out, by : 0)
MLST_out = run_7MLST(final_assembly_with_species)
parse_7MLST(MLST_out)


cgMLST_out = run_cgMLST(final_assembly_with_species)
(parse_cgMLST_out, parse_cgMLST_only_output) = parse_cgMLST(cgMLST_out)
run_pHierCC_local(parse_cgMLST_out)

run_pHierCC_enterobase_out = run_pHierCC_enterobase(parse_cgMLST_out) // for Salmo and Escher
run_pHierCC_enterobase_out_pubmlst = run_pHierCC_pubmlst(parse_cgMLST_out) // only for jejuni

extract_historical_data_enterobase_out = extract_historical_data_enterobase(run_pHierCC_enterobase_out)
plot_historical_data_enterobase(extract_historical_data_enterobase_out)


extract_historical_data_pubmlst_out = extract_historical_data_pubmlst(run_pHierCC_enterobase_out_pubmlst)
plot_historical_data_pubmlst(extract_historical_data_pubmlst_out)

// AMR predictions
run_resfinder(final_assembly_with_species)
run_amrfinder(final_assembly_with_species)

// plasmids
run_plasmidfinder(final_assembly_with_species)

// Virulence
prokka_out = run_prokka(final_assembly_with_species)
VFDB_out=run_VFDB(prokka_out)
run_virulencefinder(final_assembly_with_species)

// These three modules have "when" instructions and works only for Escher

run_ectyper(final_assembly_with_species)
parse_VFDB_ecoli(VFDB_out.ecoli)

// These three modules have "when" instructions and works only for Salmonella
run_spifinder(final_assembly_with_species)
run_Seqsero(final_assembly_with_species)
run_sistr(final_assembly_with_species)

} 

