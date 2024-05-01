params.cpus = 24
params.min_coverage_ratio = 0.1 // Etoki odrzuca contigi ktorych srednie pokrycie jest mniejsze niz 0.2 globalnego pokrycia
// my nie jestesmy tak ostrzy dajemy 0.1 co tlumaczymy uzyciem innego alignera
params.species = 's.enterica' // inne mozliwosci to c.jejuni c.coli e.coli . Gatungi do resfindera


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

process clean_fastq {
  // Prosta funkcja zaimplementowa w etoki do czyszczenia plikow fastq, ma na celu rozdzielenie odczytow ktore sa sparowane
  // od tych ktore pary nie maja, trimmowanie odczytow na podstawie jakosci, zmiana nazwy odczytow .
  // Ta funkcja generalnie mimikuje dzialanie trimmomatic-a. Wiec teoretycznie mozna uzyc jego
  
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
  python /opt/docker/EToKi/EToKi.py prepare --pe ${read_1},${read_2} -p prep_out
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
  python /opt/docker/EToKi/externals/spades.py -t ${task.cpus}  --pe-1 1 R1.fastq.gz --pe-2 1 R2.fastq.gz --pe-s 2 SE.fastq.gz -o spades_manual
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
  tag "Lacznie bam-ow dla sample $x"
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
  tag "Pilon dla sample $x"
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

  container  = 'salmonella_illumina:2.0'
  // w kontenerze 2.0 dodalem pilona ale nie chce kasowac 1.0
  tag "Pilon dla sample $x"
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
  java -jar /opt/docker/EToKi/externals/pilon.jar  --genome genomic_fasta.fasta --frags $bam1 --unpaired $bam2 --output pilon_polish --vcf --changes --outdir pilon_bwa
 
  # przygotowanie outputu
  # cat pilon_bwa/pilon_polish.fasta | awk  -v ALA=${x} 'BEGIN{OFS=""}; {if(\$0 ~ />/) print \$0,"_", ALA; else print \$0}' >> last_pilon.fasta
  cat pilon_bwa/pilon_polish.fasta >> latest_pilon.fasta
  cp  pilon_bwa/pilon_polish.changes latest_pilon.changes
  """
}

process run_coverage {
  // Liczenie pokrycia dla kazdego contiga polaczona z filtorwanie odczytow
  container  = 'salmonella_illumina:2.0'
  tag "Coverage-based filtering dla sample $x"
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
  tag "Calculating basic statistics dla sample $x"
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
  tag "Predicting MLST dla sample $x"
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

process run_Seqsero {
  container  = 'salmonella_illumina:2.0'
  tag "Predicting OH dla sample $x with Seqsero"
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
  tag "Predicting OH dla sample $x with Sistr"
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
  tag "Predicting microbial resistance dla sample $x"
  publishDir "pipeline_wyniki/${x}", mode: 'copy'
  input:
  tuple val(x), path(fasta)
  output:
  // tak naprawde pheno_table jest istotne ale pozostale pliki tlumacza jakie geny niosace opornosci znaleziono i jakie mutacje
  // w genach wywolujace opornosc znaleziono
  tuple val(x), path('resfinder_out/pheno_table.txt'), path('resfinder_out/ResFinder_results.txt'), path('resfinder_out/PointFinder_results.txt')
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
  python -m resfinder -o resfinder_out -s ${params.species}  -l 0.6 -t 0.8 --acquired --point -k /opt/docker/kma/kma -db_disinf /opt/docker/disinfinder_db/ -db_res /opt/docker/resfinder_db/ -db_point /opt/docker/pointfinder_db/ -ifa ${fasta}
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
  containerOptions '--volume /mnt/sda1/michall/db/cgMLST_30042024:/cgMLST2_entero'
  tag "Predicting cgMLST dla sample $x"
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
  /data/run_blastn_ver6.sh $fasta ${task.cpus}
  cp log.log cgMLST.txt
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
  tag "Predicting genes dla sample $x"
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
  tag "Predicting VirulenceFactors dla sample $x"
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


workflow sub1 {
take:
data1
main:
data1 | view

}

workflow {

Channel
  .fromFilePairs('/mnt/sda1/michall/Salmonella/*_R{1,2}_001.fastq.gz')
  .set {initial_fastq}

//initial_fastq | check_etoki | view

processed_fastq = clean_fastq(initial_fastq)

// processed_fastq.All_path | view_output | view()

initial_scaffold = spades(processed_fastq.All_path)

// for_remaping_PE = initial_scaffold.join(processed_fastq.PE_path, by : 0, remainder : true)
// for_remaping_SE = initial_scaffold.join(processed_fastq.SE_path, by : 0, remainder : true)
//for_remaping_PE.view()
// paired_bams = bwa_paired(for_remaping_PE)
// single_bams = bwa_single(for_remaping_SE)
// merged_bams = paired_bams.join(single_bams, by : 0, remainder : true)
// merged_bams_and_scaffold = merged_bams.join(initial_scaffold,  by : 0, remainder : true)
// merged_bams_and_scaffold | view
// pilon_initial = run_pilon(merged_bams_and_scaffold)
// laczymy wyniki pilona z kana


// Puszczamy 3-krotne wygladzanie pilon-em
first_polish_run = pilon_first(initial_scaffold, processed_fastq.PE_path, processed_fastq.SE_path)
second_polish_run = pilon_second(first_polish_run, processed_fastq.PE_path, processed_fastq.SE_path)
third_polish_run = pilon_third(second_polish_run, processed_fastq.PE_path, processed_fastq.SE_path)

// Liczymy pokrycie dla ostatecznie policzonych contigow, 
// usuwamy contigi o pokryciu ponizej 0.2 sredniego globalnego pokrycia zgodnie z kodem z etoki
// params.min_coverage_ratio = 0.2
// Analize robimy jednak moim kodem

final_assembly = calculate_coverage(third_polish_run, processed_fastq.PE_path, processed_fastq.SE_path)
// final_assembly | view


// Post-analizy
extract_final_stats(final_assembly)
run_7MLST(final_assembly)
run_Seqsero(final_assembly)
run_sistr(final_assembly)
run_pointfinder(final_assembly)
run_cgMLST(final_assembly)
prokka_out = run_prokka(final_assembly)
run_VFDB(prokka_out)
}
