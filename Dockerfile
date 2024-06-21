FROM nvidia/cuda:11.2.2-cudnn8-runtime-ubuntu20.04
LABEL maintainer "Michal Lazniewski <mlazniewski@pzh.gov.pl>"

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
#ENV PATH /opt/conda/bin:$PATH



RUN apt-get update --fix-missing && apt upgrade -y && \
    apt-get install build-essential pkgconf python3.8 python3-pip default-jre curl unzip wget vim bc htop git gcc zlib1g-dev libbz2-dev libcurl4-gnutls-dev libssl-dev liblzma-dev python3-pycurl screen iqtree -y

RUN pip install ete3 numba numpy==1.23.4 pandas==1.0.5 biopython==1.73 psutil pysam cgecore packaging tables==3.7.0 h5py
RUN pip install scikit-learn
# bez podania konkretnego numpy nie instaluje sie numba, sklearn w trakcie budowania zwraca error, a nastepnie i tak sie poprawnie instaluje
WORKDIR /usr/bin
RUN ln -s python3 python

RUN mkdir -p /opt/docker
WORKDIR /opt/docker

RUN wget https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz
RUN gunzip usearch11.0.667_i86linux32.gz
RUN chmod a+x usearch11.0.667_i86linux32

RUN git clone https://github.com/zheminzhou/EToKi.git
WORKDIR /opt/docker/EToKi

# podmieniamy bledna linijke przy instalacji pilon przez etoki
#RUN cat modules/configure.py | sed s/"subprocess.Popen('curl -Lo pilon-1.23.jar {0}'.split(), stderr=subprocess.PIPE).communicate()"/"subprocess.Popen('curl -Lo pilon-1.23.jar {0}'.format(url).split(), stderr=subprocess.PIPE).communicate()"/g > modules/configure2.py
#RUN cp modules/configure2.py modules/configure.py

# nie sciagamy kraken-a tutaj ale w wersji finalnej bez tego sie nie obejdzie 
# sciaganie bazy kraken-owej opisane jest w /home/michall/kraken2/README_KRAKEN_INSTALL
# tutaj tworzymy tylko katalog i podczas konfigurowania Etoki mowimy mu aby tam zagladal w trakcie pracy
# przy starcie kontenera po prostu mapujemy baze z dysku na kontener przez "-v"
# baza krakenowa wazy 240 Gb  
RUN mkdir -p /kraken2_sdb
RUN python EToKi.py configure --install --usearch /opt/docker/usearch11.0.667_i86linux32 --link_krakenDB /kraken2_sdb/ --path iqtree=/usr/bin/iqtree

#SEqsero
WORKDIR /opt/docker
RUN git clone https://github.com/denglab/SeqSero2.git
WORKDIR  /opt/docker/SeqSero2
RUN python3 -m pip install --user .


#SISTR
WORKDIR /opt/docker
RUN git clone https://github.com/phac-nml/sistr_cmd
WORKDIR /opt/docker/sistr_cmd
RUN python setup.py install 

#MAFFT
WORKDIR /opt/docker
RUN wget https://mafft.cbrc.jp/alignment/software/mafft-7.505-with-extensions-src.tgz
RUN tar -zxf mafft-7.505-with-extensions-src.tgz
WORKDIR /opt/docker/mafft-7.505-with-extensions/core
RUN make
RUN make install
WORKDIR /opt/docker/mafft-7.505-with-extensions/extensions
RUN make
RUN make install

# Dodanie do PATH i source 
RUN echo 'export PATH="/opt/docker/EToKi:$PATH"' >> ~/.bashrc
RUN echo 'export PATH="/opt/docker/SeqSero2/bin:$PATH"'  >> ~/.bashrc
RUN echo 'export PATH="/opt/docker/EToKi/externals/ncbi-blast-2.8.1+/bin:$PATH"' >> ~/.bashrc
RUN echo 'export PATH="/opt/docker/EToKi/externals/:$PATH"' >> ~/.bashrc

# Sciaganie MLST 7 genowy ze strony enterobase - https://enterobase.warwick.ac.uk/schemes/Salmonella.Achtman7GeneMLST/
#RUN mkdir -p /Achtman7GeneMLST_entero
# WORKDIR /Achtman7GeneMLST_entero

# ENTERO zmineilo aPI na razie korzystam z wersji jako mam lokalnie i ja kopiuje

#COPY Achtman7GeneMLST_entero /Achtman7GeneMLST_entero

# uzanjemy ze schemat jest NA ZEWNATRZ, latwiej bedzie updatowac plik z profilami

#RUN curl -s  https://enterobase.warwick.ac.uk/schemes/Salmonella.Achtman7GeneMLST/ >> log
#RUN cat log  | grep fasta | grep -v MLST | sed s'/<\|>/ /'g | awk '{print $3}' >> to_download.txt
#RUN for K in `cat to_download.txt`; do wget https://enterobase.warwick.ac.uk/schemes/Salmonella.Achtman7GeneMLST/${K}; gunzip ${K}; done
#RUN rm log to_download.txt
#RUN cat *fasta >> all_allels.fasta
# Tworzenie i indeksowanie bazy bazy 
#RUN wget https://enterobase.warwick.ac.uk/schemes/Salmonella.Achtman7GeneMLST/MLST_Achtman_ref.fasta

# Moja baza juz jest przygotowana nie ma potrzeby jej tworzyc 

# RUN /opt/docker/EToKi/EToKi.py MLSTdb -i all_allels.fasta -s MLST_Achtman_ref.fasta -r references.fasta -d MLST_database.tab

#Sciaganie cgMLST ze strony enterobase - https://enterobase.warwick.ac.uk/schemes/Salmonella.cgMLSTv2/

### cgMLST tez bedzie montowane z zewnatrz ### 

# RUN mkdir -p /cgMLST2_entero
# WORKDIR /cgMLST2_entero

#RUN wget https://enterobase.warwick.ac.uk/schemes/Salmonella.cgMLSTv2/cgMLST_salmonella.tar # nie korzystamy archiwum nie bylo
# updatowane od 2017 roku
#RUN tar -xf cgMLST_salmonella.tar

### Na razie z tego nie korzystam wiec nie sciagam bazy
# RUN curl -s  https://enterobase.warwick.ac.uk/schemes/Salmonella.cgMLSTv2/ >> log
# RUN cat log  | grep fasta | grep -v MLST | sed s'/<\|>/ /'g | awk '{print $3}' >> to_download.txt
# RUN echo profiles.list.gz >> to_download.txt
# RUN for K in `cat to_download.txt`; do wget https://enterobase.warwick.ac.uk/schemes/Salmonella.cgMLSTv2/${K}; gunzip ${K}; done

# RUN for K in `ls *fasta`; do /opt/docker/EToKi/externals/ncbi-blast-2.8.1+/bin/makeblastdb -in ${K} -dbtype nucl; done
# RUN cat to_download.txt | grep fasta |   sed s'/.gz//'g >> all_3002_allele.txt
# RUN sed -i  s'/.fasta//' all_3002_allele.txt
# RUN rm log to_download.txt
# RUN cat *fasta >> all_allels.fasta

#COPY all_3002_allele.txt /cgMLST2_entero/
#RUN wget https://enterobase.warwick.ac.uk/schemes/Salmonella.cgMLSTv2/profiles.list.gz
#RUN gunzip profiles.list.gz

#instalacja baz point/res finder
WORKDIR /opt/docker
## kma sluzy do indeksowania baz w tym soft
RUN git clone https://bitbucket.org/genomicepidemiology/kma.git
RUN cd kma && make
## bazy
### Pointfinder
WORKDIR /opt/docker
RUN git clone https://bitbucket.org/genomicepidemiology/pointfinder_db/
WORKDIR /opt/docker/pointfinder_db
RUN ls -lrt
RUN python3 INSTALL.py /opt/docker/kma/kma_index non_interactive
WORKDIR /opt/docker

### disinfinde
RUN git clone https://bitbucket.org/genomicepidemiology/disinfinder_db/
WORKDIR /opt/docker/disinfinder_db
RUN python3 INSTALL.py /opt/docker/kma/kma_index non_interactive
WORKDIR /opt/docker
### Resfinder

RUN git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git
WORKDIR /opt/docker/resfinder_db
RUN python3 INSTALL.py /opt/docker/kma/kma_index non_interactive
WORKDIR /opt/docker

### Plasmidfinder
RUN git clone https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git
WORKDIR /opt/docker/plasmidfinder_db
RUN python INSTALL.py /opt/docker/kma/kma_index non_interactive
WORKDIR /opt/docker

### Wyspy wirulencji
RUN git clone https://bitbucket.org/genomicepidemiology/spifinder_db.git
WORKDIR /opt/docker/spifinder_db
RUN python3 INSTALL.py /opt/docker/kma/kma_index non_interactive
WORKDIR /opt/docker

# instalacja samego resfindera
RUN pip install resfinder
## dodawanie sciezki do baz resfindera i do kma
RUN echo 'export PATH="/opt/docker/kma:$PATH"' >> ~/.bashrc
RUN echo 'export CGE_RESFINDER_RESGENE_DB="/opt/docker/resfinder_db"' >> ~/.bashrc
RUN echo 'export CGE_RESFINDER_RESPOINT_DB="/opt/docker/pointfinder_db"' >>  ~/.bashrc
RUN echo 'export CGE_DISINFINDER_DB="/opt/docker/disinfinder_db"' >> ~/.bashrc

# instalacja spifinder-a 
RUN git clone https://bitbucket.org/genomicepidemiology/spifinder.git
RUN chmod a+x spifinder/spifinder.py
## dodanie spifindera do PATH
RUN echo 'export PATH="/opt/docker/spifinder:$PATH"' >> ~/.bashrc
# Hiercc
RUN pip install pHierCC

## modyfikacja pHierCC w celu szybkiego doliczania profili
## skrypty poki co w tym samym katalogu co dockerfile
##  uwaga przy zmienionej dystrybucji python moze byc problem z ustaleniem gdzie pip unstaluje hierCC, wiec ustawile na sztywno sciaganie python 3.8

RUN mkdir -p /opt/docker/pHierCC_tmp
WORKDIR /opt/docker/pHierCC_tmp
COPY getDistance.py pHierCC.py  /opt/docker/pHierCC_tmp/
RUN cp pHierCC.py /usr/local/lib/python3.8/dist-packages/pHierCC/pHierCC.py
RUN cp getDistance.py /usr/local/lib/python3.8/dist-packages/pHierCC/getDistance.py

### budowanie profilu na podstawie cgMLST 
### UWAGA na razie pHierCC nie dziala dla wiecej niz 100k profili (teraz w entero jest ponad 300k), dla testow puszczam 2 wersje 
### w celu weryfikacji odtwarzania kodu po przepisaniu go na dask-a (jesli sie uda)

### NIE ROBIMY DBAZA JEST ZA DUZA 
WORKDIR /cgMLST2_entero
# RUN head -70000 profiles.list >> profiles_70k.list
# RUN head -20 profiles.list >> profiles_20.list
# RUN pHierCC -p profiles_70k.list -o profile_output_70k -n 30
# RUN pHierCC -p profiles_20.list -o profile_output_20 -n 1

### Metaphlan
### NIE INSTALUJEMY POKI CO ZA DUZO DANYCH 
# WORKDIR /opt/docker
# RUN git clone https://github.com/biobakery/MetaPhlAn.git 
# WORKDIR /opt/docker/MetaPhlAn
# RUN pip install .
### Baza do Metaphlan /uwaga 24 Gb/
### Uwaga bazy mozna pobraz tez bezposrednio ze strony http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/
### baza buduje sie dobre 2h
# RUN mkdir -p /bowtie_db
### dodajemy do PATH w trakcie budowania kontenera ta sciezke aby widziany byl bowtie-build
# ENV PATH /opt/docker/EToKi/externals/:$PATH
# RUN metaphlan --install --nproc 60 --bowtie2db /bowtie_db


# cgMLST z CGE ? wtedy musimy zainstalowac prodigal github czy apt ? budowanie z githuba zadzialalo bez problemu
WORKDIR /opt/docker
RUN git clone https://bitbucket.org/genomicepidemiology/cgmlstfinder.git
RUN git clone https://bitbucket.org/genomicepidemiology/cgmlstfinder_db.git
RUN git clone https://github.com/hyattpd/Prodigal.git

# budowanie bazy dla cgmlstfindera w oparciu o najnowsze pliki z enterobase
WORKDIR /opt/docker/cgmlstfinder_db/scripts
RUN mkdir -p /opt/docker/cgmlstfinder_db/scripts/custom_db
RUN mkdir -p /opt/docker/cgmlstfinder_db/scripts/custom_out

# Tego NIE ROBIMY POKI CO BO ZZERA 80GB
# RUN python3 cgmlst_dl.py -o /opt/docker/cgmlstfinder_db/scripts/custom_out -db /opt/docker/cgmlstfinder_db/scripts/custom_db  -k /opt/docker/kma/kma_index

## to sciaga baze ze strony CGE
## tak tylko kod sciagnal baze z 18 pazdziernika, wiec lepiej budowac baze samemu na podstawie najnowszej wersji z enterobase co zrobiono wyzej
# RUN python3 INSTALL.py -s salmonella 

# instalcja prodigal
WORKDIR /opt/docker/Prodigal
RUN make install
RUN echo 'export PATH="/opt/docker/cgmlstfinder:$PATH"' >> ~/.bashrc

WORKDIR /data
# COPY all_functions_salmonella.py run_blastn_ver6.sh  run_blastn_ver7.sh master_script_kontener.sh prep_hierCC.py /data/
#ENTRYPOINT [ "/usr/bin/tini", "--" ]
RUN pip install pandas==1.0.5 xlrd
# to wyzej powinno rozwalic resfindera ( resfinder 4.2.5 has requirement pandas>=1.4.2 ), ale dziala tak i tak za to sistr bez wersji 1.4.2 nie pojdzie
# xlrd potrebne jest to czytania xls a baza VFDB ma plik z opisami w arkuszu excelowym

#pobieranie pilon-a
RUN wget -O /opt/docker/EToKi/externals/pilon.jar https://github.com/broadinstitute/pilon/releases/download/v1.24/pilon-1.24.jar

# dociaganie bazy sistr

RUN wget https://sairidapublic.blob.core.windows.net/downloads/sistr/database/SISTR_V_1.1_db.tar.gz
RUN tar -zxvf SISTR_V_1.1_db.tar.gz  -C /usr/local/lib/python3.8/dist-packages/sistr; rm SISTR_V_1.1_db.tar.gz
RUN touch  /usr/local/lib/python3.8/dist-packages/sistr/dbstatus.txt
ENV PATH $PATH:/opt/docker/EToKi/externals/

# Etoki ma wlasnego blasta ale sciagamy nowa werszje
# 
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.15.0/ncbi-blast-2.15.0+-x64-linux.tar.gz
RUN mkdir -p /blast
RUN tar -zxf ncbi-blast-2.15.0+-x64-linux.tar.gz -C /blast --strip-components=1
# moj skrypt run_blast explicite bedzie puszczal blasta z /blast/bin wiec nie ma potrzeby dodawania do PATH
# poza tym nie wiem jak inne kompnenty zareaguja na 2 wersje blasta sciagana przeze mnie i ta z Etoki/externals
# btw na pewno podczas instalacji etoki mozna wskazac gdzie jest blast wiec 
# docelowo blast z /externals moze byc blastem systemowym
# ENV PATH="/blast/bin:$PATH"

# Etoki sciaga tez sam flye, ale wole 3mac najnowsza wersje
RUN git clone https://github.com/fenderglass/Flye
WORKDIR /data/Flye
RUN make 
# NIE uzywac -j  przy make bo nie buduja sie wszystkie binarki
# NIE dodajemy do path bo jest juz flye z Etoki w PATH 
#ENV PATH $PATH:/data/Flye/bin  
WORKDIR /data
# Pozostale skrypty
COPY coverage_filter.py calculate_stats.py run_VFDB.sh /opt/docker/EToKi/externals

## DO INSTALACJI CANU POTRZEBA DOINSTALOWAC apt install pkgconf
## Nanopolish odrzucamy

WORKDIR /data
RUN git clone https://github.com/marbl/canu.git
## RUN git clone --recursive https://github.com/jts/nanopolish.git

WORKDIR /data/canu/src
RUN make -j 20
## WORKDIR /data/nanopolish
## RUN pip install -r /data/nanopolish/scripts/requirements.txt
## Nanopolish ma specjalnie bez opcji -j bo potrafi sie nie budowac obraz przy wielu procesorach
## RUN make 

COPY all_functions_salmonella.py run_blastn_ver8.sh  run_blastn_ver7.sh master_script_kontener.sh prep_hierCC.py /data/
WORKDIR /data
CMD ["/bin/bash"]
#ENTRYPOINT ["/bin/bash", "/SARS-CoV2/scripts/master_sript_docker.sh"]
