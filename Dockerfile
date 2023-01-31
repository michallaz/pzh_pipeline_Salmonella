FROM nvidia/cuda:11.2.2-cudnn8-runtime-ubuntu20.04
LABEL maintainer "Michal Lazniewski <mlazniewski@pzh.gov.pl>"

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
#ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && \
    apt-get install python3.8 python3-pip default-jre curl unzip wget vim bc htop git gcc zlib1g-dev libbz2-dev libcurl4-gnutls-dev libssl-dev liblzma-dev python3-pycurl screen -y

RUN pip install ete3 numba numpy==1.23.4 pandas==1.0.5 biopython==1.73 psutil pysam
RUN pip install sklearn
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
RUN cat modules/configure.py | sed s/"subprocess.Popen('curl -Lo pilon-1.23.jar {0}'.split(), stderr=subprocess.PIPE).communicate()"/"subprocess.Popen('curl -Lo pilon-1.23.jar {0}'.format(url).split(), stderr=subprocess.PIPE).communicate()"/g > modules/configure2.py
RUN cp modules/configure2.py modules/configure.py

# nie sciagamy kraken-a tutaj ale w wersji finalnej bez tego sie nie obejdzie 
# sciaganie bazy kraken-owej opisane jest w /home/michall/kraken2/README_KRAKEN_INSTALL
# tutaj tworzymy tylko katalog i podczas konfigurowania Etoki mowimy mu aby tam zagladal w trakcie pracy
# przy starcie kontenera po prostu mapujemy baze z dysku na kontener przez "-v"
# baza krakenowa wazy 240 Gb  
RUN mkdir -p /kraken2_sdb
RUN python EToKi.py configure --install --usearch /opt/docker/usearch11.0.667_i86linux32 --link_krakenDB /kraken2_sdb/

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
RUN . ~/.bashrc
# Sciaganie MLST 7 genowy ze strony enterobase - https://enterobase.warwick.ac.uk/schemes/Salmonella.Achtman7GeneMLST/
RUN mkdir -p /Achtman7GeneMLST_entero
WORKDIR /Achtman7GeneMLST_entero

RUN curl -s  https://enterobase.warwick.ac.uk/schemes/Salmonella.Achtman7GeneMLST/ >> log
RUN cat log  | grep fasta | grep -v MLST | sed s'/<\|>/ /'g | awk '{print $3}' >> to_download.txt
RUN for K in `cat to_download.txt`; do wget https://enterobase.warwick.ac.uk/schemes/Salmonella.Achtman7GeneMLST/${K}; gunzip ${K}; done
RUN rm log to_download.txt
RUN cat *fasta >> all_allels.fasta
# Tworzenie i indeksowanie bazy bazy 
RUN wget https://enterobase.warwick.ac.uk/schemes/Salmonella.Achtman7GeneMLST/MLST_Achtman_ref.fasta
RUN /opt/docker/EToKi/EToKi.py MLSTdb -i all_allels.fasta -s MLST_Achtman_ref.fasta -r references.fasta -d MLST_database.tab

#Sciaganie cgMLST ze strony enterobase - https://enterobase.warwick.ac.uk/schemes/Salmonella.cgMLSTv2/

RUN mkdir -p /cgMLST2_entero
WORKDIR /cgMLST2_entero

RUN wget https://enterobase.warwick.ac.uk/schemes/Salmonella.cgMLSTv2/cgMLST_salmonella.tar
RUN tar -xf cgMLST_salmonella.tar
RUN gunzip *gz
RUN for K in `ls *fasta`; do /opt/docker/EToKi/externals/ncbi-blast-2.8.1+/bin/makeblastdb -in ${K} -dbtype nucl; done
COPY all_3002_allele.txt /cgMLST2_entero/
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

## budowanie profilu na podstawie cgMLST 
WORKDIR /cgMLST2_entero
RUN pHierCC -p profiles.list -o profile_output -n 60

# Metaphlan
WORKDIR /opt/docker
RUN git clone https://github.com/biobakery/MetaPhlAn.git 
WORKDIR /opt/docker/MetaPhlAn
RUN pip install .
## Baza do Metaphlan /uwaga 24 Gb/
### Uwaga bazy mozna pobraz tez bezposrednio ze strony http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/
### baza buduje sie dobre 2h
RUN mkdir -p /bowtie_db
RUN metaphlan --install --nproc 60 --bowtie2db /bowtie_db

WORKDIR /data
COPY all_functions_salmonella.py run_blastn_ver5.sh master_script_kontener.sh prep_hierCC.py /data/
#ENTRYPOINT [ "/usr/bin/tini", "--" ]
RUN pip install pandas==1.0.5 
# to wyzej powinno rozwalic resfindera ( resfinder 4.2.5 has requirement pandas>=1.4.2 ), ale dziala tak i tak za to sistr bez wersji 1.4.2 nie pojdzie
CMD ["/bin/bash"]
#ENTRYPOINT ["/bin/bash", "/SARS-CoV2/scripts/master_sript_docker.sh"]
