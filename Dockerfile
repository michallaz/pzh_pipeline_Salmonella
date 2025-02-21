FROM nvidia/cuda:11.2.2-cudnn8-runtime-ubuntu20.04
LABEL maintainer "Michal Lazniewski <mlazniewski@pzh.gov.pl>"

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN apt-get update --fix-missing && apt upgrade -y && \
    apt-get install build-essential pkgconf python3.8 python3-pip default-jre curl unzip wget vim bc htop git gcc zlib1g-dev libbz2-dev libcurl4-gnutls-dev libssl-dev liblzma-dev python3-pycurl screen iqtree mash seqtk -y

RUN pip install ete3 numba numpy==1.23.4 biopython==1.73 psutil pysam cgecore packaging tables==3.7.0 h5py xlrd 
### xlrd - allows to read excels, need that for VFDB
### Not sure why I choose these version of some libraries  

RUN pip install scikit-learn 
RUN pip install fiona==1.9.6 geodatasets geopandas plotly

WORKDIR /usr/bin
RUN ln -s python3 python

RUN mkdir -p /opt/docker
WORKDIR /opt/docker

RUN wget https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz ;\
    gunzip usearch11.0.667_i86linux32.gz ;\
    chmod a+x usearch11.0.667_i86linux32


# EToKI

### Szczatkowo, ale ciagle z niego korzystamy
RUN git clone https://github.com/zheminzhou/EToKi.git
WORKDIR /opt/docker/EToKi

### nie sciagamy kraken-a tutaj ale w wersji finalnej bez tego sie nie obejdzie 
### sciaganie bazy kraken-owej opisane jest w /home/michall/kraken2/README_KRAKEN_INSTALL
### tutaj tworzymy tylko katalog i podczas konfigurowania Etoki mowimy mu aby tam zagladal w trakcie pracy
### przy starcie kontenera po prostu mapujemy baze z dysku na kontener przez "-v"
### baza krakenowa wazy 240 Gb  

RUN mkdir -p /kraken2_sdb
RUN python EToKi.py configure --install --usearch /opt/docker/usearch11.0.667_i86linux32 --link_krakenDB /kraken2_sdb/ --path iqtree=/usr/bin/iqtree

ENV PATH="/opt/docker/EToKi:/opt/docker/EToKi/externals/ncbi-blast-2.8.1+/bin:/opt/docker/EToKi/externals/:$PATH"

## Updating pilon
RUN wget -O /opt/docker/EToKi/externals/pilon.jar https://github.com/broadinstitute/pilon/releases/download/v1.24/pilon-1.24.jar


# SeqSero2
WORKDIR /opt/docker
RUN git clone https://github.com/denglab/SeqSero2.git; cd /opt/docker/SeqSero2; python3 -m pip install --user . ; cd /opt/docker

# MAFFT
RUN wget https://mafft.cbrc.jp/alignment/software/mafft-7.505-with-extensions-src.tgz ;\
    tar -zxf mafft-7.505-with-extensions-src.tgz; \
    cd /opt/docker/mafft-7.505-with-extensions/core ;\
    make; make install; \
    cd /opt/docker/mafft-7.505-with-extensions/extensions ;\
    make; make install ;\
    cd /opt/docker


# Various CGE programs

## kma (required by most CGE databases)

RUN git clone https://bitbucket.org/genomicepidemiology/kma.git; cd kma ; make; cd /opt/docker

## CGE databases


#RUN git clone https://bitbucket.org/genomicepidemiology/pointfinder_db/; \
#    git clone https://bitbucket.org/genomicepidemiology/disinfinder_db/; \
#    git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git; \
#    git clone https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git; \
#    git clone --depth 1 https://git@bitbucket.org/genomicepidemiology/virulencefinder_db.git; \
#    git clone https://bitbucket.org/genomicepidemiology/spifinder_db.git; \
#    git clone https://bitbucket.org/genomicepidemiology/cgmlstfinder_db.git ;\
#    git clone https://bitbucket.org/genomicepidemiology/mlst_db.git
    

## "Installing" database (save cgmlstfinder_db that we dont use for now )

#RUN cd /opt/docker/pointfinder_db; python3 INSTALL.py /opt/docker/kma/kma_index non_interactive; \
#    cd /opt/docker/disinfinder_db; python3 INSTALL.py /opt/docker/kma/kma_index non_interactive; \
#    cd /opt/docker/resfinder_db; python3 INSTALL.py /opt/docker/kma/kma_index non_interactive; \
#    cd /opt/docker/plasmidfinder_db; python INSTALL.py /opt/docker/kma/kma_index non_interactive; \
#    cd /opt/docker/virulencefinder_db; python INSTALL.py /opt/docker/kma/kma_index non_interactive; \
#    cd /opt/docker/spifinder_db; python3 INSTALL.py /opt/docker/kma/kma_index non_interactive; \
#    cd /opt/docker/mlst_db; python INSTALL.py /opt/docker/kma/kma_index non_interactive; \
#    cd /opt/docker

## CGE programs

### Most are just python scripts
ENV VIRULENCEFINDER_VERSION 2.0.5

RUN git clone https://bitbucket.org/genomicepidemiology/plasmidfinder.git; \
    git clone https://bitbucket.org/genomicepidemiology/kmerfinder.git; \
    git clone https://bitbucket.org/genomicepidemiology/spifinder.git; \
    git clone https://bitbucket.org/genomicepidemiology/cgmlstfinder.git; \
    git clone https://bitbucket.org/genomicepidemiology/mlst.git ;\
    pip install resfinder ;\
    git clone https://bitbucket.org/genomicepidemiology/virulencefinder.git; cd virulencefinder; git checkout tags/${VIRULENCEFINDER_VERSION} ;\
    cd /opt/docker
 
ENV PATH="/opt/docker/kma:$PATH"
#ENV CGE_RESFINDER_RESGENE_DB="/opt/docker/resfinder_db"
#ENV CGE_RESFINDER_RESPOINT_DB="/opt/docker/pointfinder_db"
#ENV CGE_DISINFINDER_DB="/opt/docker/disinfinder_db"

# Metaphlan

### Metaphlan database is downloaded elsewhere
RUN git clone https://github.com/biobakery/MetaPhlAn.git; cd /opt/docker/MetaPhlAn; pip install .; cd /opt/docker 

# Prodigal (alternative to Prokka)
RUN git clone https://github.com/hyattpd/Prodigal.git; cd /opt/docker/Prodigal; \
    make install; \
    cd /opt/docker

# Sistr

RUN git clone https://github.com/phac-nml/sistr_cmd; cd /opt/docker/sistr_cmd; python setup.py install; cd /opt/docker
 
## Updating Sistr database
### Sistr only reads database from specific loaction within filesystem amd this cannot be changed

RUN wget https://sairidapublic.blob.core.windows.net/downloads/sistr/database/SISTR_V_1.1_db.tar.gz ;\
    tar -zxvf SISTR_V_1.1_db.tar.gz -C /usr/local/lib/python3.8/dist-packages/sistr ;\
    rm SISTR_V_1.1_db.tar.gz ;\
    touch  /usr/local/lib/python3.8/dist-packages/sistr/dbstatus.txt ;\
    cd /opt/docker


#  NCBI-Blast

### To ensure that this blast is used (not the one installed via EToki) one must explixitly call binaries from /blast/bin/
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.15.0/ncbi-blast-2.15.0+-x64-linux.tar.gz ;\
    mkdir -p /blast ;\
    tar -zxf ncbi-blast-2.15.0+-x64-linux.tar.gz -C /blast --strip-components=1 ;\
    cd /opt/docker

# Flye

### Do not use make -j, some binaries are not properly created
### Like blast to use this version of Flye (not one from EToki) call /opt/docker/Flye/bin/Flye explicitly

RUN git clone https://github.com/fenderglass/Flye ;\
    cd /opt/docker/Flye ;\
    make ;\
    cd /opt/docker

# Canu and nanopolish

###  NOT used by the pipeline, at least for now, but maybe in the future

# WORKDIR /data
# RUN git clone https://github.com/marbl/canu.git ;\
#    cd /opt/docker/canu/src/ ;\
#    make ;\
#    sed -i 's/\$java = `command -v \$java`/\$java = `which \$java`/'  /opt/docker/canu/build/bin/../lib/perl5/site_perl/canu/Defaults.pm ;\
#    sed -i 's/--threads 2/--threads 1/g'  /opt/docker/canu/build/bin/../lib/perl5/site_perl/canu/Consensus.pm ;\
#    cd /opt/docker 

### canu ma problem z wykryciem javy przez blad w jego skrypcie
### /opt/docker/canu/build/bin/../lib/perl5/site_perl/canu/Defaults.pm
### linijka 1082 jest 
### $java = `command -v $java`;  #  See Execution.pm getBinDirectoryShellCode()
### powoduje ze polecenie "java" nie jest wywolywane, bo "command", polecenie shellai, u mnie w docker z jakiego powodu nie dziala
### wystarczy ta linijke  zamienic na
### $java = `which $java`
### Przynajmniej poki autorzy canu nie zauwaza bledu, choc nie wiem czy nie jest to kwestia mojego dockera



# RUN git clone --recursive https://github.com/jts/nanopolish.git

# WORKDIR /data/canu/src
# RUN make -j 20
# WORKDIR /data/nanopolish
# RUN pip install -r /data/nanopolish/scripts/requirements.txt
# Nanopolish ma specjalnie bez opcji -j bo potrafi sie nie budowac obraz przy wielu procesorach
# RUN make 

# Medaka

RUN pip install urllib3==1.26.16 pysam==0.21.0 parasail==1.3.4 ;\
    pip install pyabpoa==1.4.0 ;\
    apt update ;\
    apt install python3.8-venv libncurses5-dev -y

RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash
RUN apt install git-lfs ;\
    git lfs install ;\
    git clone https://github.com/nanoporetech/medaka.git ;\
    cd /opt/docker/medaka ;\
    make install_root ;\
    cd /opt/docker


# FASTQC

### Pasted from SARS-CoV2 dockerfile

ARG FASTQC_VERSION="0.12.1"
WORKDIR /opt/docker
RUN curl -fsSL "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip" -o "fastqc_v${FASTQC_VERSION}.zip" && \
    unzip "fastqc_v${FASTQC_VERSION}.zip" && \
    rm "fastqc_v${FASTQC_VERSION}.zip"
ENV PATH="/opt/docker/FastQC:$PATH"

# AMRFIDER

## HAMMER

WORKDIR /opt/docker
RUN cd /opt/docker ;\
    wget http://eddylab.org/software/hmmer/hmmer-3.3.2.tar.gz; tar -zxf hmmer-3.3.2.tar.gz ;\
    cd /opt/docker/hmmer-3.3.2 ;\
    ./configure; make; make install ;\
    cd /opt/docker 
 
## AMRFIDER
RUN git clone https://github.com/ncbi/amr.git ;\
    cd /opt/docker/amr ;\
    git reset --hard 25a3690 ;\
    git checkout stxtype ;\
    git submodule update --init ;\
    make; make install ;\
    cd /opt/docker

# ECTyper
# Fixed version 2.0 commit 21f2cbd
RUN git clone https://github.com/phac-nml/ecoli_serotyping.git ;\
    cd /opt/docker/ecoli_serotyping ;\
    git reset --hard 21f2cbd ;\
    pip3 install . ;\
    ectyper_init ;\
    cd /opt/docker

## ECTyper databases, obsolete used by ectyper v.1.0

#RUN wget -O /usr/local/lib/python3.8/dist-packages/ectyper-1.0.0-py3.8.egg/ectyper/Data/assembly_summary_refseq.txt http://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt ;\  
#    wget -O  /usr/local/lib/python3.8/dist-packages/ectyper-1.0.0-py3.8.egg/ectyper/Data/refseq.genomes.k21s1000.msh https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh


# Custom scripts 
COPY bin/all_functions_salmonella.py bin/run_blastn_ver11.sh bin/parse_fastqc_output.py geojeson/ne_10m_admin_0_countries.geojson bin/plot_historical_data_plotly.py  bin/parse_kmerfinder.py /data/
COPY bin/json_output_contaminations.py bin/coverage_filter.py bin/calculate_stats.py bin/run_VFDB.sh bin/run_fastqc_and_generate_json.py bin/initial_mlst_parser.py bin/extract_final_stats_parser.py bin/resfinder_parser.py bin/amrfinder_parser.py bin/plasmidfinder_parser.py bin/sistr_parser.py bin/seqsero_parser.py bin/ectyper_parser.py bin/spifinder_parser.py bin/vfdb_parser.py /bin/virulencefinder_parser.py bin/prepare_full_json.py /opt/docker/EToKi/externals

CMD ["/bin/bash"]
