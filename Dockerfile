FROM python:3.6-jessie

# install packages
RUN apt-get update && apt-get install -y \
    vim

# need to use pip for transvar
COPY ./requirements.txt /requirements.txt
RUN pip install -r /requirements.txt

# set locales or else transvar wont work correctly
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y locales

RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
    dpkg-reconfigure --frontend=noninteractive locales && \
    update-locale LANG=en_US.UTF-8

ENV LANG en_US.UTF-8 

# get conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
RUN bash ~/miniconda.sh -b -p ./miniconda
ENV PATH="/usr/local/bin:/miniconda/bin:$PATH"
# add channels
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# # get blast and dependencies
RUN conda install -y samtools
# ENV BLASTDB /annotation-station/annotation-station/data/blast_databases

# get blat
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat && mv blat /bin/
RUN chmod u+x /bin/blat

# set up working directory
COPY . /annotation-station
WORKDIR /annotation-station
RUN chmod -R 777 /annotation-station/annotation-station/data

CMD /bin/bash
