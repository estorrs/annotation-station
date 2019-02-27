FROM python:3.6-jessie

# install packages
RUN apt-get update && apt-get install -y \
    vim

COPY ./requirements.txt /requirements.txt
RUN pip install -r /requirements.txt

# set locales or else transvar wont work correctly
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y locales

RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
    dpkg-reconfigure --frontend=noninteractive locales && \
    update-locale LANG=en_US.UTF-8

ENV LANG en_US.UTF-8 

# set up working directory
#USER root
#RUN useradd foo
COPY . /annotation-station
WORKDIR /annotation-station
RUN chmod -R 777 /annotation-station/annotation-station/data
#RUN chown -R root:root /annotation-station

CMD /bin/bash
