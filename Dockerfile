#chrisw 20190520
#build an image for running AKLIMATE

FROM ubuntu:18.04

ENV tbb_os linux

RUN ["mkdir", "-p", "/dockerfile_scripts/"]
COPY ./dockerfile_scripts /dockerfile_scripts

# aklimate_lib contains files from:
# https://github.com/VladoUzunangelov/junkle
# https://github.com/VladoUzunangelov/Spicer
# https://github.com/VladoUzunangelov/tcga_scripts

#RUN ["mkdir", "-p", "/aklimate_lib/"]
#COPY ./aklimate_lib /aklimate_lib

#RUN ["ln", "-s", "/aklimate_lib/repos", "/root/repos"]
RUN ["ln", "-s", "/data/repos", "/root/repos"]

RUN ["bash", "/dockerfile_scripts/install_apt_stuff.sh"]
RUN ["bash", "/dockerfile_scripts/install_r_stuff.sh"]
