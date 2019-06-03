#chrisw 20190520
#build an image for running AKLIMATE

FROM ubuntu:18.04

ENV tbb_os linux

COPY ./dockerfile_scripts/install_apt_stuff.sh /
COPY ./dockerfile_scripts/install_r_stuff.sh /
COPY ./run_aklimate.R /

# These are from https://github.com/VladoUzunangelov/junkle
COPY ./junkle-utils.R /
COPY ./junkle.R /

# need to also copy utils.R, which was at /home/ubuntu/repos/tcga_scripts/utils.R in the original VM
# It is found at https://github.com/VladoUzunangelov/tcga_scripts
RUN ["mkdir", "-p", "/home/ubuntu/repos/tcga_scripts/"]
COPY ./utils.R /home/ubuntu/repos/tcga_scripts/

# Spicer repo is at https://github.com/VladoUzunangelov/Spicer
RUN ["mkdir", "-p", "/home/ubuntu/repos/Spicer/"]
COPY ./Spicer /home/ubuntu/repos/Spicer

RUN ["bash", "install_apt_stuff.sh"]
RUN ["bash", "install_r_stuff.sh"]
