#chrisw 20190520
#build an image for running AKLIMATE

FROM ubuntu:18.04

ENV tbb_os linux

COPY ./dockerfile_scripts/install_apt_stuff.sh /
COPY ./dockerfile_scripts/install_r_stuff.sh /

RUN ["bash", "install_apt_stuff.sh"]
RUN ["bash", "install_r_stuff.sh"]
