#!/bin/bash

#chrisw 20190520
#install various packages for running AKLIMATE

apt update
apt upgrade -y
apt dist-upgrade -y
apt autoremove -y


# to avoid input prompt for package tzdata
DEBIAN_FRONTEND='noninteractive' apt install -y r-base

apt-get install -y xorg xvfb xauth xfonts-base

# install INTEL MKL library
git clone https://github.com/eddelbuettel/mkl4deb.git
bash ./mkl4deb/script.sh

apt install -y libcurl4-openssl-dev libssl-dev libxml2-dev

# required for Similarity R package
apt install -y libtbb-dev
# export tbb_os=linux
