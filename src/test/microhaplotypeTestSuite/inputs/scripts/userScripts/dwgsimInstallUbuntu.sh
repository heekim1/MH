#!/bin/sh
echo 'Installing DWGSIM'
git clone https://github.com/nh13/DWGSIM.git
cd DWGSIM
git submodule init
git submodule update
echo 'Installing dependencies'
sudo apt-get install libncurses5-dev
make
echo 'Done'
