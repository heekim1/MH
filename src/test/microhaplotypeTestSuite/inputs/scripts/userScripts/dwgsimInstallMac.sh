#!/bin/sh
echo 'Installing DWGSIM'
git clone https://github.com/nh13/DWGSIM.git
cd DWGSIM
git submodule init
git submodule update
make
echo 'Done'
