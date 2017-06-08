#!/bin/sh
echo "Generating dictionary\n"
rootcint -f dissecter_dict.cc -c dissecter.hh FEMInfo.hh XMITInfo.hh LinkDef.h
echo "Compiling dissecter.cc\n"
g++ dissecter_dict.cc dissecter.cc -Wall -o dissecter.exe `root-config --cflags  --glibs`
