#!/bin/sh
echo "Generating dictionary of decoder.hh and HeaderInfo.hh\n"
#rootcint -f decoder_dict.cc -c decoder.hh HeaderInfo.hh LinkDef.h
echo "Compiling decoder.cc\n"
#g++ decoder_dict.cc decoder.cc -Wall -o decoder.exe `root-config --cflags  --glibs`
echo "Compiling plotter.cc\n"
g++ plotter.cc -Wall -o plotter.exe `root-config --cflags  --glibs`

