#!/bin/sh
echo -e "Generating dictionary of decoder.hh and HeaderInfo.hh\n"
rootcint -f decoder_dict.cc -c decoder.hh HeaderInfo.hh LinkDef.h
echo -e "Compiling decoder.cc\n"
g++ decoder_dict.cc decoder.cc -Wall -o decoder.exe `root-config --cflags  --glibs`
echo -e "Compiling plotter.cc\n"
g++ decoder_dict.cc plotter.cc -Wall -o plotter.exe `root-config --cflags  --glibs`
echo -e "Compiling channel_mapper.cc\n"
g++ decoder_dict.cc channel_mapper.cc -Wall -o channel_mapper.exe `root-config --cflags  --glibs`
