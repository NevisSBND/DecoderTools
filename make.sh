#!/bin/sh
rootcint -f decoder_dict.cc -c decoder.hh LinkDef.h
g++ decoder_dict.cc decoder.cc -Wall -o decoder.exe `root-config --cflags  --glibs`

