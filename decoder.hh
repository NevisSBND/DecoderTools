#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

// Types of words
enum WordType{
  kHeaderFirst,
  kHeaderIDSlot,
  kHeaderNWordsMSB,
  kHeaderNWordsLSB,
  kHeaderEventMSB,
  kHeaderEventLSB,
  kHeaderFrameMSB,
  kHeaderFrameLSB,
  kHeaderChecksumMSB,
  kHeaderChecksumLSB,
  kHeaderSampleMSB,
  kHeaderSampleLSB,
  // All kHeader... must be consecutive and follow the order they data format prescribes
  kChannelHeader,
  kADC,
  kChannelEnding,
  kUnknown
};

// Classify words into one of the types
WordType get_word_type( uint16_t word );

// Loop over a binary file, interpret words and write them to a ROOT file
int decoder( const char* argv );
