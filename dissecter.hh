#include "stdint.h"

// Types of XMIT words (32 bits)
enum XMITWordType{
  kXMITHeader,
  kXMITTailer,
  kXMITFEMWord,
};

// Types of FEM words (16 bits)
enum FEMWordType{
  kFEMHeaderFirst,
  kFEMHeaderIDSlot,
  kFEMHeaderNWordsMSB,
  kFEMHeaderNWordsLSB,
  kFEMHeaderEventMSB,
  kFEMHeaderEventLSB,
  kFEMHeaderFrameMSB,
  kFEMHeaderFrameLSB,
  kFEMHeaderChecksumMSB,
  kFEMHeaderChecksumLSB,
  kFEMHeaderSampleMSB,
  kFEMHeaderSampleLSB,
  // All kFEMHeader... must be consecutive and follow the order the data format prescribes
  kFEMChannelHeader,
  kFEMUnknown // To be used also as wildcard until waveform decoding is implemented
};

// Classify XMIT (32-bit) words into one of the types
XMITWordType get_xmit_word_type( uint32_t word );

// Classify FEM (16-bit) words into one of the types
FEMWordType get_fem_word_type( uint16_t word );

// Loop over a binary file, interpret words and write them to a ROOT file
int dissecter( const char* argv );
