#include "stdint.h"
#include <deque>

const std::string redtxt("\033[1;31m");
const std::string blacktxt("\033[0m");

// Types of Nevis readout words
// Both 16-bits and 32-bits ones in a single enum to use them as states for a state machine
enum NevisWordType{
  // Types of FEM words (16 bits)
  kFEMHeaderFirst = 0,
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
  kFEMTimeHeader,
  kFEMADC,
  kFEMADCHuffman,
  kFEMADCTailer,
  kFEMUnknown,

  // Types of 32-bit words
  kXMITHeader,
  kXMITTailer,
  k2FEMWords, // Generic type for two 16-bit FEM words concatenated
};

// Classify XMIT (32-bit) words into one of the Nevis types
NevisWordType get_xmit_word_type( uint32_t word );

// Classify FEM (16-bit) words into one of the Nevis types
NevisWordType get_fem_word_type( uint16_t word );

// Convert Nevis type into string
std::string word_type_to_string( NevisWordType type );

// Loop over a binary file, interpret words and write them to a ROOT file
int dissecter( const char* argv );

// Add new 16-bit word to FIFO and remove oldest one
void add_to_fifo( std::deque<uint16_t>& fifo, uint16_t word );

// Dump content of FIFO
void dump_fifo( std::deque<uint16_t>& fifo );

// Return adjacency matrix element
bool is_adjacent( NevisWordType row, NevisWordType col );

// Color in red non-adjacent elements
std::string color_element( NevisWordType row, NevisWordType col );

//Decode Huffman code
int decode_huffman(int zeros);

//Calculate the average difference between neighboring samples
double FlippingBitMetric( std::vector<uint16_t> ROIs);

struct channel_ROI_info;
