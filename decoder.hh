
// Types of words
enum WordType{
  kHeaderFirst = 0,
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
  // All kHeader... must be consecutive and follow the order the data format prescribes
  kChannelHeader,
  kADC,
  kADCHuffman,
  kChannelEnding,
  kUnknown
};

// Classify words into one of the types
WordType get_word_type( uint16_t word );

// Decode Huffman code
int decode_huffman( int zeros );

// Loop over a binary file, interpret words and write them to a ROOT file
int decoder( const char* argv );
