#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

#include "decoder.hh"
#include "HeaderInfo.hh"

// Classify word
WordType get_word_type( uint16_t word ){
  static int headerCounter = 0; // Initialized only the first time
  if( ((word & 0xF000) == 0xF000)
      || headerCounter < 12 // Hack: since 6th header pair does not start with "F". Bug?
      ){
    headerCounter++;
    return static_cast<WordType>( kHeaderFirst + headerCounter - 1 );
  }
  if( (word & 0xF000) == 0x4000 ) return kChannelHeader;
  if( (word & 0xF000) == 0x0000 ) return kADC;
  if( (word & 0xC000) == 0x8000 ) return kADCHuffman; // Look for headers first
  if( (word & 0xF000) == 0x5000 ){
    if( (word & 0x3F) == 63 ) headerCounter = 0; // Hack: look for header only after channel 63 is finished
    return kChannelEnding;
  }
  else return kUnknown;
}

// Decode Huffman
int decode_huffman( int zeros ){
  switch( zeros ){
  case 0: return 0;
  case 1: return -1;
  case 2: return 1;
  case 3: return -2;
  case 4: return 2;
  case 5: return -3;
  case 6: return 3;
  default:
    std::cerr << "ERROR: Number of zeros (" << zeros << ") in Huffman word is out of bounds" << std::endl;
    return 4096;
  }
}

// Loop over a binary file, interpret words and write them to a ROOT file
//int main( int argc, char** argv ){
int decoder( const char* argv ){
  //  std::string inFileName = argv[1];
  std::string inFileName(argv);

  std::ifstream binFile;
  binFile.open( inFileName.c_str(), std::ios::binary );
  if( !binFile.is_open() ){
    std::cerr << "ERROR: Could not open file " << inFileName << std::endl;
    return 0;
  }

  std::string outFileName = inFileName.substr(0, inFileName.find_last_of(".")) + ".root";
  TFile rootFile( outFileName.c_str(), "RECREATE" );

  // Header information
  HeaderInfo header;
  // Matrix of waveforms: 64 channels x N samples
  std::vector< std::vector<uint16_t> > waveform;
  waveform.resize(64);
  size_t currentChannel;
  std::vector<uint16_t> currentWaveform;

  TTree* outTree = new TTree("decoderTree", "Decoder output tree");
  outTree->Branch("header", &header );
  outTree->Branch("waveform", &waveform );

  int entry = 0;
  while( binFile.peek() != EOF ){
    uint32_t word32b;
    binFile.read( reinterpret_cast<char*>(&word32b), sizeof(word32b) );
    if( (word32b == 0xFFFFFFFF) || (word32b == 0xE0000000) ){ // Temporary: ignore XMIT words
      std::cout << std::setfill('0');
      std::cout << "INFO: XMIT word " << std::hex << std::setw(8) << word32b << " found and ignored" <<  std::endl;
      std::cout << std::dec; // revert to decimal
      continue;
    }
    uint16_t first16b = word32b & 0xFFFF;
    uint16_t last16b = (word32b>>16) & 0xFFFF;
    uint16_t words16b[2] = {first16b, last16b};

    for(size_t i = 0; i < 2; i++){
      uint16_t word = words16b[i];
      std::string type;
      switch( get_word_type( word ) ){
      case kHeaderFirst:
	type = "Header First";
	if(entry > 0){
	  // The beginning of a new entry triggers the filling of the last one
	  header.wordcount = (header.wordcount & 0xFFFFFF);
	  int32_t nwords_diff = header.nwords - header.wordcount;
	  if( nwords_diff != 0 ) std::cerr << "WARNING: Word count difference: " << nwords_diff << std::endl;
	  header.mychecksum = (header.mychecksum & 0xFFFFFF);
	  uint32_t checksum_diff = header.checksum - header.mychecksum;
	  if( checksum_diff != 0 ) std::cerr << "WARNING: Checksum difference: " << checksum_diff << std::endl;
	  outTree->Fill();
	  std::cout << "Entry " << entry << " written to TTree" <<  std::endl;
	}
	// Reset
	header.clear();
	for(size_t ch = 0; ch < 64; ch++){
	  waveform[ch].clear();
	}
	entry++;
	std::cout << "Beginning to process entry " << entry << std::endl;
	break;
      case kHeaderIDSlot:
	type = "Header ID+Slot";
	header.slot = (word & 0x1F);
	header.id = ((word>>5) & 0xF);
	header.test = ((word>>9) & 0x1);
	header.overflow = ((word>>10) & 0x1);
	header.full = ((word>>11) & 0x1);
	break;
      case kHeaderNWordsMSB:
	type = "Header N Words MSB";
	header.nwords += ((word & 0xFFF)<<12);
	break;
      case kHeaderNWordsLSB:
	type = "Header N Words LSB";
	header.nwords += (word & 0xFFF);
	break;
      case kHeaderEventMSB:
	type = "Header Event MSB";
	header.event += ((word & 0xFFF)<<12);
	break;
      case kHeaderEventLSB:
	type = "Header Event LSB";
	header.event += (word & 0xFFF);
	break;
      case kHeaderFrameMSB:
	type = "Header Frame MSB";
	header.frame += ((word & 0xFFF)<<12);
	break;
      case kHeaderFrameLSB:
	type = "Header Frame LSB";
	header.frame += (word & 0xFFF);
	break;
      case kHeaderChecksumMSB:
	type = "Header Checksum MSB";
	header.checksum += ((word & 0xFFF)<<12);
	break;
      case kHeaderChecksumLSB:
	type = "Header Checksum LSB";
	header.checksum += (word & 0xFFF);
	break;
      case kHeaderSampleMSB:
	type = "Header Sample MSB";
	// Do nothing for the 6th header pair since it may be buggy/undefined
	break;
      case kHeaderSampleLSB:
	type = "Header Sample LSB";
	// Do nothing for the 6th header pair since it may be buggy/undefined
	break;
      case kChannelHeader:
	type = "Channel Header";
	header.wordcount++;
	header.mychecksum += word;
	currentChannel = (word & 0x3F);
	currentWaveform.clear();
	std::cout << "Reading channel " << currentChannel << std::endl;
	break;
      case kADC:
	type = "ADC";
	header.wordcount++;
	header.mychecksum += word;
	currentWaveform.push_back( (word & 0xFFF) );
	//std::cout << "ADC value " << (word & 0xFFF) << std::endl;
	break;
      case kADCHuffman: { // Brackets needed to declare variables
	type = "ADC Huffman";
	header.wordcount++;
	header.mychecksum += word;
	int zeros = 0; // Counter of 0s interleaved between 1s
	std::vector<int> differences; // Vector with Huffman-decoded differences
	differences.reserve(14); // Speed up: reserve for the maximum number of differences stored in a Huffman word
	// Read the lowest 14 bits from left to right
	for(uint16_t mask = 0x2000; mask > 0x0; mask = (mask>>1)){
	  if( (word & mask) == mask ){ // Found 1
	    differences.push_back( decode_huffman(zeros) );
	    zeros = 0; // Reset counter
	  } else zeros++;
	}
	// Differences are time-ordered from right to left
	for(size_t i = differences.size(); i > 0; i--){
	  currentWaveform.push_back( currentWaveform.back() + differences[i - 1] );
	}
	differences.clear();
      }	break;
      case kChannelEnding:
	type = "Channel Ending";
	header.wordcount++;
	header.mychecksum += word;
	if( (word & 0x3F) == currentChannel ){
	  waveform[currentChannel] = currentWaveform;
	  std::cout << "Finished reading channel " << currentChannel << std::endl;
	}
	else std::cerr << "ERROR: Channel header " << currentChannel << " and ending " << (word & 0x3F) << " do not match" << std::endl;
	// Reset
	currentChannel = 999;
	currentWaveform.clear();
	break;
      case kUnknown:
	header.wordcount++;
	header.mychecksum += word;
	type = "Unknown";
	break;
      } // end of switch
      //std::cout << std::setfill('0');
      //std::cout << std::hex << std::setw(4) << words16b[i] << " is " << type <<  std::endl;
      //std::cout << std::dec; // revert to decimal
    } // end of loop over 2 words
  } // end of reading the file
  // Write the last entry (might be incomplete)
  header.wordcount = (header.wordcount & 0xFFFFFF);
  int32_t nwords_diff = header.nwords - header.wordcount;
  if( nwords_diff != 0 ) std::cerr << "WARNING: Word count difference: " << nwords_diff << std::endl;
  header.mychecksum = (header.mychecksum & 0xFFFFFF);
  uint32_t checksum_diff = header.checksum - header.mychecksum;
  if( checksum_diff != 0 ) std::cerr << "WARNING: Checksum difference: " << checksum_diff << std::endl;
  outTree->Fill();

  binFile.close();
  outTree->Write();
  rootFile.Close();
  return 1;
}

// To run as a standalone application
# ifndef __CINT__
int main( int argc, char** argv ){
  return decoder( argv[1] );
}
# endif
