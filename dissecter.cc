#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

#include "dissecter.hh"
#include "XMITInfo.hh"

bool verbose = false;
//bool verbose = true;

// Classify XMIT word
XMITWordType get_xmit_word_type( uint32_t word ){
  if( word == 0xFFFFFFFF ) return kXMITHeader;
  if( word == 0xE0000000 ) return kXMITTailer;
  else return kXMITFEMWord;
}

// Classify FEM word
FEMWordType get_fem_word_type( uint16_t word, size_t order ){
  static size_t headerCounter = 0; // Initialized only the first time
  // Extra checks on order and headerCounter are needed because in SN stream Huffman words can look like FEM headers
  if( (word == 0xFFFF) && (order == 0) && headerCounter == 0 ){ // The first FEM header appears always in the LSBs of the 32-bit word
    headerCounter = 1;
    return kFEMHeaderFirst;
  }
  if( ((word & 0xF000) == 0xF000)
      && headerCounter > 0 && headerCounter < 12 // Expect 6 16-bit header pairs
      && (headerCounter % 2) == order // Each header should be in the expected upper or lower 16 bits of the 32-bit word
      ){
    FEMWordType type = static_cast<FEMWordType>( kFEMHeaderFirst + headerCounter );
    headerCounter++;
    if( headerCounter == 12 ) headerCounter = 0;
    return type;
  }
  if( (word & 0xF000) == 0x1000 ){ // SN stream format
    headerCounter = 0;
    return kFEMChannelHeader;
  } else {
    headerCounter = 0;
    return kFEMUnknown; // For dissecter we are not interested in anything else
  }
}

// Loop over a binary file, interpret words and write them to a ROOT file
int dissecter( const char* argv ){

  std::string inFileName(argv);

  std::ifstream binFile;
  binFile.open( inFileName.c_str(), std::ios::binary );
  if( !binFile.is_open() ){
    std::cerr << "ERROR: Could not open file " << inFileName << std::endl;
    return 0;
  }

  std::string outFileName = inFileName.substr(0, inFileName.find_last_of(".")) + "_dissect_v2.root";
  TFile rootFile( outFileName.c_str(), "RECREATE" );

  // XMIT information
  XMITInfo xmitinfo;

  // FEM information
  FEMInfo feminfo;

  TTree* outTree = new TTree("dissecterTree", "Dissecter output tree");
  outTree->Branch("xmitinfo", &xmitinfo );

  int xmitpacket = 0;
  size_t femcount = 0;
  bool headerCompleted = false;
  FEMWordType prev_femword = kFEMUnknown;

  while( binFile.peek() != EOF ){
    uint32_t word32b;
    binFile.read( reinterpret_cast<char*>(&word32b), sizeof(word32b) );

    switch( get_xmit_word_type(word32b) ){

    case kXMITHeader:
      if( xmitpacket > 0 && xmitinfo.size() > 0 ){
	// The beginning of a new XMIT packet triggers the filling of the last FEM
	xmitinfo.push_back(feminfo);
	if( verbose ) std::cout << "FEM at slot " << (int)feminfo.slot << " processed" <<  std::endl;
	// The beginning of a new XMIT packet triggers the filling of the last one
	outTree->Fill();
	if( verbose ) std::cout << "XMIT packet " << xmitpacket << " written to TTree" <<  std::endl;
      }
      // Reset
      xmitinfo.clear();
      feminfo.clear_data();
      headerCompleted = false;
      femcount = 0;
      xmitpacket++;
      if( verbose ) std::cout << "Beginning to process XMIT packet " << xmitpacket << std::endl;
      break;

    case kXMITTailer:
      if( verbose ) std::cout << "Finishing to process XMIT packet " << xmitpacket << std::endl;
      break;

    case kXMITFEMWord:
      uint16_t first16b = word32b & 0xFFFF;
      uint16_t last16b = (word32b>>16) & 0xFFFF;
      uint16_t words16b[2] = {first16b, last16b};

      for(size_t i = 0; i < 2; i++){
	uint16_t word = words16b[i];
	std::string type;

	FEMWordType femword = get_fem_word_type(word, i);

	switch( femword ){

	case kFEMHeaderFirst:
	  type = "FEM Header First";
	  if( headerCompleted ){
	    // The beginning of a new FEM triggers the filling of the last one
	    xmitinfo.push_back(feminfo);
	    if( verbose ) std::cout << "FEM at slot " << (int)feminfo.slot << " processed" <<  std::endl;
	    // Reset
	    feminfo.clear_data();
	    headerCompleted = false;
	    femcount++;
	    if( verbose ) std::cout << "Beginning to process FEM " << femcount << std::endl;
	  }
	  prev_femword = femword;
	  break;

	case kFEMHeaderIDSlot:
	  type = "FEM Header ID+Slot";
	  if( prev_femword == kFEMHeaderFirst && !headerCompleted ){
	    feminfo.slot = (word & 0x1F);
	    feminfo.id = ((word>>5) & 0xF);
	    feminfo.test = ((word>>9) & 0x1);
	    feminfo.overflow = ((word>>10) & 0x1);
	    feminfo.full = ((word>>11) & 0x1);
	    prev_femword = femword;
	  }
	  break;

	case kFEMHeaderNWordsMSB:
	  if( prev_femword == kFEMHeaderIDSlot && !headerCompleted ){
	    type = "FEM Header N Words MSB";
	    feminfo.nwords += ((word & 0xFFF)<<12);
	    prev_femword = femword;
	  }
	  break;

	case kFEMHeaderNWordsLSB:
	  if( prev_femword == kFEMHeaderNWordsMSB && !headerCompleted ){
	    type = "FEM Header N Words LSB";
	    feminfo.nwords += (word & 0xFFF);
	    prev_femword = femword;
	  }
	  break;

	case kFEMHeaderEventMSB:
	  if( prev_femword == kFEMHeaderNWordsLSB && !headerCompleted ){
	    type = "FEM Header Event MSB";
	    feminfo.event += ((word & 0xFFF)<<12);
	    prev_femword = femword;
	  }
	  break;

	case kFEMHeaderEventLSB:
	  if( prev_femword == kFEMHeaderEventMSB && !headerCompleted ){
	    type = "FEM Header Event LSB";
	    feminfo.event += (word & 0xFFF);
	    prev_femword = femword;
	  }
	  break;

	case kFEMHeaderFrameMSB:
	  if( prev_femword == kFEMHeaderEventLSB && !headerCompleted ){
	    type = "FEM Header Frame MSB";
	    feminfo.frame += ((word & 0xFFF)<<12);
	    prev_femword = femword;
	  }
	  break;

	case kFEMHeaderFrameLSB:
	  if( prev_femword == kFEMHeaderFrameMSB && !headerCompleted ){
	    type = "FEM Header Frame LSB";
	    feminfo.frame += (word & 0xFFF);
	    prev_femword = femword;
	  }
	  break;

	case kFEMHeaderChecksumMSB:
	  if( prev_femword == kFEMHeaderFrameLSB && !headerCompleted ){
	    type = "FEM Header Checksum MSB";
	    feminfo.checksum += ((word & 0xFFF)<<12);
	    prev_femword = femword;
	  }
	  break;

	case kFEMHeaderChecksumLSB:
	  if( prev_femword == kFEMHeaderChecksumMSB && !headerCompleted ){
	    type = "FEM Header Checksum LSB";
	    feminfo.checksum += (word & 0xFFF);
	    prev_femword = femword;
	  }
	  break;

	case kFEMHeaderSampleMSB:
	  if( prev_femword == kFEMHeaderChecksumLSB && !headerCompleted ){
	    type = "FEM Header Sample MSB";
	    // Do nothing for the 6th header pair since it may be buggy/undefined
	    prev_femword = femword;
	  }
	  break;

	case kFEMHeaderSampleLSB:
	  if( prev_femword == kFEMHeaderSampleMSB && !headerCompleted ){
	    type = "FEM Header Sample LSB";
	    headerCompleted = true;
	    // Do nothing for the 6th header pair since it may be buggy/undefined
	    prev_femword = femword;
	  }
	  break;

	case kFEMChannelHeader:
	  type = "FEM Channel Header";
	  feminfo.channelcount++;
	  feminfo.wordcount++;
	  if( verbose ) std::cout << "Reading channel " << (word & 0x3F) << std::endl;
	  if( ((word >> 6) & 0x3F) != (feminfo.frame & 0x3F) ){
	    if( verbose ) std::cout << "WARNING: Frame LSBs from FEM Channel Header do not match the FEM Header Frame. Data is corrupt!" << std::endl;
	    feminfo.badframecount++;
	  }
	  break;

	case kFEMUnknown:
	  type = "FEM Unknown";
	  feminfo.wordcount++;
	  break;
	} // end of FEM switch
	//if( verbose ) std::cout << std::setfill('0');
	//if( verbose ) std::cout << std::hex << std::setw(4) << words16b[i] << " is " << type <<  std::endl;
	//if( verbose ) std::cout << std::dec; // revert to decimal
      } // end of loop over 2 words
      break;
    } // end of XMIT switch
  } // end of reading the file
  // Write the last FEM and XMIT packet (might be incomplete)
  xmitinfo.push_back(feminfo);
  outTree->Fill();

  binFile.close();
  outTree->Write();
  rootFile.Close();
  return 0;
}

// To run as a standalone application
# ifndef __CINT__
int main( int argc, char** argv ){
  return dissecter( argv[1] );
}
# endif
