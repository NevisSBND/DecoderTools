#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <deque>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

#include "dissecter.hh"
#include "XMITInfo.hh"

//bool verbose = false;
bool verbose = true;
bool print_matrix = true; // Print summary matrix
bool print_labels = true; // Print summary matrix labels

// Classify 32-bit word (XMIT or 2 FEMs)
NevisWordType get_32bit_word_type( uint32_t word ){
  if( word == 0xFFFFFFFF ) return kXMITHeader;
  if( word == 0xE0000000 ) return kXMITTailer;
  else return k2FEMWords;
}

// Classify FEM word
NevisWordType get_fem_word_type( uint16_t word, size_t order ){
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
    NevisWordType type = static_cast<NevisWordType>( kFEMHeaderFirst + headerCounter );
    headerCounter++;
    if( headerCounter == 12 ) headerCounter = 0;
    return type;
  }
  if( (word & 0xF000) == 0x1000 ){ // SN stream format
    headerCounter = 0;
    return kFEMChannelHeader;
  } else if( (word & 0xC000) == 0x4000 ){
    headerCounter = 0;
    return kFEMTimeHeader;
  } else if( (word & 0xF000) == 0x2000 ){
    headerCounter = 0;
    return kFEMADC;
  } else if( (word & 0x8000) == 0x8000 ){
    headerCounter = 0;
    return kFEMADCHuffman;
  } else if( (word & 0xF000) == 0x3000 ){
    headerCounter = 0;
    return kFEMADCTailer;
  } else {
    headerCounter = 0;
    return kFEMUnknown;
  }
}

// Convert Nevis type into string
std::string word_type_to_string( NevisWordType type ){
  std::string output;
  switch( type ){
  case kXMITHeader:
    output = "XMIT Header";
    break;
  case kXMITTailer:
    output = "XMIT Tailer";
    break;
  case kFEMHeaderFirst:
    output = "FEM Header First";
    break;
  case kFEMHeaderIDSlot:
    output = "FEM Header ID+Slot";
    break;
  case kFEMHeaderNWordsMSB:
    output = "FEM Header N Words MSB";
    break;
  case kFEMHeaderNWordsLSB:
    output = "FEM Header N Words LSB";
    break;
  case kFEMHeaderEventMSB:
    output = "FEM Header Event MSB";
    break;
  case kFEMHeaderEventLSB:
    output = "FEM Header Event LSB";
    break;
  case kFEMHeaderFrameMSB:
    output = "FEM Header Frame MSB";
    break;
  case kFEMHeaderFrameLSB:
    output = "FEM Header Frame LSB";
    break;
  case kFEMHeaderChecksumMSB:
    output = "FEM Header Checksum MSB";
    break;
  case kFEMHeaderChecksumLSB:
    output = "FEM Header Checksum LSB";
    break;
  case kFEMHeaderSampleMSB:
    output = "FEM Header Sample MSB";
    break;
  case kFEMHeaderSampleLSB:
    output = "FEM Header Sample LSB";
    break;
  case kFEMChannelHeader:
    output = "FEM Channel Header";
    break;
  case kFEMTimeHeader:
    output = "FEM Time Header";
    break;
  case kFEMADC:
    output = "FEM ADC";
    break;
  case kFEMADCHuffman:
    output = "FEM ADC Huffman";
    break;
  case kFEMADCTailer:
    output = "FEM ADC Tailer";
    break;
  default:
    output = "FEM Unknown";
    break;
  }
  return output;
}

// Add new 16-bit word to FIFO and remove oldest one
void add_to_fifo( std::deque<uint16_t>& fifo, uint16_t word ){
  fifo.push_back(word);
  fifo.pop_front();
}

// Dump content of FIFO
void dump_fifo( std::deque<uint16_t>& fifo ){
  if(verbose){
    for(size_t w = 0; w < fifo.size(); w++){
      std::cout << std::uppercase << std::setfill('0') << std::setw(4) << std::hex << fifo[w] << " ";
      if( (w + 1) % 16 == 0 ) std::cout << std::dec << "\n";
    }
    std::cout << std::dec << std::setfill(' ') << std::endl;
  }
}

// Return adjacency matrix element
bool is_adjacent( NevisWordType row, NevisWordType col ){
  if( row == kXMITHeader ){
    if( col != kXMITTailer ) return false;
    else return true;
  } else if( row == kXMITTailer ){
    if( (col != kFEMChannelHeader) && (col != kFEMADCTailer) && (col != kFEMADC) && (col != kFEMADCHuffman) ) return false;
    else return true;
  } else if( row == kFEMHeaderFirst ){
    if( (col != kFEMChannelHeader) && (col != kXMITHeader) && (col != kFEMADCTailer) && (col != kFEMADC) && (col != kFEMADCHuffman) ) return false;
    else return true;
  } else if( row == kFEMHeaderIDSlot ){
    if( col != kFEMHeaderFirst ) return false;
    else return true;
  } else if( row == kFEMHeaderNWordsMSB ){
    if( col != kFEMHeaderIDSlot ) return false;
    else return true;
  } else if( row == kFEMHeaderNWordsLSB ){
    if( col != kFEMHeaderNWordsMSB ) return false;
    else return true;
  } else if( row == kFEMHeaderEventMSB ){
    if( col != kFEMHeaderNWordsLSB ) return false;
    else return true;
  } else if( row == kFEMHeaderEventLSB ){
    if( col != kFEMHeaderEventMSB ) return false;
    else return true;
  } else if( row == kFEMHeaderFrameMSB ){
    if( col != kFEMHeaderEventLSB ) return false;
    else return true;
  } else if( row == kFEMHeaderFrameLSB ){
    if( col != kFEMHeaderFrameMSB ) return false;
    else return true;
  } else if( row == kFEMHeaderChecksumMSB ){
    if( col != kFEMHeaderFrameLSB ) return false;
    else return true;
  } else if( row == kFEMHeaderChecksumLSB ){
    if( col != kFEMHeaderChecksumMSB ) return false;
    else return true;
  } else if( row == kFEMHeaderSampleMSB ){
    if( col != kFEMHeaderChecksumLSB ) return false;
    else return true;
  } else if( row == kFEMHeaderSampleLSB ){
    if( col != kFEMHeaderSampleMSB ) return false;
    else return true;
  } else if( row == kFEMChannelHeader ){
    if( (col != kFEMChannelHeader) && (col != kFEMADCTailer) && (col != kFEMHeaderSampleLSB) ) return false;
    else return true;
  } else if( row == kFEMTimeHeader ){
    if( (col != kFEMChannelHeader) && (col != kFEMADCTailer) ) return false;
    else return true;
  } else if( row == kFEMADC ){
    if( (col != kFEMADC) && (col != kFEMADCHuffman) && (col != kFEMTimeHeader) && (col != kFEMChannelHeader) ) return false;
    else return true;
  } else if( row == kFEMADCHuffman ){
    if( (col != kFEMADC) && (col != kFEMADCHuffman) ) return false;
    else return true;
  } else if( row == kFEMADCTailer ){
    if( (col != kFEMADCHuffman) && (col != kFEMADC) ) return false;
    else return true;
  } else { // Default case
    return false;
  }
}

// Color in red non-adjacent elements
std::string color_element( NevisWordType row, NevisWordType col, int val ){
  if( is_adjacent(row, col) ) return blacktxt;
  else return ((val == 0) ? blacktxt : redtxt);
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

  std::string outFileName = inFileName.substr(0, inFileName.find_last_of(".")) + "_dissect_v3.root";
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
  NevisWordType prev_word = kFEMUnknown;

  std::vector< std::vector<int> > graphMatrix; // Summary matrix for the graph describing the Nevis Words relations
  graphMatrix.resize(21); // There are 21 Nevis words
  for(size_t i = 0; i < graphMatrix.size(); i++){
    graphMatrix[i].resize(21);
  }

  std::deque<uint16_t> fifo (32, 0xBAD);

  // Process binary file
  while( binFile.peek() != EOF ){
    uint32_t word32b;
    binFile.read( reinterpret_cast<char*>(&word32b), sizeof(word32b) );

    NevisWordType word_type_32b = get_32bit_word_type(word32b);

    switch( word_type_32b ){

    case kXMITHeader:
      add_to_fifo( fifo, word32b & 0xFFFF );
      add_to_fifo( fifo, (word32b>>16) & 0xFFFF );
      if( xmitpacket > 0 && xmitinfo.size() > 0 ){
	// The beginning of a new XMIT packet triggers the filling of the last FEM
	feminfo.mychecksum = (feminfo.mychecksum & 0xFFFFFF);
	uint32_t checksum_diff = feminfo.checksum - feminfo.mychecksum;
	if( verbose ) if( checksum_diff != 0 || checksum_diff != 28762 ) std::cout << "Checksum difference: " << checksum_diff << std::endl;
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

      if( prev_word != kXMITTailer ){ // To do: the function is_adjacent could replace this check
	if( verbose ) std::cout << "Format violation: current word " << word_type_to_string( word_type_32b ) << " preceded by " << word_type_to_string( prev_word ) << std::endl;
	dump_fifo(fifo);
      }

      graphMatrix[word_type_32b][prev_word]++;
      prev_word = word_type_32b;
      break;

    case kXMITTailer:
      add_to_fifo( fifo, word32b & 0xFFFF );
      add_to_fifo( fifo, (word32b>>16) & 0xFFFF );
      if( verbose ) std::cout << "Finishing to process XMIT packet " << xmitpacket << std::endl;

      if( (prev_word != kFEMChannelHeader) && (prev_word != kFEMADCTailer) ){
	if( verbose ){
	  std::cout << "Format violation: current word " << word_type_to_string( word_type_32b ) << " preceded by " << word_type_to_string( prev_word ) << std::endl;
	  if( (prev_word == kFEMADC) || (prev_word == kFEMADCHuffman) ) std::cout << "Lost FEM ADC Tailer" << std::endl;
	}
	dump_fifo(fifo);
      }

      graphMatrix[word_type_32b][prev_word]++;
      prev_word = word_type_32b;
      break;

    default: //case k2FEMWords:
      uint16_t first16b = word32b & 0xFFFF;
      uint16_t last16b = (word32b>>16) & 0xFFFF;
      uint16_t words16b[2] = {first16b, last16b};

      for(size_t i = 0; i < 2; i++){
	uint16_t word = words16b[i];
	add_to_fifo( fifo, word );

	NevisWordType femword = get_fem_word_type(word, i);

	switch( femword ){

	case kFEMHeaderFirst:
	  if( headerCompleted ){
	    // The beginning of a new FEM triggers the filling of the last one
	    feminfo.mychecksum = (feminfo.mychecksum & 0xFFFFFF);
	    uint32_t checksum_diff = feminfo.checksum - feminfo.mychecksum;
	    if( verbose ) if( checksum_diff != 0 || checksum_diff != 28762 ) std::cout << "Checksum difference: " << checksum_diff << std::endl;
	    xmitinfo.push_back(feminfo);
	    if( verbose ) std::cout << "FEM at slot " << (int)feminfo.slot << " processed" <<  std::endl;
	    // Reset
	    feminfo.clear_data();
	    headerCompleted = false;
	    femcount++;
	    if( verbose ) std::cout << "Beginning to process FEM " << femcount << std::endl;
	  }

	  if( (prev_word != kFEMChannelHeader) && (prev_word != kXMITHeader) && (prev_word != kFEMADCTailer) ){
	    if( verbose ){
	      std::cout << "Format violation: current word " << word_type_to_string( femword ) << " preceded by " << word_type_to_string( prev_word ) << std::endl;
	      if( (prev_word == kFEMADC) || (prev_word == kFEMADCHuffman) ) std::cout << "\tLost FEM ADC Tailer?" << std::endl;
	    }
	    dump_fifo(fifo);
	  }

	  graphMatrix[femword][prev_word]++;
	  prev_word = femword;
	  break;

	case kFEMHeaderIDSlot:
	  if( prev_word == kFEMHeaderFirst && !headerCompleted ){
	    feminfo.slot = (word & 0x1F);
	    feminfo.id = ((word>>5) & 0xF);
	    feminfo.test = ((word>>9) & 0x1);
	    feminfo.overflow = ((word>>10) & 0x1);
	    feminfo.full = ((word>>11) & 0x1);

	    if( prev_word != kFEMHeaderFirst ){
	      if( verbose ) std::cout << "Format violation: current word " << word_type_to_string( femword )
				      << " preceded by " << word_type_to_string( prev_word ) << std::endl;
	      dump_fifo(fifo);
	    }

	    graphMatrix[femword][prev_word]++;
	    prev_word = femword;
	  }
	  break;

	case kFEMHeaderNWordsMSB:
	  if( prev_word == kFEMHeaderIDSlot && !headerCompleted ){
	    feminfo.nwords += ((word & 0xFFF)<<12);

	    if( prev_word != kFEMHeaderIDSlot ){
	      if( verbose ) std::cout << "Format violation: current word " << word_type_to_string( femword )
				      << " preceded by " << word_type_to_string( prev_word ) << std::endl;
	      dump_fifo(fifo);
	    }

	    graphMatrix[femword][prev_word]++;
	    prev_word = femword;
	  }
	  break;

	case kFEMHeaderNWordsLSB:
	  if( prev_word == kFEMHeaderNWordsMSB && !headerCompleted ){
	    feminfo.nwords += (word & 0xFFF);

	    if( prev_word != kFEMHeaderNWordsMSB ){
	      if( verbose ) std::cout << "Format violation: current word " << word_type_to_string( femword )
				      << " preceded by " << word_type_to_string( prev_word ) << std::endl;
	      dump_fifo(fifo);
	    }

	    graphMatrix[femword][prev_word]++;
	    prev_word = femword;
	  }
	  break;

	case kFEMHeaderEventMSB:
	  if( prev_word == kFEMHeaderNWordsLSB && !headerCompleted ){
	    feminfo.event += ((word & 0xFFF)<<12);

	    if( prev_word != kFEMHeaderNWordsLSB ){
	      if( verbose ) std::cout << "Format violation: current word " << word_type_to_string( femword )
				      << " preceded by " << word_type_to_string( prev_word ) << std::endl;
	      dump_fifo(fifo);
	    }

	    graphMatrix[femword][prev_word]++;
	    prev_word = femword;
	  }
	  break;

	case kFEMHeaderEventLSB:
	  if( prev_word == kFEMHeaderEventMSB && !headerCompleted ){
	    feminfo.event += (word & 0xFFF);

	    if( prev_word != kFEMHeaderEventMSB ){
	      if( verbose ) std::cout << "Format violation: current word " << word_type_to_string( femword )
				      << " preceded by " << word_type_to_string( prev_word ) << std::endl;
	      dump_fifo(fifo);
	    }

	    graphMatrix[femword][prev_word]++;
	    prev_word = femword;
	  }
	  break;

	case kFEMHeaderFrameMSB:
	  if( prev_word == kFEMHeaderEventLSB && !headerCompleted ){
	    feminfo.frame += ((word & 0xFFF)<<12);

	    if( prev_word != kFEMHeaderEventLSB ){
	      if( verbose ) std::cout << "Format violation: current word " << word_type_to_string( femword )
				      << " preceded by " << word_type_to_string( prev_word ) << std::endl;
	      dump_fifo(fifo);
	    }

	    graphMatrix[femword][prev_word]++;
	    prev_word = femword;
	  }
	  break;

	case kFEMHeaderFrameLSB:
	  if( prev_word == kFEMHeaderFrameMSB && !headerCompleted ){
	    feminfo.frame += (word & 0xFFF);

	    if( prev_word != kFEMHeaderFrameMSB ){
	      if( verbose ) std::cout << "Format violation: current word " << word_type_to_string( femword )
				      << " preceded by " << word_type_to_string( prev_word ) << std::endl;
	      dump_fifo(fifo);
	    }

	    graphMatrix[femword][prev_word]++;
	    prev_word = femword;
	  }
	  break;

	case kFEMHeaderChecksumMSB:
	  if( prev_word == kFEMHeaderFrameLSB && !headerCompleted ){
	    feminfo.checksum += ((word & 0xFFF)<<12);

	    if( prev_word != kFEMHeaderFrameLSB ){
	      if( verbose ) std::cout << "Format violation: current word " << word_type_to_string( femword )
				      << " preceded by " << word_type_to_string( prev_word ) << std::endl;
	      dump_fifo(fifo);
	    }

	    graphMatrix[femword][prev_word]++;
	    prev_word = femword;
	  }
	  break;

	case kFEMHeaderChecksumLSB:
	  if( prev_word == kFEMHeaderChecksumMSB && !headerCompleted ){
	    feminfo.checksum += (word & 0xFFF);

	    if( prev_word != kFEMHeaderChecksumMSB ){
	      if( verbose ) std::cout << "Format violation: current word " << word_type_to_string( femword )
				      << " preceded by " << word_type_to_string( prev_word ) << std::endl;
	      dump_fifo(fifo);
	    }

	    graphMatrix[femword][prev_word]++;
	    prev_word = femword;
	  }
	  break;

	case kFEMHeaderSampleMSB:
	  if( prev_word == kFEMHeaderChecksumLSB && !headerCompleted ){
	    // Do nothing for the 6th header pair since it may be buggy/undefined

	    if( prev_word != kFEMHeaderChecksumLSB ){
	      if( verbose ) std::cout << "Format violation: current word " << word_type_to_string( femword )
				      << " preceded by " << word_type_to_string( prev_word ) << std::endl;
	      dump_fifo(fifo);
	    }

	    graphMatrix[femword][prev_word]++;
	    prev_word = femword;
	  }
	  break;

	case kFEMHeaderSampleLSB:
	  if( prev_word == kFEMHeaderSampleMSB && !headerCompleted ){
	    headerCompleted = true;
	    // Do nothing for the 6th header pair since it may be buggy/undefined

	    if( prev_word != kFEMHeaderSampleMSB ){
	      if( verbose ) std::cout << "Format violation: current word " << word_type_to_string( femword )
				      << " preceded by " << word_type_to_string( prev_word ) << std::endl;
	      dump_fifo(fifo);
	    }

	    graphMatrix[femword][prev_word]++;
	    prev_word = femword;
	  }
	  break;

	case kFEMChannelHeader:
	  feminfo.channelcount++;
	  feminfo.wordcount++;
	  feminfo.mychecksum += word;
	  if( verbose ) std::cout << "Reading FEM@slot " << (int)feminfo.slot << " channel " << (word & 0x3F) << std::endl;
	  if( ((word >> 6) & 0x3F) != (feminfo.frame & 0x3F) ){
	    if( verbose ) std::cout << "WARNING: Frame LSBs from FEM Channel Header do not match the FEM Header Frame. Data is corrupt!" << std::endl;
	    feminfo.badframecount++;
	  }

	  if( (prev_word != kFEMChannelHeader) && (prev_word != kFEMADCTailer) && (prev_word != kFEMHeaderSampleLSB) ){
	    if( verbose ) std::cout << "Format violation: current word " << word_type_to_string( femword )
				    << " preceded by " << word_type_to_string( prev_word ) << std::endl;
	    dump_fifo(fifo);
	  }

	  graphMatrix[femword][prev_word]++;
	  prev_word = femword;
	  break;

	case kFEMTimeHeader:
	  feminfo.wordcount++;
	  feminfo.mychecksum += word;

	  if( (prev_word != kFEMChannelHeader) && (prev_word != kFEMADCTailer) ){
	    if( verbose ) std::cout << "Format violation: current word " << word_type_to_string( femword )
				    << " preceded by " << word_type_to_string( prev_word ) << std::endl;
	    dump_fifo(fifo);
	  }

	  graphMatrix[femword][prev_word]++;
	  prev_word = femword;
	  break;

	case kFEMADC:
	  feminfo.wordcount++;
	  feminfo.mychecksum += word;

	  if( (prev_word != kFEMADC) && (prev_word != kFEMADCHuffman) && (prev_word != kFEMTimeHeader) ){
	    if( verbose ){
	      std::cout << "Format violation: current word " << word_type_to_string( femword )
			<< " preceded by " << word_type_to_string( prev_word ) << std::endl;
	      if( prev_word == kFEMChannelHeader ) std::cout << "\tRollover?" << std::endl;
	    }
	    dump_fifo(fifo);
	  }

	  graphMatrix[femword][prev_word]++;
	  prev_word = femword;
	  break;

	case kFEMADCHuffman:
	  feminfo.wordcount++;
	  feminfo.mychecksum += word;

	  if( (prev_word != kFEMADC) && (prev_word != kFEMADCHuffman) ){
	    if( verbose ) std::cout << "Format violation: current word " << word_type_to_string( femword )
				    << " preceded by " << word_type_to_string( prev_word ) << std::endl;
	    dump_fifo(fifo);
	  }

	  graphMatrix[femword][prev_word]++;
	  prev_word = femword;
	  break;

	case kFEMADCTailer:
	  feminfo.wordcount++;
	  feminfo.mychecksum += word;

	  if( (prev_word != kFEMADCHuffman) && (prev_word != kFEMADC) ){
	    if( verbose ) std::cout << "Format violation: current word " << word_type_to_string( femword )
				    << " preceded by " << word_type_to_string( prev_word ) << std::endl;
	    dump_fifo(fifo);
	  }

	  graphMatrix[femword][prev_word]++;
	  prev_word = femword;
	  break;

	default: //case kFEMUnknown:
	  feminfo.wordcount++;
	  feminfo.mychecksum += word;
	  graphMatrix[femword][prev_word]++;
	  prev_word = femword;
	  break;
	} // end of FEM switch
	//if( verbose ) std::cout << std::setfill('0');
	//if( verbose ) std::cout << std::hex << std::setw(4) << words16b[i] << " is " << word_type_to_string( femword ) <<  std::endl;
	//if( verbose ) std::cout << std::dec; // revert to decimal
      } // end of loop over 2 words
      break;
    } // end of XMIT switch
  } // end of reading the file
  // Write the last FEM and XMIT packet (might be incomplete)
  feminfo.mychecksum = (feminfo.mychecksum & 0xFFFFFF);
  uint32_t checksum_diff = feminfo.checksum - feminfo.mychecksum;
  if( verbose ) if( checksum_diff != 0 || checksum_diff != 28762 ) std::cout << "Checksum difference: " << checksum_diff << std::endl;
  xmitinfo.push_back(feminfo);
  outTree->Fill();

  if( print_matrix ){
    // Print summary matrix
    float total_words;
    for(size_t i = 0; i < graphMatrix.size(); i++) for( size_t j = 0; j < graphMatrix[i].size(); j++) total_words += graphMatrix[i][j];
    std::cout << "*** Summary matrix ***\n";
    for(size_t i = 0; i < graphMatrix.size(); i++){
      // Print column labels
      if( print_labels ){
	std::cout << std::setw(23 + 3) << " "; // Margin for row labels
	for( size_t j = 0; j < graphMatrix[i].size(); j++ ){
	  std::cout << " " << std::setw(2) << j;
	}
	std::cout << "\n";
      }
      // Print row label
      std::cout << std::setw(23) << word_type_to_string( static_cast<NevisWordType>(i) );
      if( print_labels ) std::cout << " " << std::setw(2) << i;
      float good_words = 0, bad_words = 0, row_words = 0;
      // Print matrix row
      for( size_t j = 0; j < graphMatrix[i].size(); j++){
	std::cout << " " << color_element( static_cast<NevisWordType>(i), static_cast<NevisWordType>(j), graphMatrix[i][j] ) << std::setw(2) << graphMatrix[i][j];
	(is_adjacent( static_cast<NevisWordType>(i), static_cast<NevisWordType>(j) ) || (graphMatrix[i][j] == 0)) ? good_words += graphMatrix[i][j] : bad_words += graphMatrix[i][j];
	row_words += graphMatrix[i][j];
      }
      std::cout << blacktxt << " Good (row): " << 100.*good_words/row_words << "% Bad (row): " << 100.*bad_words/row_words
		<< "% Good (tot): " << 100.*good_words/total_words << "% Bad (tot): " << 100.*bad_words/total_words << "%\n";
    }
    std::cout << std::endl;
  }

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
