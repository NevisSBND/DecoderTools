#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>

#include <TRint.h>
#include <TFile.h>
#include <TTree.h>

#include "HeaderInfo.hh"

// Struct holding all values that define a channel
struct channel_address{

  int fWire; // TPC wire within plane
  std::string fPlane; // Plane withing TPC
  std::string fAdapterConnector; // BNL adapter board connector
  int fAdapterPin; // BNL adapter board pin
  int fAnalogConnector; // BNL analog board connector
  int fAnalogPin; // BNL analog board pin
  int fASIC; // BNL ASIC number
  int fASICch; // BNL ASIC channel
  int fFEMBch; // BNL FEM channel
  int fFEMB; // BNL FEMB number
  int fWIB; // WIB number
  int fSlot; // Nevis FEM slot within crate
  int fFEMch; // Channel within Nevis FEM

  channel_address() 
    : fWire(-1), fPlane("NotUsed"), fAdapterConnector("N/A"), fAdapterPin(-1), fAnalogConnector(-1), fAnalogPin(-1), 
      fASIC(-1), fASICch(-1), fFEMBch(-1), fFEMB(-1), fWIB(-1),
      fSlot(-1), fFEMch(-1) {}

  channel_address( int newWire, std::string newPlane, std::string newAdapterConnector, int newAdapterPin, int newAnalogConnector, int newAnalogPin,
		   int newAsic, int newAsicch, int newFembch, int newFemb, int newWib,
		   int newSlot = -1, int newFemch = -1 )
    : fWire(newWire), fPlane(newPlane), fAdapterConnector(newAdapterConnector), fAdapterPin(newAdapterPin), 
      fAnalogConnector(newAnalogConnector), fAnalogPin(newAnalogPin),
      fASIC(newAsic), fASICch(newAsicch), fFEMBch(newFembch), fFEMB(newFemb), fWIB(newWib),
      fSlot(newSlot), fFEMch(newFemch) {}
};

int channel_mapper( const char* mapper_run, const char* bnl_map_filename ){

  std::map< int, channel_address > channelMap;

  // Read text file with BNL pin mapping
  std::ifstream BNLFile;
  // Variables to read the text file
  int wire, adapterPin, analogConnector, analogPin, asic, asicch, fembch, femb, wib;
  std::string plane, adapterConnector;
  BNLFile.open( bnl_map_filename );
  // Expects LArIAT_Pin_Mapping.xlsx to be saved as plain text file with columns separated by spaces 
  while( BNLFile >> wire >> plane >> adapterConnector >> adapterPin >> analogConnector >> analogPin >> asic >> asicch >> fembch >> femb >> wib ){
    // Build same key as the one produced by the BNL electronics running in channel-map mode
    // ASIC number goes from 1-8 in text file. We subtract 1 to make it 0-7 (used by the binary file)
    channel_address thisChannelAddress( wire, plane, adapterConnector, adapterPin, analogConnector, analogPin, asic -1, asicch, fembch, femb, wib ); 
    int BNLKey = ( (thisChannelAddress.fFEMB & 0xF) << 8 ) + ( (thisChannelAddress.fASIC & 0xF) << 4 ) + (thisChannelAddress.fASICch & 0xF);

    // Test if the key already existed when filling the map
    std::pair< std::map< int, channel_address >::iterator, bool > retVal;
    retVal = channelMap.insert( std::pair< int, channel_address >(BNLKey, thisChannelAddress) );
    if( retVal.second == false ){
      std::cerr << "ERROR: the triplet FEMB " << (thisChannelAddress.fFEMB & 0xF)
		<< " ASIC " << ((thisChannelAddress.fASIC & 0xF) + 1) // ASIC number goes from 1-8 in text file
		<< " ASICch " << (thisChannelAddress.fASICch & 0xF) << " is not unique in " << bnl_map_filename << std::endl;
      exit(1);
    }
  }
  BNLFile.close();

  /*
  // Sanity check
  for( std::map< int, channel_address >::iterator it = channelMap.begin(); it != channelMap.end(); ++it){
    channel_address& temp = it->second;
    if( it->first != (( (temp.fFEMB & 0xF) << 8 ) + ( (temp.fASIC & 0xF) << 4 ) + ( temp.fASICch & 0xF )) ){
      std::cerr << "ERROR: Channel map key " << it->first << " and values are not consistent" <<std::endl;
      exit(1);
    }
  }
  */
  // Done with the BNL pin mapping

  // ROOT file with the mapper run
  TFile inFile( mapper_run, "READ" );
  if( !inFile.IsOpen() ){
    std::cerr << "Unable to open file: " << mapper_run << std::endl;
    exit(1);
  }
  else std::cout << "Opening file: " << mapper_run << std::endl;

  // Get the input tree
  const char* inTreeName = "decoderTree";
  TTree *inTree = (TTree*)inFile.Get(inTreeName);
  if( !inTree ){
    std::cerr << "Tree not found: " << inTreeName << std::endl;
    exit(1);
  }
  else std::cout << "Tree found: " << inTreeName << std::endl;

  std::vector< std::vector<uint16_t> >* waveform = NULL;
  inTree->SetBranchAddress("waveform", &waveform);
  HeaderInfo* hinfo = NULL;
  inTree->SetBranchAddress("header", &hinfo);

  uint16_t useEvent = 1; // Event number to be used for building the map

  int bASICch = -1; // BNL's ASIC channel
  int bASICno = -1; // BNL's ASIC number
  int bFEMBno = -1; // BNL's FEMB number

  // Loop over FEMs in one event
  int entry = 0;
  while( inTree->GetEntry(entry) && hinfo->event == useEvent ){
    std::cout << "Processing entry " << entry << ": event " << hinfo->event << ", FEM " << (int)(hinfo->slot) << std::endl;

    if( (*waveform).size() != 64 ){
      std::cerr << "ERROR: less than 64 channels found in FEM" << std::endl; 
      exit(1);
    }
    
    // Loop over channels in one FEM
    for(size_t ich = 0; ich < (*waveform).size(); ich++){
      std::cout << "\tReading channel " << ich << std::endl;
      int lastADC = 0;
      // Loop over ADCs in one channel
      // Check all ADCs have the same value
      for(size_t t = 0; t < (*waveform)[ich].size(); t++){
	if( t > 0 && lastADC != (int)(*waveform)[ich][t] ){
	  std::cerr << "\tERROR: ADC do not remain constant in channel" << std::endl; 
	  exit(1);
	}
	lastADC = (int)(*waveform)[ich][t];
      }
      bASICch = lastADC & 0xF;
      bASICno = (lastADC >> 4) & 0xF;
      // The maximum number of ASICs per FEMB is 8
      if( bASICno > 7 ){
	std::cerr << "\tERROR: ADC data is wrong (ASIC number is out of bounds)" << std::endl; 
	exit(1);
      }
      bFEMBno = (lastADC >> 8) & 0xF;
      // The maximum number of FEMBs in LArIAT VST is 5
      if( bFEMBno > 5 ){
	std::cerr << "\tERROR: ADC data is wrong (FEMB number is out of bounds)" << std::endl; 
	exit(1);	
      }
      // std:: cout << "FEM " << std::setw(2) << (int)(hinfo->slot) << " Channel " << std::setw(2) << ich 
      // 		 << " ASIC channel " << std::setw(2) << bASICch << " ASIC number " << bASICno << " FEMB number " << bFEMBno << std::endl; 
      
      // Add the Nevis FEM slot and channel to the channel map
      // Do not check for existence since some FEM channels are not connected to a wire and do not appear in the text file
      (channelMap[lastADC]).fSlot = (int)(hinfo->slot);
      (channelMap[lastADC]).fFEMch = (int)ich;
      // Write the ASIC channel, number and FEMB number to make sure all channels have them
      (channelMap[lastADC]).fASICch = bASICch;
      (channelMap[lastADC]).fASIC = bASICno;
      (channelMap[lastADC]).fFEMB = bFEMBno;

    } // end of loop over channels in one FEM
    entry++;
  } // end of loop over FEMs in one event
  // Done with the ROOT file with the mapper run

  // Print the channel map
  std::ofstream channelMap_file;

  std::string name1(mapper_run);
  std::string name2(bnl_map_filename);
  // Base name. Add extension later
  std::string channelMap_filename = name1.substr(0, name1.find_last_of("/")) + "/"
    + "channelMap_" + name1.substr( name1.find_last_of("/") + 1, name1.find_last_of(".") - name1.find_last_of("/") - 1 ) + "_"
    + name2.substr( name2.find_last_of("/") + 1, name2.find_last_of(".") - name2.find_last_of("/") - 1 );

  std::string channelMap_filename_txt = channelMap_filename + ".txt";
  std::cout << "Writing channel map as text file to " << channelMap_filename_txt << std::endl;
  channelMap_file.open( channelMap_filename_txt, std::ios::out );
  for( std::map< int, channel_address >::iterator it = channelMap.begin(); it != channelMap.end(); ++it){
    channel_address& val = it->second;
    channelMap_file << std::right << std::setw(3) << val.fWire << " " << std::setw(10) << val.fPlane << " " 
		     << std::setw(3) << val.fAdapterConnector << " " << std::setw(2) << val.fAdapterPin << " " 
		     << std::setw(2) << val.fAnalogConnector << " " << std::setw(2) << val.fAnalogPin << " "
		     << std::setw(2) << val.fASIC << " " << std::setw(2) << val.fASICch << " "
		     << std::setw(3) << val.fFEMBch << " " << std::setw(2) << val.fFEMB << " " << std::setw(2) << val.fWIB << " "
		     << std::setw(2) << val.fSlot << " " << std::setw(2) << val.fFEMch << std::endl;
  }
  channelMap_file.close();

  std::string channelMap_filename_root = channelMap_filename + ".root";
  std::cout << "Writing channel map as ROOT file to " << channelMap_filename_root << std::endl;
  TFile outFile( channelMap_filename_root.c_str(), "RECREATE" );
  TTree* mapTree = new TTree( "channelMapTree", "Tree with the LArIAT VST channel map" );
  mapTree->ReadFile( channelMap_filename_txt.c_str(), 
		     "wire/I:plane/C:adapterconnector/C:adapterpin/I:analogconnector:analogpin:asic:asicch:fembch:femb:wib:slot:femch" );
  mapTree->Write();
  outFile.Close();

  return 0;
}

// To run as a standalone application
# ifndef __CINT__
int main( int argc, char** argv ){
  if( argc != 3 ){
    std::cerr << "Usage ./channel_mapper.exe DECODED_MAPPER_RUN.root BNL_PIN_MAP.txt" << std::endl;
    exit(1);
  }
  // To create interactive windows to see the plots
  TRint theApp( "tapp", &argc, argv );
  int status = channel_mapper( theApp.Argv(1), theApp.Argv(2) ); // TRint modifies argc & argv!
  theApp.Run();
  return status;
}
# endif
