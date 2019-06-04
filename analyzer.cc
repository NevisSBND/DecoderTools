#include <iostream>

#include <TROOT.h>
#include <TRint.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>

#include "HeaderInfo.hh"

int analyzer( const char* runFile ){

  double triggerRate = 0.2; // in Hz
  double frameLength = 2560.; // in 2 MHz samples
  int NFEMs = 10; // number of FEMs
  int firstFEM = 4; // slot of the first FEM

  // Input ROOT file
  TFile inFile( runFile, "READ" );
  if( !inFile.IsOpen() ){
    std::cerr << "Unable to open file: " << runFile << std::endl;
    exit(1);
  }
  else std::cout << "Opening file: " << runFile << std::endl;

  // Get the input tree
  const char* inTreeName = "decoderTree";
  TTree *inTree = (TTree*)inFile.Get(inTreeName);
  if( !inTree ){
    std::cerr << "Tree not found: " << inTreeName << std::endl;
    exit(1);
  }
  else std::cout << "Tree found: " << inTreeName << std::endl;

  inTree->SetBranchStatus("waveform", 0); // Speed up by not reading the waveform
  HeaderInfo* hinfo = NULL;
  inTree->SetBranchAddress("header", &hinfo);

  int entries = inTree->GetEntries();
  if( entries < 2 ){
    std::cerr << "Analyzer needs more than one entry to compute time interval" << std::endl;
    exit(1);
  }

  TH1D *hFrames = new TH1D("hFrames", "Event frames;Frame;Entries/1000 frames", 16778, 0, 16778000);
  hFrames->SetDirectory(gROOT);
  TH1D *hDeltaFrame = new TH1D("hDeltaFrame", "Frame difference;#DeltaFrame;Entries/frame", 9999, 1, 10000);
  hDeltaFrame->SetDirectory(gROOT);
  TH1D *hEvents = new TH1D("hEvents", "Event numbers;Event;Entries/1000 events", 16778, 0, 16778000);
  hEvents->SetDirectory(gROOT);
  TH1D *hDeltaEvent = new TH1D("hDeltaEvent", "Event no. difference;#DeltaEvent;Entries/events", 100, 0, 100);
  hDeltaEvent->SetDirectory(gROOT);

  inTree->GetEntry(0);
  hFrames->Fill((int)hinfo->frame);
  int prevFrame = (int)hinfo->frame;
  hEvents->Fill((int)hinfo->event);
  int prevEvent = (int)hinfo->event;

  int prevSlot = firstFEM;
  int thisSlot = firstFEM;
  int aux = 0;

  double triggerGaps = 0; // Unexpected gap between triggers (missing triggers?)
  double eventJumps = 0; // FEM event number header jumps by >=2
  double  totalTriggers = 0; // Total triggers

  // Loop over 1 to N entries
  for( int i = 1; i < entries; i++ ){
    inTree->GetEntry(i);
    int thisFrame = (int)hinfo->frame;
    hFrames->Fill(thisFrame);
    int deltaFrame = thisFrame - prevFrame;
    hDeltaFrame->Fill(deltaFrame);
    if( prevFrame != thisFrame ) totalTriggers++;
    prevFrame = thisFrame;

    int thisEvent = (int)hinfo->event;
    hEvents->Fill(thisEvent);
    int deltaEvent = thisEvent - prevEvent;
    hDeltaEvent->Fill(deltaEvent);
    prevEvent = thisEvent;

    thisSlot = (int)hinfo->slot;

    // Look for missed triggers when the frame difference between FEMs is bigger than the tolerance
    if( deltaFrame > 1.05/(triggerRate * frameLength * 0.5e-6) ){ // 1.05 --> add 5% tolerance
      std::cout << "\n\nTRIGGER GAP" << std::endl;
      std::cout << "\tPREVIOUS EVENT" << std::endl;
      for( int j = i - NFEMs; j < i; j++ ) inTree->Show(j);
      std::cout << "\tTHIS EVENT" << std::endl;
      for( int j = i; j < i + NFEMs; j++ ) inTree->Show(j);
      triggerGaps++;
      // If the frame difference occurs in the middle of a crate-readout, FEMs are desynchronized
      if( thisSlot != firstFEM && prevSlot != firstFEM + NFEMs - 1 ){
	std::cout << "DESYNC. Enter anything to continue" << std::endl;
	std::cin >> aux;
      }
    }

    // Look for event jumps
    if(deltaEvent > 1){
      std::cout << "\n\nMISSING EVENT" << std::endl;
      std::cout << "\tPREVIOUS EVENT" << std::endl;
      for( int j = i - NFEMs; j < i; j++ ) inTree->Show(j);
      std::cout << "\tTHIS EVENT" << std::endl;
      for( int j = i; j < i + NFEMs; j++ ) inTree->Show(j);
      eventJumps++;
      // If the event difference occurs in the middle of a crate-readout, FEMs are desynchronized
      if( thisSlot != firstFEM && prevSlot != firstFEM + NFEMs - 1 ){
	std::cout << "DESYNC. Enter anything to continue" << std::endl;
	std::cin >> aux;
      }
    }

    inTree->GetEntry(i);
    prevSlot = (int)hinfo->slot;
  } // end of loop over N entries

  std::cout << "--- Input parameters ---" << std::endl;
  std::cout << "Trigger rate: " << triggerRate << " Hz" << std::endl;
  std::cout << "Frame length: " << frameLength << " samples" << std::endl;
  std::cout << "Number of FEMs: " << NFEMs << std::endl;

  std::cout << "--- Output parameters ---" << std::endl;
  std::cout << "Missed triggers: " << triggerGaps << std::endl;
  std::cout << "Event jumps: " << eventJumps << std::endl;
  std::cout << "Total triggers " << totalTriggers << std::endl;
  std::cout << "Fraction of missed triggers: " << triggerGaps/totalTriggers << std::endl;
  std::cout << "Fraction of event jumps: " << eventJumps/totalTriggers << std::endl;

  inTree->ResetBranchAddresses(); // "detach" from local variables

  std::string outFileName(runFile);
  outFileName = outFileName.substr(0, outFileName.find_last_of(".")) + "_ana.root";
  TFile rootFile( outFileName.c_str(), "RECREATE" );
  
  TCanvas *cFrame = new TCanvas("cFrame", "cFrame");
  cFrame->Divide(1, 2);
  cFrame->cd(1);
  hFrames->Draw();
  cFrame->cd(2);
  hDeltaFrame->Draw();
  cFrame->cd(0);
  cFrame->Modified(); cFrame->Update();
  cFrame->Write();

  TCanvas *cEvent = new TCanvas("cEvent", "cEvent");
  cEvent->Divide(1, 2);
  cEvent->cd(1);
  hEvents->Draw();
  cEvent->cd(2);
  hDeltaEvent->Draw();
  cEvent->cd(0);
  cEvent->Modified(); cEvent->Update();
  cEvent->Write();

  rootFile.Close();

  return 0;
}

// To run as a standalone application
# ifndef __CINT__
int main( int argc, char** argv ){
  if( argc != 2 ){
    std::cerr << "Usage ./analyzer.exe DECODED_RUN.root" << std::endl;
    exit(1);
  }
  // To create interactive windows to see the plots
  TRint theApp( "tapp", &argc, argv );
  int status = analyzer( theApp.Argv(1) ); // TRint modifies argc & argv!
  theApp.Run();
  return status;
}
# endif
