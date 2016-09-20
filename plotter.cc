#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>

#include <TStyle.h>
#include <TRint.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>

int plotter( const char* argv ){
  gStyle->SetOptTitle(1); // Title in canvas

  TFile inFile( argv, "READ" );
  if( !inFile.IsOpen() ){
    std::cerr << "Unable to open file: " << argv << std::endl; 
    exit(0);
  }
  else std::cout << "Opening file: " << argv << std::endl;
 
  // Get the input tree
  const char* inTreeName = "decoderTree";
  TTree *inTree = (TTree*)inFile.Get(inTreeName);
  if( !inTree ){
    std::cerr << "Tree not found: " << inTreeName << std::endl; 
    exit(0); 
  }
  else std::cout << "Tree found: " << inTreeName << std::endl;

  TCanvas* canvas = new TCanvas("wfdisplay", "Waveform display", 1600, 700);
  canvas->Divide(4,2);
  std::string prompt;
  int event = 0;
  int maxEvent = inTree->GetEntries();
  char cmd1[512], cmd2[512]; // Buffers for commands in TTree::Draw
  TMultiGraph* mgr;
  
  while( prompt != "exit" && event < maxEvent ){
    sprintf(cmd2, "Entry$==%i", event);
    for(int ich = 0; ich < 64; ich++){
      canvas->cd((int)floor(ich/8) + 1); // Plot 8 channels per pad
      sprintf(cmd1, "waveform[%i][]:Iteration$", ich);
      //inTree->Draw("waveform[0][]:Iteration$","Entry$==0");
      if( (ich % 8) == 0 ){ 
	inTree->Draw(cmd1, cmd2);
	mgr = new TMultiGraph( Form("ev%ich%i_%i", event, ich, ich + 7), 
			       Form("Event %i channels %i-%i; Ticks; ADC", event, ich, ich + 7) );
      }
      else inTree->Draw(cmd1, cmd2, "same");
      // Capture output of TTree::Draw
      TGraph *graph = new TGraph( inTree->GetSelectedRows(), inTree->GetV2(), inTree->GetV1() );
      graph->SetMarkerColor(ich % 8 + 1);
      graph->SetLineColor(ich % 8 + 1);
      graph->SetNameTitle(Form("ch%i", ich), Form("Channel %i", ich));
      mgr->Add(graph);
      if( ((ich + 1) % 8) == 0 ) mgr->Draw("APL");
      canvas->Update();
    }
    event++;
    std::cout << "Enter \"exit\" to quit" << std::endl;
    std::cin >> prompt;
  }
  return 1;
}

// To run as a standalone application
# ifndef __CINT__
int main( int argc, char** argv ){
  // To create interactive windows to see the plots
  TRint theApp( "tapp", &argc, argv );
  int status = plotter( theApp.Argv(1) ); // TRint modifies argc & argv!
  theApp.Run();
  return status;
}
# endif

