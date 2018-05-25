#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <vector>

#include <TStyle.h>
#include <TRint.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>

#include "HeaderInfo.hh"

int plotter( const char* argv ){
  gStyle->SetOptTitle(1); // Title in canvas
  gStyle->SetLineWidth(1); // Thinnest lines
  gStyle->SetTitleOffset(1.2*(gStyle->GetTitleOffset("y")), "y"); // Add 20% offset to the vertical title
  gStyle->SetPadRightMargin(0.75*(gStyle->GetPadRightMargin()));
  gStyle->SetPadLeftMargin(1.17*(gStyle->GetPadLeftMargin()));

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

  std::vector< std::vector<uint16_t> >* waveform = NULL;
  inTree->SetBranchAddress("waveform", &waveform);
  HeaderInfo* hinfo = NULL;
  inTree->SetBranchAddress("header", &hinfo);

  // Canvas showing 64 channels
  TCanvas* canvas64 = new TCanvas("wfdisplay64", "Waveform display: 64 channels", 1024, 700);
  int colPad64 = 4; // Number of columns of pads
  int rowPad64 = 2; // Number of rows of pads
  canvas64->Divide(colPad64, rowPad64);

  // Four canvases showing 16 channels each
  std::vector<TCanvas*> vcanvas16(4);
  int colPad16 = 4; // Number of columns of pads
  int rowPad16 = 4; // Number of rows of pads
  for(int c = 0; c < (int)vcanvas16.size(); c++){
    vcanvas16.at(c) = new TCanvas( Form("wfdisplay16_%i", c), Form("Waveform display: %i-%i channels", c*16, c*16 + 15), 1024, 700);
    vcanvas16.at(c)->Divide(colPad16, rowPad16);
  }

  std::string prompt;
  int entry = 0;
  int maxEntry = inTree->GetEntries() - 1;
  int chPad64 = 64/(colPad64 * rowPad64); // Number of channels per pad
  // int chPad16 = 64/(int(vcanvas16.size()) * colPad16 * rowPad16); // Number of channels per pad
  std::vector<TMultiGraph*> vmgr(colPad64 * rowPad64); // Vector of pointers for memory management
  //  std::vector<TGraph*> vgr2(int(vcanvas16.size()) * colPad16 * rowPad16);
  std::vector<TGraph*> vgr2; // Vector of pointers for memory management
  // By using TGraph, it is already assumed that chPad16 will be 1
  // For a more generic approach, TMultiGraphs should be used

  while( entry >= 0 && entry <= maxEntry ){
    inTree->GetEntry(entry);
    std::cout << "Processing entry " << entry << ": event " << hinfo->event << ", FEM " << (int)(hinfo->slot) << std::endl;

    for(int ich = 0; ich < 64; ich++){
      TMultiGraph*& mgr = vmgr.at( (size_t)floor(ich/chPad64) ); // Reference pointer
      if( (ich % chPad64) == 0 ){
	mgr = new TMultiGraph( Form("ev%ifem%ich%i_%i", hinfo->event, (int)(hinfo->slot), ich, ich + chPad64 -1),
			       Form("Event %i FEM %i channels %i-%i; Time (#times 500 ns); ADC + 100 #times channel #", 
				    hinfo->event, (int)(hinfo->slot), ich, ich + chPad64 - 1) );
      }
      size_t ichsize = (*waveform)[ich].size();
      std::vector<int> time(ichsize);
      std::vector<int> wf(ichsize); // Because TGraph does not accept uint16_t
      std::vector<int> wfShifted(ichsize); // Shifted copy of the waveform so channels do not overlap
      if( ichsize > 0 ){
	// Fill the time vector
	for(size_t t = 0; t < ichsize; t++){
	  wf[t] = (int)(*waveform)[ich][t];
	  wfShifted[t] = wf[t] + ich*100; // Add offset (100.*ich) so channels do not overlap
	  time[t] = (int)t;
	}
	TGraph *graph = new TGraph( wfShifted.size(), &(time[0]), &(wfShifted[0]) );
	graph->SetMarkerColor(ich % 8 + 1); // Use ROOT colors 1 - 8
	graph->SetLineColor(ich % 8 + 1);
	graph->SetNameTitle(Form("ch%i", ich), Form("Channel %i", ich));
	mgr->Add(graph);

	// Subset of channels in canvas
	size_t canvasIndex = (size_t)floor(ich/(colPad16 * rowPad16));
	//vcanvas16.at( canvasIndex )->cd( colPad16*(ich%rowPad16) + (int)floor(ich/rowPad16) + 1 ); // Plot 64 channels in columns first
	vcanvas16.at( canvasIndex )->cd( ich%(colPad16 * rowPad16) + 1 );
	vgr2.push_back( new TGraph( wf.size(), &(time[0]), &(wf[0]) ) );
	TGraph* graph2 = vgr2.back();
	graph2->SetMarkerColor(ich % 8 + 1); // Use ROOT colors 1 - 8
	graph2->SetLineColor(ich % 8 + 1);
	graph2->SetNameTitle(Form("ch%i", ich), Form("Event %i FEM %i Channel %i; Time (#times 500 ns); ADC", hinfo->event, (int)(hinfo->slot), ich));
	graph2->Draw("APL");
	//vcanvas16.at( canvasIndex )->Update();
      }
      if( ((ich + 1) % chPad64) == 0 ){
	canvas64->cd( (int)floor(ich/chPad64) + 1 ); // Plot "chPad64" channels per pad
	mgr->Draw("APL");
	//canvas64->Update();
      }
    }
    canvas64->Update();
    //canvas64->Print(".png");
    for(size_t c = 0; c < vcanvas16.size(); c++){
      vcanvas16.at(c)->Update();
    }
    //vcanvas16.at(0)->SaveAs(Form("wfdisplay16_1_ev%i.png", hinfo->event ));

    std::cout << "\nEnter \"n\" for next entry\n";
    std::cout << "Enter \"exit\" to return to ROOT command line\n";
    std::cout << "Enter \"p\" for previous entryt\n";
    std::cout << "Enter entry number in [0, " << maxEntry << "] to go to that entry" << std::endl;

    std::cin >> prompt;
    if( prompt == "n" ){
      if( entry < maxEntry ) entry++;
      else std::cout << "You are displaying the last entry" << std::endl;
    }
    else if( prompt == "exit" ) break;
    else if( prompt == "p" ){
      if( entry > 1 ) entry--;
      else std::cout << "You are displaying the first entry" << std::endl;
    }
    else{
      try{
	int aux = std::stoi(prompt);
	if( aux >= 0 && aux <= maxEntry ) entry = aux;
	else std::cout << "Entry number out of bounds" << std::endl;
      }catch(std::invalid_argument& e){
	std::cout << "Input not recognized" << std::endl;
      }
    }

    // Clean memory
    for(size_t m = 0; m < vmgr.size(); m++){
      delete vmgr[m]; // Deleting the TMultiGraph deletes all the TGraphs associated
    }
    for(size_t g = 0; g < vgr2.size(); g++){
      delete vgr2[g];
    }
    vgr2.clear();

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
