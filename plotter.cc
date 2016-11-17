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

int plotter( const char* argv ){
  gStyle->SetOptTitle(1); // Title in canvas
  gStyle->SetLineWidth(1); // Thinnest lines

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
  int event = 1;
  int maxEvent = inTree->GetEntries();
  char cmd1[512], cmd2[512]; // Buffers for commands in TTree::Draw
  int chPad64 = 64/(colPad64 * rowPad64); // Number of channels per pad
  int chPad16 = 64/(int(vcanvas16.size()) * colPad16 * rowPad16); // Number of channels per pad
  std::vector<TMultiGraph*> vmgr(colPad64 * rowPad64); // Vector of pointers for memory management
  //  std::vector<TGraph*> vgr2(int(vcanvas16.size()) * colPad16 * rowPad16);
  std::vector<TGraph*> vgr2; // Vector of pointers for memory management
  // By using TGraph, it is already assumed that chPad16 will be 1
  // For a more generic approach, TMultiGraphs should be used

  while( event >= 1 && event <= maxEvent ){
    sprintf(cmd2, "Entry$==%i", event - 1); // FEM starts counting events from 1 but TTree starts counting entries from 0
    for(int ich = 0; ich < 64; ich++){
      // All channels in same canvas
      canvas64->cd( (int)floor(ich/chPad64) + 1 ); // Plot "chPad64" channels per pad
      TMultiGraph*& mgr = vmgr.at( (size_t)floor(ich/chPad64) ); // Reference pointer
      sprintf(cmd1, "waveform[%i][] + 100.*%i:Iteration$", ich, ich); // Add offset (100.*ich) so channels do not overlap
      //inTree->Draw("waveform[0][]:Iteration$","Entry$==0");
      if( (ich % chPad64) == 0 ){
	inTree->Draw(cmd1, cmd2, "goff");
	mgr = new TMultiGraph( Form("ev%ich%i_%i", event, ich, ich + chPad64 -1),
			       Form("Event %i channels %i-%i; Time (#times 500 ns); ADC + 100 #times channel #", event, ich, ich + chPad64 -1) );
      }
      else inTree->Draw(cmd1, cmd2, "samegoff");
      // Capture output of TTree::Draw
      if( inTree->GetSelectedRows() ){
	TGraph *graph = new TGraph( inTree->GetSelectedRows(), inTree->GetV2(), inTree->GetV1() );
	graph->SetMarkerColor(ich % 8 + 1); // Use ROOT colors 1 - 8
	graph->SetLineColor(ich % 8 + 1);
	graph->SetNameTitle(Form("ch%i", ich), Form("Channel %i", ich));
	mgr->Add(graph);
      }
      if( ((ich + 1) % chPad64) == 0 ) mgr->Draw("APL");
      //canvas64->Update();

      // Subset of channels in canvas
      size_t canvasIndex = (size_t)floor(ich/(colPad16 * rowPad16));
      //vcanvas16.at( canvasIndex )->cd( colPad16*(ich%rowPad16) + (int)floor(ich/rowPad16) + 1 ); // Plot 64 channels in columns first
      vcanvas16.at( canvasIndex )->cd( ich%(colPad16 * rowPad16) + 1 );
      sprintf(cmd1, "waveform[%i][]:Iteration$", ich);
      if( (ich % chPad16) == 0 ){
	inTree->Draw(cmd1, cmd2, "goff");
      }
      else inTree->Draw(cmd1, cmd2, "samegoff");
      // Capture output of TTree::Draw
      if( inTree->GetSelectedRows() ){
	//TGraph *graph2 = new TGraph( inTree->GetSelectedRows(), inTree->GetV2(), inTree->GetV1() );
	vgr2.push_back( new TGraph( inTree->GetSelectedRows(), inTree->GetV2(), inTree->GetV1() ) );
	TGraph* graph2 = vgr2.back();
	graph2->SetMarkerColor(ich % 8 + 1); // Use ROOT colors 1 - 8
	graph2->SetLineColor(ich % 8 + 1);
	graph2->SetNameTitle(Form("ch%i", ich), Form("Event %i Channel %i; Time (#times 500 ns); ADC", event, ich));
	graph2->Draw("APL");
      //vcanvas16.at( canvasIndex )->Update();
      }
    }
    canvas64->Update();
    for(size_t c = 0; c < vcanvas16.size(); c++){
      vcanvas16.at(c)->Update();
    }
    //vcanvas16.at(0)->SaveAs(Form("wfdisplay16_1_ev%i.png", event ));

    std::cout << "\nEnter \"n\" for next event\n";
    std::cout << "Enter \"exit\" to return to ROOT command line\n";
    std::cout << "Enter \"p\" for previous event\n";
    std::cout << "Enter event number in [1, " << maxEvent << "] to go to that event" << std::endl;

    std::cin >> prompt;
    if( prompt == "n" ){
      if( event < maxEvent ) event++;
      else std::cout << "You are displaying the last event" << std::endl;
    }
    else if( prompt == "exit" ) break;
    else if( prompt == "p" ){
      if( event > 1 ) event--;
      else std::cout << "You are displaying the first event" << std::endl;
    }
    else{
      try{
	int aux = std::stoi(prompt);
	if( aux >= 1 && aux <= maxEvent ) event = aux;
	else std::cout << "Event number out of bounds" << std::endl;
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
