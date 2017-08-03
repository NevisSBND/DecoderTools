#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>

#include <TROOT.h>
#include <TStyle.h>
#include <TRint.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TProfile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TPaveStats.h>


#include "XMITInfo.hh"

int analyzer( int argc, char** argv ){

  // To create interactive windows to see the plots
  TRint theApp( "tapp", &argc, argv );

  /////////////
  // Initialize
  /////////////

  //gStyle->SetOptTitle(1); // Title in canvas
  gStyle->SetOptStat(1); // Stats in canvas
  gStyle->SetLineWidth(1); // Thinnest lines
  gStyle->SetHistLineColor(kRed);
  gStyle->SetTitleOffset(1.2*(gStyle->GetTitleOffset("y")), "y"); // Add 20% offset to the vertical title
  gStyle->SetPadRightMargin(0.75*(gStyle->GetPadRightMargin()));
  gStyle->SetPadLeftMargin(1.17*(gStyle->GetPadLeftMargin()));
  gROOT->SetBatch(false);
  bool verbose = false;
  //bool verbose = true;

  if(theApp.Argc() != 3 ){
    std::cerr << "Run it as ./analyzer.exe $SEB_NUMBER $RUN_NUMBER " << std::endl;
    exit(1);
  }

  std::string sebname(theApp.Argv(1));
  std::string runname(theApp.Argv(2));

  std::string inFileName = "seb-" + sebname + "-" + runname + "-hadd_dissect_v2.root";

  TFile inFile( inFileName.c_str(), "READ" );
  if( !inFile.IsOpen() ){
    std::cerr << "Unable to open file: " << argv << std::endl;
    exit(1);
  }
  else std::cout << "Opening file: " << argv << std::endl;

  // Get the input tree
  const char* inTreeName = "dissecterTree";
  TTree *inTree = (TTree*)inFile.Get(inTreeName);
  if( !inTree ){
    std::cerr << "Tree not found: " << inTreeName << std::endl;
    exit(1);
  }
  else std::cout << "Tree found: " << inTreeName << std::endl;

  std::string outFileName = inFileName.substr(0, inFileName.find_last_of(".")) + "_analyzed.root";
  TFile outFile( outFileName.c_str(), "RECREATE" );

  TCanvas c1("c1", "c1");
  inTree->Draw("slot:(frame*1.6e-3)>>hOverflowSlotFrac(3600,0,3600,15,4,19)","overflow*1.6e-3","colz"); // Fraction of overflows per second vs slot

  TCanvas c2("c2", "c2");
  inTree->Draw("slot:(frame*1.6e-3)>>hFullSlotFrac(3600,0,3600,15,4,19)","full*1.6e-3","colz"); // Fraction of fulls per second vs slot

  TCanvas c3("c3", "c3");
  inTree->Draw("slot:(frame*1.6e-3)>>hDataSlotRate(3600,0,3600,15,4,19)","2.*(wordcount+12.)*1e-6","colz"); // Rate of data vs slot

  TCanvas c4("c4", "c4");
  inTree->Draw("(frame*1.6e-3)>>hOverflowFrac(3600,0,3600)","overflow*1.6e-3",""); // Fraction of overflows per second vs slot

  TCanvas c5("c5", "c5");
  inTree->Draw("(frame*1.6e-3)>>hFullFrac(3600,0,3600)","full*1.6e-3",""); // Fraction of full per second vs slot

  TCanvas c6("c6", "c6");
  inTree->Draw("(frame*1.6e-3)>>hDataRate(3600,0,3600)","2.*(wordcount+12.)*1e-6",""); // Rate of data

  TCanvas c7("c7", "c7");
  inTree->Draw("slot:(frame*1.6e-3)>>hWordAsymSlot(3600,0,3600,15,4,19)","0.5*(wordcount-nwords)/(wordcount+nwords)","colz"); // Word asymmetry per second vs slot

  TCanvas c8("c8", "c8");
  inTree->Draw("(frame*1.6e-3)>>hWordAsym(3600,0,3600)","0.5*(wordcount-nwords)/(wordcount+nwords)",""); // Word asymmetry per second

  TCanvas c9("c9", "c9");
  inTree->Draw("slot:(frame*1.6e-3)>>hEventAsymSlot(3600,0,3600,15,4,19)","0.5*(frame-event)/(frame+event)","colz"); // Frame-event asymmetry per second vs slot

  TCanvas c10("c10", "c10");
  inTree->Draw("(frame*1.6e-3)>>hEventAsym(3600,0,3600)","0.5*(frame-event)/(frame+event)",""); // Frame-event asymmetry per second

  TCanvas c11("c11", "c11");
  inTree->Draw("slot:(frame*1.6e-3)>>hBadFrameSlot(3600,0,3600,15,4,19)","badframecount","colz"); // Frames per second with mismatched frame per slot

  TCanvas c12("c12", "c12");
  inTree->Draw("(frame*1.6e-3)>>hBadFrame(3600,0,3600)","badframecount",""); // Frames per second with mismatched frame

  // To do: give them titles, axis titles...

  outFile.Write();

  XMITInfo* xmitinfo = NULL;
  inTree->SetBranchAddress("xmitinfo", &xmitinfo);
  XMITInfo* last_xmitinfo = new XMITInfo();

  // Retrieve first entry to extract information for booking ROOT objects
  inTree->GetEntry(0);
  // Vectors of vectors for TGraphs: first index denotes FEM board, second index denotes frame
  std::vector< std::vector< double > > vvFEMFrameTime;
  vvFEMFrameTime.resize( xmitinfo->size() );
  std::vector< std::vector< double > > vvFEMDataRate;
  vvFEMDataRate.resize( xmitinfo->size() );
  std::vector<bool> first_overflow, first_full, first_overflow_full;
  for( size_t ifem = 0; ifem < xmitinfo->size(); ifem++ ){
    vvFEMFrameTime.at(ifem).reserve(inTree->GetEntries());
    vvFEMDataRate.at(ifem).reserve(inTree->GetEntries());
    first_overflow.push_back(false);
    first_full.push_back(false);
    first_overflow_full.push_back(false);
  }
  const size_t firstFEMSlot = (size_t)((*xmitinfo)[0].slot);
  // Vectors for total (crate, all FEM boards)
  std::vector< double > vTotalFrameTime;
  vTotalFrameTime.reserve(inTree->GetEntries());
  std::vector< double > vTotalDataRate;
  vTotalDataRate.reserve(inTree->GetEntries());

  TProfile* pCrateDataRate = new TProfile("pCrateDataRate", 
					  Form("Crate %s data date; Time (s); Crate %s data rate (MB/s)", sebname.c_str(), sebname.c_str()), 3600, 0., 3600.);

  TH1D* hCrateCompressionFactor = new TH1D("hCrateCompressionFactor",
					   Form("Crate %s compression factor; Crate %s compression factor; Entries", sebname.c_str(), sebname.c_str()), 1000, 0, 1.);

  /////////////
  // Event loop
  /////////////

  // Loop over tree entries (XMITInfo objects)
  for( int ientry = 0; ientry < inTree->GetEntries(); ientry++ ){ 
  //  for( int ientry = 0; ientry < 5000; ientry++ ){ //debug
    inTree->GetEntry(ientry);
    if( verbose ) std::cout << "XMIT packet " << ientry << " with " << xmitinfo->size() << " FEMs" << std::endl;
    if( ientry > 0 && (xmitinfo->size() != last_xmitinfo->size()) ){
      std::cerr << "WARNING: XMIT packet " << ientry << " has a different number of FEMs (" 
		<< xmitinfo->size() << ") from previous packet (" << last_xmitinfo->size() << ")" << std::endl; 
      // Note: Found some FEMs within the same frame splitted across two XMIT packets. Might be due to two Huffman words FFFF FFFF being interpreted as XMIT header
    }
    double crateDataRate = 0;
    double crateWordCount = 0;
    double crateNoCompressionWordCount = 0;

    uint32_t last_frame = 0;
    // Loop over FEMs
    for( size_t ifem = 0; ifem < xmitinfo->size(); ifem++ ){
      FEMInfo feminfo = (*xmitinfo)[ifem];
      uint8_t slot = feminfo.slot; // FEM slot in crate
      if( verbose ) std::cout << "FEM at slot " << (int)slot << std::endl;
      bool overflow = feminfo.overflow; // Overflow flag
      bool full = feminfo.full; // Full flag

      uint32_t nwords = feminfo.nwords; // Number of words in frame
      uint32_t event = feminfo.event; // Event number
      uint32_t frame = feminfo.frame; // Frame number
      if( ifem > 0 && (frame != last_frame) ){
	std:: cerr << "WARNING: XMIT packet " << ientry << " contains FEMs with different frames! Frame at FEM at slot " << slot << " is " << frame 
		   << " while FEM at slot " << (*xmitinfo)[ifem - 1].slot << " is " << (*xmitinfo)[ifem - 1].frame << std::endl;
      }
      last_frame = frame;
      uint32_t wordcount = feminfo.wordcount; // Actual (counted) number of words in frame
      uint32_t channelcount = feminfo.channelcount; // Actual (counted) number of channels in frame

      if( !(first_overflow[slot - firstFEMSlot]) && overflow ){
	std::cout << "FEM at slot " << (int)feminfo.slot << " first overflow at frame " << frame << std::endl;
	first_overflow[slot - firstFEMSlot] = true;
      }
      if( !(first_full[slot - firstFEMSlot]) && full ){
	std::cout << "FEM at slot " << (int)feminfo.slot << " first full at frame " << frame << std::endl;
	first_full[slot - firstFEMSlot] = true;
      }
      if( !(first_overflow_full[slot - firstFEMSlot]) && overflow && full ){
	std::cout << "FEM at slot " << (int)feminfo.slot << " first overflow && full at frame " << frame << std::endl;
	first_overflow_full[slot - firstFEMSlot] = true;
      }

      if( ientry > 0 ){
	FEMInfo last_feminfo = (*last_xmitinfo)[ifem];
	if( slot != last_feminfo.slot ){
	  std::cerr << "WARNING: FEM order in XMIT packet " << ientry << " has changed! Expected FEM at slot " << (int)last_feminfo.slot 
		    << " but found FEM at slot " << (int)slot << std::endl;
	  vvFEMDataRate.at(slot - firstFEMSlot).push_back(-1);
	  vvFEMFrameTime.at(slot - firstFEMSlot).push_back(frame * 1.6e-3);
	  continue;
	} 
	int frameDiff = (int)(frame - last_feminfo.frame);
	if( verbose ) std::cout << "Frame diff " << frameDiff << std::endl;
	double femDataRate = ((int)wordcount + 12)*(2.e-6)/(frameDiff * 1.6e-3); // MB/s; + 12 for FEM headers; 1 word = 2 B; 1 frame = 1.6 ms 
	vvFEMDataRate.at(slot - firstFEMSlot).push_back(femDataRate);
	vvFEMFrameTime.at(slot - firstFEMSlot).push_back(frame * 1.6e-3);
	crateDataRate += femDataRate;
	double femCompressionFactor = wordcount/(frameDiff * channelcount * 3200.);
	crateWordCount += (double)wordcount;
	crateNoCompressionWordCount += double(frameDiff * channelcount * 3200);
      }
    } // End of loop over FEMs
  
    vTotalDataRate.push_back(crateDataRate);
    vTotalFrameTime.push_back((int)last_frame * 1.6e-3); // Warning: using frame number from last FEM for global crate frame

    pCrateDataRate->Fill( ((int)last_frame * 1.6e-3), crateDataRate ); // Warning: using frame number from last FEM for global crate frame

    if(crateNoCompressionWordCount > 0) hCrateCompressionFactor->Fill(crateWordCount/crateNoCompressionWordCount);

    *last_xmitinfo = *xmitinfo; // Copy last XMITInfo to compare
  } // End of loop over tree entries
  
  ////////////
  // Finalize
  ////////////

  TCanvas* c_hCrateCompressionFactor = new TCanvas(Form("seb-%s-%s-c_hCrateCompressionFactor", sebname.c_str(), runname.c_str()), Form("seb-%s-%s-c_hCrateCompressionFactor", runname.c_str(), sebname.c_str()));
  gPad->SetLogy();
  gPad->SetLogx();
  hCrateCompressionFactor->Draw();
  gPad->Update();
  TLatex tx_invComp(0, 0, Form("1/Mean = %.2f", 1./(hCrateCompressionFactor->GetMean())));
  TPaveStats *ps = (TPaveStats*)gPad->GetPrimitive("stats");
  ps->SetName(hCrateCompressionFactor->GetName());
  ps->GetListOfLines()->Add(&tx_invComp);
  hCrateCompressionFactor->SetStats(0);
  c_hCrateCompressionFactor->Modified();
  c_hCrateCompressionFactor->Update();
  //  gPad->Print(0,".png");

  hCrateCompressionFactor->Write();

  TCanvas* cDataRate = new TCanvas(Form("seb-%s-%s-cDataRate", sebname.c_str(), runname.c_str()), "Instantaneous data rate");

  TGraph* gTotalDataRate = new TGraph( vTotalDataRate.size(), &(vTotalFrameTime[0]), &(vTotalDataRate[0]) );
  gTotalDataRate->Draw("ALP");
  for( size_t ifem = 0; ifem < vvFEMDataRate.size(); ifem++ ){
    TGraph* gFEMDataRate = new TGraph( vvFEMDataRate.at(ifem).size(), &(vvFEMFrameTime.at(ifem)[0]), &(vvFEMDataRate.at(ifem)[0]) );
    gFEMDataRate->SetLineColor((int)ifem);
    gFEMDataRate->Draw("L");
  }
  cDataRate->Write();

  TCanvas* c_pCrateDataRate = new TCanvas(Form("seb-%s-%s-c_pCrateDataRate", sebname.c_str(), runname.c_str()), "Crate data rate (prof)");
  pCrateDataRate->Draw();
  pCrateDataRate->Write();
  c_pCrateDataRate->Write();

  outFile.Close();

  theApp.Run();

  return 0;
}

// To run as a standalone application
# ifndef __CINT__
int main( int argc, char** argv ){

  int status = analyzer( argc, argv );

  return status;
}
# endif
