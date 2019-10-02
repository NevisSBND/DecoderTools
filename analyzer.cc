//########################################################################
// Tailored for condor jobs, only generates output root file
// Plot comparison will be done by another script.
//            07/29/2019,Guanqun
//########################################################################
#include <iostream>
#include <algorithm>    //for operations on vector element
#include <numeric>
#include <math.h>
#include <vector>
#include <TROOT.h>
#include <TRint.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TVectorD.h>
#include <TH2C.h>
#include <TH2S.h>
#include <TMath.h>
#include <TCanvas.h>


#include "HeaderInfo.hh"

using namespace std;


//for test run, `encode_Huffman` function is not needed.
/*
//returns the number of words after Huffman compression for a waveform
long int encode_Huffman(std::vector<unsigned int> w){   
	long int nword=1; //number of words in the waveform after Huffman compression
	int H_size = 14; //number of available bits to save data in Huffman word, NU stream specifically
	int word_size=0;
	int size_added=0;  //length just added to the candidate Huffman word 
	for(int i=1; i< (int)w.size() ; i++){
		switch(w[i]-w[i-1]){
			case 0: { size_added = 1; word_size+=size_added; break; }
			case -1: { size_added = 2; word_size+=size_added; break; }
			case 1: { size_added = 3; word_size+=size_added; break; }
			case -2: { size_added = 4; word_size+=size_added; break; }
			case +2: { size_added = 5; word_size+=size_added; break; }
			case -3: { size_added = 6; word_size+=size_added; break; }
			case +3: { size_added = 7; word_size+=size_added; break; }
			default: {
				 if(word_size != 0) nword+=2;
				 else nword++;
				 word_size = 0;
			}
		}
		if(word_size > H_size){
			nword++;
			word_size = size_added;
		}
	}
	if(word_size != 0) nword++; //for the last word in the waveform
	
	return nword;   
}

*/

int analyzer( const char* runFile ){

  double triggerRate = 1; // in Hz
  double frameLength = 2560.; // in 2 MHz samples   //guanqun: number of ticks in 1 frame.
  double framesize = frameLength*0.5e-6;  //the length of 1 frame (in second). The length of 1 tick is 0.5us.
  int NFEMs = 16; // number of FEMs
  int firstFEM = 3; // slot of the first FEM
  int EventBin = 10;  //binning size of the time axis("Event" axis)

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

  //inTree->SetBranchStatus("waveform", 0); // Speed up by not reading the waveform
  std::vector<std::vector<uint16_t>>* hwaveform=0;
  HeaderInfo* hinfo = NULL;
  inTree->SetBranchAddress("header", &hinfo);
  inTree->SetBranchAddress("waveform", &hwaveform);

  int entries = inTree->GetEntries();
  if( entries < 2 ){
    std::cerr << "Analyzer needs more than one entry to compute time interval" << std::endl;
    exit(1);
  }
  //grab the smallest event number in decoded file for plotting purpose
  inTree->GetEntry(0);
  double Event_min = (double)hinfo->event;
  double totEvent = ceil(entries/NFEMs); //calculate the total number events in the run based on entries
  int Event_Start = floor(Event_min/EventBin)*EventBin;  //the start of the event axis.
  int EventRange = ceil(totEvent/EventBin)+1;  //calcualte num of bins needed for the "event" axis assuming using bin size 'EventBin'.
  std::cout << "There are total " << entries << " entries in the tree"<< std::endl;
  

  std::string outFileName(runFile);
  //std::string outFileName = "tryanalyzer_ana.root";
  std::string outPDFName= outFileName.substr(0, outFileName.find_last_of(".")) + "_ana.pdf";
  outFileName = outFileName.substr(0, outFileName.find_last_of(".")) + "_ana.root";
  TFile rootFile( outFileName.c_str(), "RECREATE" ); 
  std::vector<double> max_event;  //maximum event number in the input data file.
  std::vector<double> max_frame;  //maximum frame number in the input data file.

 
 // gROOT->cd();  //make sure all newly created histos are associated with gROOT directory.

  //comparison of frame number between different FEMs within one event
  TH1D* hFrameComp[NFEMs-1];
  for(int i=0; i< NFEMs-1; i++){
	hFrameComp[i] = new TH1D(Form("hFrameComp_%d_%d", firstFEM, firstFEM+i+1), Form("Frame difference between slot %d and %d within 1 event", firstFEM, firstFEM+i+1),10,-5,5);

  }

  //internal comparison of trigger sample within 1 event, expect 0 trigger sample difference.
   TH2D* hTriggerSamComp = new TH2D("hTriggerSamComp","Trigger sample difference between first slot and other slots within 1 event;Slot;trigger sample difference; entries", NFEMs-1, firstFEM+1, firstFEM+NFEMs,100,-500,500);

  //real-time trigger rate calculated from each slot.
  TH2D *hFrameTrigger=new TH2D("hFrameTrigger", "Real-time triggerRate calculated by the frame difference; Slot; Trigger Rate; entries", NFEMs, firstFEM, firstFEM+NFEMs, (int)8*triggerRate, 0, 4*triggerRate);
  

  //simply the distribution of all the variables for difference FEMs.
  TH2D *hSlot = new TH2D("hSlot", "Slot number; Slot; Manually counted Slot number; Entries", NFEMs, firstFEM, firstFEM+NFEMs, NFEMs, firstFEM, firstFEM+NFEMs);
  TH2D *hOverflow = new TH2D("hOverflow", "Overflow distribution;Slot; Overflow flag; Entries", NFEMs, firstFEM, firstFEM+NFEMs, 2, 0,2);
  TH2D *hFull = new TH2D("hFull", "FullDistribution; Slot;Full flag; Entries",NFEMs, firstFEM, firstFEM+NFEMs, 2,0,2);
  TH2D *hId = new TH2D("hId", "ID; Slot;ID; Entries", NFEMs, firstFEM, firstFEM+NFEMs,127,0,127);
  TH2D *hnwords = new TH2D("hnwords", "distribution of number of words; Slot; number of words; entries",NFEMs, firstFEM, firstFEM+NFEMs,10000, 320000, 330000);
  TH2D *hwordcount = new TH2D("hwordcount", "distribution of manually counted number of words;Slot; wordcount; entries",NFEMs, firstFEM, firstFEM+NFEMs, 10000, 0,1000000);
  TH2D *hdiff_nword = new TH2D("hdiff_nword", "number of words - wordcount;Slot; number of words - wordcount; entries", NFEMs, firstFEM, firstFEM+NFEMs,20,-10,10);
  TH2D *hFrames = new TH2D("hFrames", "Event Frame numbers;Frame/1000;Slot;Entries",16778, 0, 16778000, NFEMs, firstFEM, firstFEM+NFEMs);   //largest frame number is 16777215
 //hFrames->SetDirectory(gROOT);
  TH2D *hDeltaFrame = new TH2D("hDeltaFrame", "Frame difference;#DeltaFrame;Slot;Entries", 2000, 0, 2000, NFEMs, firstFEM, firstFEM+NFEMs);  
  //TH2D *hDeltaFrame = new TH2D("hDeltaFrame", "Frame difference;Slot;#DeltaFrame;Entries",NFEMs, firstFEM, firstFEM+NFEMs, 9999, 1, 10000);
  TH2D *hEvents = new TH2D("hEvents", "Event numbers;Event number/10 events;Slot;Entries",EventRange, Event_Start, Event_Start+EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);  //largest event number is 16777215
  TH2D *hDeltaEvent = new TH2D("hDeltaEvent", "Event no. difference;Slot;#DeltaEvent;Entries",NFEMs, firstFEM, firstFEM+NFEMs, 5, 0, 5);
  TH2D *hChecksum = new TH2D("hChecksum", "checksum from readout FEM header words;Slot; checksum; entries",NFEMs, firstFEM, firstFEM+NFEMs, 16778, 0, 16778000);
  TH2D *hdiff_checksum = new TH2D("hdiff_checksum", "checksum - mychecksum;Slot; checksum - mychecksum; entries",NFEMs, firstFEM, firstFEM+NFEMs, 8200,-4100,4100);
  TH2D *hMychecksum = new TH2D("hMychecksum", "checksum counted manually; Slot; Mychecksum; entries", NFEMs, firstFEM, firstFEM+NFEMs, 16778, 0, 16778000);
  TH2D *hTriggerFrame = new TH2D("hTriggerFrame", "distribution of Trigger frame;Slot; Trigger frame; entries",NFEMs, firstFEM, firstFEM+NFEMs, 17,0,17); //largest trigger frame is 15
  TH2D *hTriggerSample = new TH2D("hTriggerSample", "distribution of Trigger Sample;Slot; Trigger Sample; entries",NFEMs, firstFEM, firstFEM+NFEMs, 4100 ,0,4100); //largest trigger sample is 4095
  TH2D *hNumSample = new TH2D("hNumSample", "Number of samples in waveform  per channel; channel number; number of samples; entries", NFEMs*64,0,NFEMs*64, 1000,0,10000);
  TH2D *hMean = new TH2D("hMean", "Mean value of waveform per channel; channel number; mean value of waveform; entries", NFEMs*64,0,NFEMs*64, 1000,0,5000);
  TH2D *hRMS = new TH2D("hRMS", "RMS value of waveform per channel; channel number; RMS value; entries", NFEMs*64,0,NFEMs*64, 1000,0,5000);
  TH2D *hMax = new TH2D("hMax", "Max value of waveform per channel; channel number; Max value; entries",NFEMs*64,0,NFEMs*64, 1000,0,5000);
  TH2D *hMin = new TH2D("hMin", "Min value of waveform per channel; channel number; Min value; entries", NFEMs*64,0,NFEMs*64, 201,-2,400);
 
  // how those variables changes with time for each FEM or channel.
   TH2C *hTSlot = new TH2C("hTSlot", "2D plot of Slot number changing with time(Event number); Event number/10 events;Manually counted Slot number; Slot", EventRange,Event_Start, Event_Start+EventRange*EventBin,  NFEMs, firstFEM, firstFEM+NFEMs);
   TH2C *hTOverflow = new TH2C("hTOverflow", "2D plot of overflow vs time(Event number);Event number/10 events; Slot;Overflow",EventRange, Event_Start, Event_Start+EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);    //maximum event number: 16778000
   TH2C *hTFull = new TH2C("hTFull", "2D plot of Full flag vs. time(Event number); Event number/10 events;Slot; Full flag", EventRange,Event_Start, Event_Start+EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);
   TH2C *hTId = new TH2C("hTId", "2D plot of FEM id vs. time (Event number);Event number/10 events;Slot; FEM ID",EventRange,Event_Start, Event_Start+EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);
   TH2D *hTnwords = new TH2D("hTnwords", "2D plot of nwords vs. time(Event number); Event number/10 events;Slot; nwords",EventRange,Event_Start, Event_Start+EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);
   TH2D *hTwordcount = new TH2D("hTwordcount", "2D plot of wordcount vs. time(Event number);Event number/10 events;Slot; wordcount",EventRange,Event_Start, Event_Start+EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);
   TH2D *hTdiff_nword = new TH2D("hTdiff_nword", "2D plot of (nwords - wordcount) vs. time (Event number);Event number/10 events;Slot; nwords - wordcount", EventRange,Event_Start, Event_Start+EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);
   hTdiff_nword->SetAxisRange(-10*EventBin, 10*EventBin, "Z");
   TH2D *hTFrames = new TH2D("hTFrames", "2D plot of event frame number vs. time (Event number); Event number/10 events;Slot; Frame number",EventRange,Event_Start, Event_Start+EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);
   TH2D *hTDeltaFrame = new TH2D("hTDeltaFrame", "2D plot of Frame difference;Event number/10 events;Slot; Frame difference", EventRange,Event_Start, Event_Start+EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);
   hTDeltaFrame->SetMinimum(1);
   TH2D *hTEvents = new TH2D("hTEvents", "2D plot of Event number vs. time(Event number);Event number/10 events;Slot; Event number",EventRange,Event_Start, Event_Start+EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);
   TH2D *hTDeltaEvent = new TH2D("hTDeltaEvent", "2D plot of Event no. difference vs. time (Event number);Event number/10 events;Slot; Event difference",EventRange,Event_Start, Event_Start+EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs); 
   TH2D *hTChecksum = new TH2D("hTChecksum", "2D plot of checksum vs. time (Event number);Event number/10 events;Slot;Checksum",EventRange,Event_Start, Event_Start+EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);
   TH2D *hTMychecksum = new TH2D("hTMychecksum", "2D plot of manually counted checksum vs. time(Event number);Event number/10 events;Slot; mychecksum",EventRange,Event_Start, Event_Start+EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);
   TH2D *hTdiff_checksum = new TH2D("hTdiff_checksum", "2D plot of (checksum - mychecksum) vs. time(Event number);Event number/10 events;Slot;checksum - mychecksum",EventRange,Event_Start, Event_Start+EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);
   hTdiff_checksum->SetAxisRange(-4100*EventBin, 4100*EventBin, "Z"); 
   TH2D *hTTriggerRate = new TH2D("hTTriggerRate","2D plot of trigger rate vs. time;Time(s);Slot; Trigger Rate", 10750, 0, 21500, NFEMs, firstFEM, firstFEM+NFEMs);  //based on the largest frame number that can be saved in the Header, and the frame size, the longest time is ~21475 second.
   TH2C *hTTriggerFrame = new TH2C("hTTriggerFrame", "2D plot of Trigger Frame vs. time(Event number);Event number/10 events; Slot;TriggerFrame",EventRange,Event_Start, Event_Start+EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);
   TH2S *hTTriggerSample = new TH2S("hTTriggerSample", "2D plot of Trigger Sample vs. time(Event number);Event number/10 events;Slot;Trigger Sample",EventRange,Event_Start, Event_Start+EventRange*EventBin,  NFEMs, firstFEM, firstFEM+NFEMs);
   TH2D *hTNumSample = new TH2D("hTNumSample", "Number of samples per channel vs. time(Event number);Event number/10 events;channel number; Number of Samples",EventRange,Event_Start, Event_Start+EventRange*EventBin,  NFEMs*64, 0, NFEMs*64);
   TH2D *hTMean = new TH2D("hTMean", "Mean value of waveform per channel vs. time;Event number/10 events;channel number; Mean value",EventRange,Event_Start, Event_Start+EventRange*EventBin, NFEMs*64, 0, NFEMs*64);
   TH2D *hTMax = new TH2D("hTMax", "Max value of waveform per channel vs. time;Event number/10 events;channel number; Mean value",EventRange, Event_Start, Event_Start+EventRange*EventBin, NFEMs*64, 0, NFEMs*64);
   TH2D *hTMin = new TH2D("hTMin", "Min value of waveform per channel vs. time;Event number/10 events;channel number;Mean value",EventRange, Event_Start, Event_Start+EventRange*EventBin,NFEMs*64, 0, NFEMs*64);
   TH2D *hTRMS = new TH2D("hTRMS", "RMS value of waveform per channel vs. time;Event number/10 events;channel number;Mean value",EventRange, Event_Start, Event_Start+EventRange*EventBin,NFEMs*64, 0, NFEMs*64);

  std::cout << "HISTOGRAMS ARE SUCCESSFULLY GENERATED" <<std::endl; 

  double prevEvent = -999;
  int prevSlot = firstFEM;
  int TSample_firstSlot = 0;  //trigger sample at the first physical slot

  double triggerGaps = 0; // Unexpected gap between triggers (missing triggers?)
  double eventJumps = 0; // FEM event number header jumps by >=2
  double totalTriggers = 0; // Total triggers
  double MAX_EVENT = DBL_MIN; //max event number
  std::vector<double> prev_frame(NFEMs, -999);  //a vector to save the frame info of all FEMs from last event;


  
  // Loop over entries, 1 entry is for 1 FEM
  for( int i = 0; i < entries; i++ ){
    inTree->GetEntry(i);

    int id = (int)hinfo->id;
    int slot = (int)hinfo->slot;
    int overflow = (int)hinfo->overflow;
    int full = (int)hinfo->full;
    double nwords = (double)hinfo->nwords;
    double event = (double)hinfo->event;
    double frame = (double)hinfo->frame;
    double checksum = (double)hinfo->checksum;
    int triggerframe = (int)hinfo->triggerframe;
    int triggersample = (int)hinfo->triggersample;
    double wordcount = (double)hinfo->wordcount;
    double mychecksum = (double)hinfo->mychecksum;


//    if(event > MAX_EVENT) MAX_EVENT = event;
    if(i == entries-1) MAX_EVENT = event;
 
    if(hwaveform->size() != 64) std::cout << "Event " << event << ": FEM " << slot << " only has waveforms for " << hwaveform->size() << " channels." <<std::endl; 
    //get the Max, Min, RMS of the waveform
    //loop over all the channels
    for(int k=0; k< (int)hwaveform->size(); k++){ 
	int samplenum = 0;   //size of the waveform.
        double Max, Min, Mean, RMS;
	int channelnumber = (slot - firstFEM)*64 +k;  //global channel number.
	//change unsigned int vector to double vector
	std::vector<unsigned int> mywaveform(hwaveform->at(k).begin(), hwaveform->at(k).end());
	//calculate the size, max, min and RMS of the waveform.
	samplenum = mywaveform.size();
	if(samplenum == 0) std::cout << "Channel " << channelnumber << " has no waveform, SKIP IT!!" << std::endl;
	else{
		Max = *std::max_element(mywaveform.begin(), mywaveform.end());
		Min = *std::min_element(mywaveform.begin(), mywaveform.end());
		RMS = TMath::RMS(mywaveform.begin(), mywaveform.end());
		Mean = TMath::Mean(mywaveform.begin(), mywaveform.end());
        	if(i%10000 == 0 || i>540000)  std::cout << "i =" << i << " Max " << Max << " Min " << Min << " RMS " << RMS << " Mean " << Mean <<std::endl;
	
  		//fill in the histos
		hMean->Fill(channelnumber, Mean);
		hMax->Fill(channelnumber, Max);
		hMin->Fill(channelnumber, Min);
		hRMS->Fill(channelnumber, RMS);
		hTMean->Fill(event,channelnumber, Mean);
		hTMax->Fill(event,channelnumber, Max);
		hTMin->Fill(event,channelnumber, Min);
		hTRMS->Fill(event,channelnumber, RMS);
	}//end of RMS, Min, Max operation
 	hNumSample->Fill(channelnumber, samplenum);
	hTNumSample->Fill(event, channelnumber, samplenum);

    } //end of channel waveforms operation.

    //fill in histos
    hSlot->Fill(slot,(i%NFEMs)+firstFEM);   //compare slot number with manually counted slot number.
    hTSlot->Fill(event,(i%NFEMs)+firstFEM, slot);
    hOverflow->Fill(slot,overflow);
    hTOverflow->Fill(event,slot, overflow);
    hFull->Fill(slot,full);
    hTFull->Fill(event,slot,full);
    hId->Fill(slot,id);
    hTId->Fill(event,slot, id);
    hnwords->Fill(slot,nwords);
    hTnwords->Fill(event,slot, nwords);
    hwordcount->Fill(slot,wordcount);
    hTwordcount->Fill(event,slot, wordcount);
    hdiff_nword->Fill(slot,nwords-wordcount);
    hTdiff_nword->Fill( event,slot, nwords - wordcount);
    hFrames->Fill(frame,slot);
    hTFrames->Fill(event,slot, frame);
    hEvents->Fill(event,slot);
    hTEvents->Fill( event,slot,event);
    hChecksum->Fill(slot,checksum);
    hTChecksum->Fill(event,slot, checksum);
    hMychecksum->Fill(slot,mychecksum);
    hTMychecksum->Fill( event,slot, mychecksum);
    hdiff_checksum->Fill(slot,checksum-mychecksum);
    hTdiff_checksum->Fill(event,slot, checksum - mychecksum);
    hTriggerFrame->Fill(slot,triggerframe);
    hTTriggerFrame->Fill(event,slot, triggerframe);
    hTriggerSample->Fill(slot, triggersample);
    hTTriggerSample->Fill(event,slot, triggersample);
    //the number of frames which go into 1 bin will give you the # of events in 1 second(trigger rate).
    hTTriggerRate->Fill(frame*framesize,slot);

    //if it's the first FEM, don't calculate DeltaFrame and DeltaEvent
    if(i == 0){
	prevEvent = event;
	prev_frame[i%NFEMs]=frame;
	TSample_firstSlot=triggersample;
	continue;	
    } 
   
 
   //comparison between two adjecent FEMs (peridically).
    int deltaFrame = (double)(frame - prev_frame[(i-1)%NFEMs]);
    hDeltaFrame->Fill(deltaFrame,slot);
    hTDeltaFrame->Fill(event,slot, deltaFrame);
    if( deltaFrame != 0 ) totalTriggers++;

    int deltaEvent = event - prevEvent;
    hDeltaEvent->Fill(slot,deltaEvent);
    hTDeltaEvent->Fill(event,slot, deltaEvent);
    if(deltaEvent == 0){ //for entries in the same event
	//if deltaEvent is 0 for the first FEM, something wrong must happen
	if((i%NFEMs) == 0) std::cout<< "ERROR: the first physical slot doesn't start a new event" << std::endl;
	//if this is not the first FEM.
	else{
		hFrameComp[(i-1)%NFEMs]->Fill(prev_frame[0] - frame);
		hTriggerSamComp->Fill(slot, triggersample - TSample_firstSlot);
	}
    }
    else if(slot == firstFEM) TSample_firstSlot=triggersample;
    else std::cout << "ERROR: new event start from a non-first slot!!" << std::endl;  



    //use Frame difference between this event and previous event to calculate the real-time trigger rate.
    if(prev_frame[i%NFEMs] > 0) hFrameTrigger->Fill(slot, 1/((frame-prev_frame[i%NFEMs])*framesize));


    // Look for missed triggers when the frame difference between FEMs is bigger than the tolerance
    if( deltaFrame > 1.05/(triggerRate * frameLength * 0.5e-6) ){ // 1.05 --> add 5% tolerance
      std::cout << "\n\nTRIGGER GAP" << std::endl;
      std::cout << "\tPREVIOUS EVENT" << std::endl;
      for( int j = i - NFEMs; j < i; j++ ) inTree->Show(j);
      std::cout << "\tTHIS EVENT" << std::endl;
      for( int j = i; j < i + NFEMs; j++ ) inTree->Show(j);
      triggerGaps++;
      // If the frame difference occurs in the middle of a crate-readout, FEMs are desynchronized
      if( slot != firstFEM && prevSlot != firstFEM + NFEMs - 1 ){
	std::cout << "DESYNC. Enter anything to continue" << std::endl;
	//std::cin >> aux;
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
      if( slot != firstFEM && prevSlot != firstFEM + NFEMs - 1 ){
	std::cout << "DESYNC. Enter anything to continue" << std::endl;
	//std::cin >> aux;
      }
    }


    //update previous frame, event, and slot.
    prev_frame[i%NFEMs]=frame;
    prevEvent = event;
    prevSlot = slot;
  } // end of loop over N entries

  std::cout << "Histograms filled" << std::endl;

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
  inFile.Close();


  
  //start writing to the output root file.

  rootFile.cd();

  //write the maximum of the event number in the file.
  //for future plot purpose.
  max_event.push_back(MAX_EVENT);
  rootFile.WriteObject(&max_event, "max_event");

  //write the maximum frame number info into file, for future plot purpose
  max_frame.push_back(prev_frame[(entries-1)%NFEMs]);
  rootFile.WriteObject(&max_frame, "max_frame");


  //write plots into file.
  hSlot->Write();
  hTSlot->Write();
  hOverflow->Write();
  hTOverflow->Write();
  hFull->Write();
  hTFull->Write();
  hId->Write();
  hTId->Write();
  hnwords->Write();
  hTnwords->Write();
  hwordcount->Write();
  hTwordcount->Write();
  hdiff_nword->Write();
  hTdiff_nword->Write();
  hFrames->Write();
  hTFrames->Write();
  hDeltaFrame->Write();
  hTDeltaFrame->Write();
  for(int i=0; i< NFEMs-1; i++) hFrameComp[i]->Write();
  hEvents->Write();
  hTEvents->Write();
  hDeltaEvent->Write();
  hTDeltaEvent->Write();
  hChecksum->Write();
  hTChecksum->Write();
  hdiff_checksum->Write();
  hTdiff_checksum->Write();
  hMychecksum->Write();
  hTMychecksum->Write();
  hTriggerFrame->Write();
  hTTriggerFrame->Write();
  hTriggerSample->Write();
  hTTriggerSample->Write();
  hTriggerSamComp->Write();
  hNumSample->Write();
  hTNumSample->Write();
  hMean->Write();
  hTMean->Write();
  hRMS->Write();
  hTRMS->Write();
  hMax->Write();
  hTMax->Write();
  hMin->Write();
  hTMin->Write();
  hTTriggerRate->Write();
  hFrameTrigger->Write();
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
  int status = analyzer( argv[1] ); 
  return status;
}
# endif
