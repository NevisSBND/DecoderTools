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
#include <TLine.h>
#include <TStyle.h>


#include "HeaderInfo.hh"

using namespace std;
TFile *refFile;

//returns the number of words after Huffman compression for a waveform
long int encode_Huffman(std::vector<unsigned int> w){   
	long int nword=1; //number of words in the waveform after Huffman compression
	int H_size = 14; //number of available bits to save data in Huffman word, NU stream specifically
	int word_size=0;
	int size_added=0;  //length just added to the candidate Huffman word 
	for(int i=1; i< w.size() ; i++){
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



int analyzer( const char* runFile ){

  bool is_reference = false;   //if this is for reference run or not.
  int sample_size = 256;  //length of the fake waveform from the WIB
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

  double totEvent = ceil(entries/NFEMs); //calculate the total number events in the run based on entries
  int EventRange = ceil(totEvent/EventBin);  //calcualte num of bins needed for the "event" axis assuming using bin size 'EventBin'.
  std::cout << "There are total " << entries << " entries in the tree"<< std::endl;
  

  std::string outFileName(runFile);
  //std::string outFileName = "tryanalyzer_ana.root";
  std::string outPDFName= outFileName.substr(0, outFileName.find_last_of(".")) + "_ana.pdf";
  outFileName = outFileName.substr(0, outFileName.find_last_of(".")) + "_ana.root";
  TFile rootFile( outFileName.c_str(), "RECREATE" ); 
  TVectorD max_event(1);  //maximum event number in the input data file.

 
 // gROOT->cd();  //make sure all newly created histos are associated with gROOT directory.

  //comparison of frame number between different FEMs within one event
  TH1D* hFrameComp[NFEMs-1];
  for(int i=0; i< NFEMs-1; i++){
	hFrameComp[i] = new TH1D(Form("hFrameComp_%d_%d", firstFEM, firstFEM+i+1), Form("Frame difference between slot %d and %d within 1 event", firstFEM, firstFEM+i+1),10,-5,5);

  }

  //internal comparison of trigger sample within 1 event, expect 0 trigger sample difference.
   TH2D* hTriggerSamComp = new TH2D("hTriggerSamComp","Trigger sample difference between first slot and other slots within 1 event;Slot;trigger sample difference; entries", NFEMs-1, firstFEM+1, firstFEM+NFEMs,100,-500,500);

  //real-time trigger rate calculated from each slot.
  TH2D *hFrameTrigger=new TH2D("hFrameTrigger", "Real-time triggerRate calculated by the frame difference; Slot; Trigger Rate; entries", NFEMs, firstFEM, firstFEM+NFEMs, (int)2*triggerRate, 0, 2*triggerRate);
  

  //simply the distribution of all the variables for difference FEMs.
  TH2D *hSlot = new TH2D("hSlot", "Slot number; Slot; Manually counted Slot number; Entries", NFEMs, firstFEM, firstFEM+NFEMs, NFEMs, firstFEM, firstFEM+NFEMs);
  TH2D *hOverflow = new TH2D("hOverflow", "Overflow distribution;Slot; Overflow flag; Entries", NFEMs, firstFEM, firstFEM+NFEMs, 2, 0,2);
  TH2D *hFull = new TH2D("hFull", "FullDistribution; Slot;Full flag; Entries",NFEMs, firstFEM, firstFEM+NFEMs, 2,0,2);
  TH2D *hId = new TH2D("hId", "ID; Slot;ID; Entries", NFEMs, firstFEM, firstFEM+NFEMs,127,0,127);
  TH2D *hnwords = new TH2D("hnwords", "distribution of number of words; Slot; number of words; entries",NFEMs, firstFEM, firstFEM+NFEMs,10000, 0,100000);
  TH2D *hwordcount = new TH2D("hwordcount", "distribution of manually counted number of words;Slot; wordcount; entries",NFEMs, firstFEM, firstFEM+NFEMs, 10000, 0,100000);
  TH2D *hdiff_nword = new TH2D("hdiff_nword", "number of words - wordcount;Slot; number of words - wordcount; entries", NFEMs, firstFEM, firstFEM+NFEMs,20,-10,10);
  TH2D *hFrames = new TH2D("hFrames", "Event Frame numbers;Slot;Frame;Entries", NFEMs, firstFEM, firstFEM+NFEMs,16778, 0, 16778000);   //largest frame number is 16777215
 //hFrames->SetDirectory(gROOT);  
  TH2D *hDeltaFrame = new TH2D("hDeltaFrame", "Frame difference;Slot;#DeltaFrame;Entries",NFEMs, firstFEM, firstFEM+NFEMs, 9999, 1, 10000);
  TH2D *hEvents = new TH2D("hEvents", "Event numbers;Slot;Event;Entries",NFEMs, firstFEM, firstFEM+NFEMs, EventRange,0, EventRange*EventBin);  //largest event number is 16777215
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
   TH2C *hTSlot = new TH2C("hTSlot", "2D plot of Slot number changing with time(Event number); Event number/10 events;Manually counted Slot number; Slot", EventRange,0, EventRange*EventBin,  NFEMs, firstFEM, firstFEM+NFEMs);
   TH2C *hTOverflow = new TH2C("hTOverflow", "2D plot of overflow vs time(Event number);Event number/10 events; Slot;Overflow",EventRange,0, EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);    //maximum event number: 16778000
   TH2C *hTFull = new TH2C("hTFull", "2D plot of Full flag vs. time(Event number); Event number/10 events;Slot; Full flag", EventRange,0, EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);
   TH2C *hTId = new TH2C("hTId", "2D plot of FEM id vs. time (Event number);Event number/10 events;Slot; FEM ID",EventRange,0, EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);
   TH2D *hTnwords = new TH2D("hTnwords", "2D plot of nwords vs. time(Event number); Event number/10 events;Slot; nwords",EventRange,0, EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);
   TH2D *hTwordcount = new TH2D("hTwordcount", "2D plot of wordcount vs. time(Event number);Event number/10 events;Slot; wordcount",EventRange,0, EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);
   TH2D *hTdiff_nword = new TH2D("hTdiff_nword", "2D plot of (nwords - wordcount) vs. time (Event number);Event number/10 events;Slot; nwords - wordcount", EventRange,0, EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);
   hTdiff_nword->SetAxisRange(-10*EventBin, 10*EventBin, "Z");
   TH2D *hTFrames = new TH2D("hTFrames", "2D plot of event frame number vs. time (Event number); Event number/10 events;Slot; Frame number",EventRange,0, EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);
   TH2D *hTDeltaFrame = new TH2D("hTDeltaFrame", "2D plot of Frame difference;Event number/10 events;Slot; Frame difference", EventRange,0, EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);
   hTDeltaFrame->SetMinimum(1);
   TH2D *hTEvents = new TH2D("hTEvents", "2D plot of Event number vs. time(Event number);Event number/10 events;Slot; Event number",EventRange,0, EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);
   TH2D *hTDeltaEvent = new TH2D("hTDeltaEvent", "2D plot of Event no. difference vs. time (Event number);Event number/10 events;Slot; Event difference",EventRange,0, EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs); 
   TH2D *hTChecksum = new TH2D("hTChecksum", "2D plot of checksum vs. time (Event number);Event number/10 events;Slot;Checksum",EventRange,0, EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);
   TH2D *hTMychecksum = new TH2D("hTMychecksum", "2D plot of manually counted checksum vs. time(Event number);Event number/10 events;Slot; mychecksum",EventRange,0, EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);
   TH2D *hTdiff_checksum = new TH2D("hTdiff_checksum", "2D plot of (checksum - mychecksum) vs. time(Event number);Event number/10 events;Slot;checksum - mychecksum",EventRange,0, EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);
   hTdiff_checksum->SetAxisRange(-4100*EventBin, 4100*EventBin, "Z"); 
   TH2D *hTTriggerRate = new TH2D("hTTriggerRate","2D plot of trigger rate vs. time;Time(s);Slot; Trigger Rate", 10750, 0, 21500, NFEMs, firstFEM, firstFEM+NFEMs);  //based on the largest frame number that can be saved in the Header, and the frame size, the longest time is ~21475 second.
   TH2C *hTTriggerFrame = new TH2C("hTTriggerFrame", "2D plot of Trigger Frame vs. time(Event number);Event number/10 events; Slot;TriggerFrame",EventRange,0, EventRange*EventBin, NFEMs, firstFEM, firstFEM+NFEMs);
   TH2S *hTTriggerSample = new TH2S("hTTriggerSample", "2D plot of Trigger Sample vs. time(Event number);Event number/10 events;Slot;Trigger Sample",EventRange,0, EventRange*EventBin,  NFEMs, firstFEM, firstFEM+NFEMs);
   TH2D *hTNumSample = new TH2D("hTNumSample", "Number of samples per channel vs. time(Event number);Event number/10 events;channel number; Number of Samples",EventRange,0, EventRange*EventBin,  NFEMs*64, 0, NFEMs*64);
   TH2D *hTMean = new TH2D("hTMean", "Mean value of waveform per channel vs. time;Event number/10 events;channel number; Mean value",EventRange,0, EventRange*EventBin, NFEMs*64, 0, NFEMs*64);
   TH2D *hTMax = new TH2D("hTMax", "Max value of waveform per channel vs. time;Event number/10 events;channel number; Mean value",EventRange,0, EventRange*EventBin, NFEMs*64, 0, NFEMs*64);
   TH2D *hTMin = new TH2D("hTMin", "Min value of waveform per channel vs. time;Event number/10 events;channel number;Mean value",EventRange,0, EventRange*EventBin,NFEMs*64, 0, NFEMs*64);
   TH2D *hTRMS = new TH2D("hTRMS", "RMS value of waveform per channel vs. time;Event number/10 events;channel number;Mean value",EventRange,0, EventRange*EventBin,NFEMs*64, 0, NFEMs*64);

  //nwords under Huffman compression for 64 channels (for plotting).
   std::vector<long int> nword_channel_Huff_max(64,0);
   std::vector<long int> nword_channel_Huff_min(64,0);
   std::vector<long int> nword_slot(2,0);  //max nword and min nword for 1 slot. max->0 index, min:1 index.
   std::vector<long int> checksum_channel(64,0); //checksum of 256 samples in 64 channels
   if(!is_reference){
	//if it's not a reference run, delete the memory
	std::vector<long int>().swap(nword_channel_Huff_max);
	std::vector<long int>().swap(nword_channel_Huff_min);
	std::vector<long int>().swap(nword_slot); 
	std::vector<long int>().swap(checksum_channel);
   }
   std::cout << "HISTOGRAMS ARE SUCCESSFULLY GENERATED" <<std::endl; 

  double prevEvent = -999;
  int prevSlot = firstFEM;
  int TSample_firstSlot = 0;  //trigger sample at the first physical slot

  double triggerGaps = 0; // Unexpected gap between triggers (missing triggers?)
  double eventJumps = 0; // FEM event number header jumps by >=2
  double totalTriggers = 0; // Total triggers
  double MAX_EVENT = DBL_MIN; //max event number
  std::vector<double> prev_frame(NFEMs, -999);  //a vector to save the frame info of all FEMs from last event;
  std::vector<long int> ref_checksum_channel;
  std::vector<long int>* pointer_checksum_channel;


  
  //if it's not a reference run, open the reference file.
  if(!is_reference){
	  refFile = new TFile("output_NevisTPCNUXMIT_generator_ref.root", "READ" );
	  if( !refFile->IsOpen() ){
	    std::cerr << "Unable to open the reference root file: output_NevisTPCNUXMIT_generator_ref.root" << std::endl;
	    exit(1);
	  }

	  //grab the checksum of 256 samples in 64 channels.
	  refFile->GetObject("checksum_channel", pointer_checksum_channel);
	  ref_checksum_channel = *pointer_checksum_channel;
  }
  // if it is a reference run, you don't need to use `ref_checksum_channel`, clear its memory and delete the pointer.
  else {
	delete pointer_checksum_channel;
	std::vector<long int>().swap(ref_checksum_channel);
  }



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
		//if it's reference run, calculate the range of nwords with Huffman compression for each channel
		if(is_reference && i==0){
			//get the number of word after Huffman compression
			std::vector<long int> Huffman_size;
			for(int l=0; l<samplenum && l<sample_size; l++){
				//permute the waveform to get range of total number of words
				std::vector<unsigned int> sort_wave(mywaveform.begin()+l, mywaveform.end());	
				sort_wave.insert(sort_wave.end(), mywaveform.begin(), mywaveform.begin()+l);
				Huffman_size.push_back(encode_Huffman(sort_wave));
			}
			//for waveforms with different starting point, get the maximum nwords and minimum nwords under Huffman compression.
			nword_channel_Huff_max[k] = *std::max_element(Huffman_size.begin(),Huffman_size.end());		
		   	nword_channel_Huff_min[k] = *std::min_element(Huffman_size.begin(),Huffman_size.end());	
			
			//get the checksum of selected number (sample_size) of samples in each channel
			checksum_channel[k]=std::accumulate(mywaveform.begin(), mywaveform.begin()+sample_size, 0L);
		}
		else if(!is_reference){
			//calculate the checksum of 'sample_size' ADC words for test run..
			for(int l=0; l<samplenum+1-sample_size; l++){
				long int test_checksum = std::accumulate(mywaveform.begin()+l, mywaveform.begin()+l+sample_size,0L);
				if(test_checksum != ref_checksum_channel[k]) std::cout << "Waveform Different!! \t|| reference: " << ref_checksum_channel[k] << "\t test: " << test_checksum << std::endl;
			}
		}
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
    hFrames->Fill(slot,frame);
    hTFrames->Fill(event,slot, frame);
    hEvents->Fill(slot,event);
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
    int deltaFrame = frame - prev_frame[(i-1)%NFEMs];
    hDeltaFrame->Fill(slot,deltaFrame);
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


 

//***************REFERENCE RUN: write histograms into the output file******************************************.
  if(is_reference){
	  std::cout << "This is reference run!!" << std::endl;
	  std::cout << "Writing histograms into file!!" << std::endl;
	  rootFile.cd();
	  //write the maximum of the event number in the file.
	  max_event[0] = MAX_EVENT;
	  max_event.Write("max_event");

	  //write the range of nword for one slot: 0->Max, 1->Min
	  nword_slot[0]=std::accumulate(nword_channel_Huff_max.begin(), nword_channel_Huff_max.end(),0L);
	  nword_slot[1]=std::accumulate(nword_channel_Huff_min.begin(), nword_channel_Huff_min.end(),0L);
	  rootFile.WriteObject(&nword_slot, "nword_slot");  //write this vector to output file.

	  //write the checksum of certain number of ADC words for every channel (64 channels) 
	  rootFile.WriteObject(&checksum_channel, "checksum_channel");

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
	  TCanvas * c = new TCanvas("reference", "reference plots");
	  gPad->SetLogz(1);  
	  gStyle->SetOptStat(0);
	  gStyle->SetNumberContours(40);  //set the number of colors
	  //gStyle->SetPalette(56);
	  gStyle->SetPalette(kPastel);
	  c->Print(".pdf[");
	  c->cd();
	  hSlot->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hSlot;
	  c->Clear();
	  c->Draw();
	  hTSlot->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hTSlot;
	  c->Clear();
	  c->Draw();
	  hOverflow->Draw("colz");
	  c->Update();
	  c->Print(".pdf"); 
          delete hOverflow; 
	  c->Clear();
	  c->Draw(); 
	  gPad->SetLogz(0);
	  hTOverflow->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
          delete hTOverflow;  
	  c->Clear();
	  c->Draw(); 
	  gPad->SetLogz(1);
	  hFull->Draw("colz");
	  c->Update();
	  c->Print(".pdf"); 
	  delete hFull; 
	  c->Clear();
	  c->Draw(); 
	  gPad->SetLogz(0);
	  hTFull->Draw("colz"); 
	  c->Update();
	  c->Print(".pdf");
	  delete hTFull;
	  c->Clear();
	  c->Draw();
	  gPad->SetLogz(1);
	  hId->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hId;
	  c->Clear();
	  c->Draw();
	  hTId->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hTId;
	  c->Clear();
	  c->Draw();
	  hnwords->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hnwords;
	  c->Clear();
	  c->Draw();
	  hTnwords->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hTnwords;
	  c->Clear();
	  c->Draw();
	  hwordcount->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hwordcount;
	  c->Clear();
	  c->Draw();
	  hTwordcount->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hTwordcount;
	  c->Clear();
	  c->Draw();
	  hdiff_nword->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hdiff_nword;
	  c->Clear();
	  c->Draw();
	  gPad->SetLogz(0);
	  hTdiff_nword->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hTdiff_nword;
	  c->Clear();
	  c->Draw();
	  gPad->SetLogz(1);
	  hFrames->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hFrames;
	  c->Clear();
	  c->Draw();
	  hTFrames->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hTFrames;
	  c->Clear();
	  c->Draw();
	  hDeltaFrame->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hDeltaFrame;
	  c->Clear();
	  c->Draw();
	  hTDeltaFrame->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hTDeltaFrame;
	  c->Clear();
	  c->Draw();
	  hEvents->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hEvents;
	  c->Clear();
	  c->Draw();
	  hTEvents->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hTEvents;
	  c->Clear();
	  c->Draw();
	  hDeltaEvent->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hDeltaEvent;
	  c->Clear();
	  c->Draw();
	  hTDeltaEvent->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hTDeltaEvent;
	  c->Clear();
	  c->Draw();
	  hChecksum->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hChecksum;
	  c->Clear();
	  c->Draw();
	  hTChecksum->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hTChecksum;
	  c->Clear();
	  c->Draw();
	  hdiff_checksum->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hdiff_checksum;
	  c->Clear();
	  c->Draw();
	  gPad->SetLogz(0);
	  hTdiff_checksum->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hTdiff_checksum;
	  c->Clear();
	  c->Draw();
	  gPad->SetLogz(1);
	  hMychecksum->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hMychecksum;
	  c->Clear();
	  c->Draw();
	  hTMychecksum->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hTMychecksum;
	  c->Clear();
	  c->Draw();
	  hTriggerFrame->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hTriggerFrame;
	  c->Clear();
	  c->Draw();
	  hTTriggerFrame->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hTTriggerFrame;
	  c->Clear();
	  c->Draw();
	  hTriggerSample->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hTriggerSample;
	  c->Clear();
	  c->Draw();
	  hTTriggerSample->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hTTriggerSample;
	  c->Clear();
	  c->Draw();
	  hNumSample->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hNumSample;
	  c->Clear();
	  c->Draw();
	  hTNumSample->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hTNumSample;
	  c->Clear();
	  c->Draw();
	  hMean->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hMean;
	  c->Clear();
	  c->Draw();
	  hTMean->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hTMean;
	  c->Clear();
	  c->Draw();
	  hRMS->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hRMS;
	  c->Clear();
	  c->Draw();
	  hTRMS->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hTRMS;
	  c->Clear();
	  c->Draw();
	  hMax->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hMax;
	  c->Clear();
	  c->Draw();
	  hTMax->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hTMax;
	  c->Clear();
	  c->Draw();
	  hMin->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hMin;
	  c->Clear();
	  c->Draw();
	  hTMin->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hTMin;
	  c->Clear();
	  c->Draw();
	  hTTriggerRate->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hTTriggerRate;
	  c->Clear();
	  c->Draw();
	  hFrameTrigger->Draw("colz");
	  c->Update();
	  c->Print(".pdf");
	  delete hFrameTrigger;
	  c->Print(".pdf]");

	  delete c;
	  rootFile.Close();
  }//*************************REFERENCE RUN: end of writing into file*****************************************
  else{


	/*
	  //open the reference root file.
	  TFile *refFile = new TFile("output_NevisTPCNUXMIT_generator_ref.root", "READ" );
	  if( !refFile->IsOpen() ){
	    std::cerr << "Unable to open the reference root file: output_NevisTPCNUXMIT_generator_ref.root" << std::endl;
	    exit(1);
	  }
	 */
	  //define pointers to grab info from reference root file
	  TH2D *hRef;
	  TH2D *hTRef;
	  TVectorD* vRef = (TVectorD*)refFile->Get("max_event"); //(*vRef)[0] will be the largest event number in reference run 
	  double ref_event = ((*vRef))[0];
	  std::cout << "ref_event is " << ref_event << std::endl;
	  std::vector<long int>* ref_Huff;
	  refFile->GetObject("nword_slot",ref_Huff);
	  std::vector<long int> nword_Huff = *ref_Huff;
	  std::cout<< "Max nword is " << (*ref_Huff)[0] << " || Min nword is " << (*ref_Huff)[1] <<std::endl; 


	  double z_max, z_min, z_maxT, z_minT; //maximum and minimum of z values for two plots

	  rootFile.cd();
	  TCanvas *c = new TCanvas("c", "c");
	  
	  gStyle->SetOptStat("ou");   //show overflow and undeflow in stat box
	  gStyle->SetNumberContours(40);
	  gStyle->SetPalette(kPastel);
	  

	  c->Print((outPDFName+"[").c_str());
	  //draw Frame info.
	  hRef = (TH2D*)refFile->Get("hFrames");
	  hTRef = (TH2D*)refFile->Get("hTFrames");
	  hRef->Scale(MAX_EVENT/ref_event); //scale the reference plots with the same total number of events
	  hTRef->SetAxisRange(0,EventRange*EventBin, "X");   //cut the X range so that reference plot and test run plot match.
	  z_max = hRef->GetMaximum();
	  z_min = hRef->GetMinimum();
	  z_maxT = hTRef->GetMaximum();
	  z_minT = hTRef->GetMinimum();
	  c->Divide(1, 2);
	  c->cd(1);
	  gPad->SetLogz(1);    //set log scale on z.
	  if(hFrames->GetMaximum() > z_max) z_max = hFrames->GetMaximum();
	  if(hFrames->GetMinimum() < z_min) z_min = hFrames->GetMinimum();
	  hFrames->SetAxisRange(z_min-1, z_max+1, "Z");
	  hFrames->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);    //set log scale on z.
	  hRef->SetAxisRange(z_min-1, z_max+1, "Z"); 
	  hRef->SetTitle("reference run");
	  hRef->Draw("colz");
	  c->Update();
	  //std::cout << "GetEntries is " << hFrames->GetEntries() << " " << hRef->GetEntries()<<std::endl;
	  //rootFile.cd();
	  c->Write("cFrame");
	  c->Print(outPDFName.c_str());
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);    //set log scale on z.
	  if(hTFrames->GetMaximum() > z_maxT) z_maxT = hTFrames->GetMaximum();
	  if(hTFrames->GetMinimum() < z_minT) z_minT = hTFrames->GetMinimum();
	  hTFrames->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTFrames->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);    //set log scale on z.
	  hTRef->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTRef->SetTitle("reference run");
	  hTRef->Draw("colz");
	  c->Update();
	  c->Write("cTFrame");
	  c->Print(outPDFName.c_str());
	  delete hFrames; delete hTFrames;


	  hRef = (TH2D*)refFile->Get("hEvents");
	  hTRef = (TH2D*)refFile->Get("hTEvents");
	  hRef->Scale(MAX_EVENT/ref_event); //scale the reference plots with the same total number of events
	  hTRef->SetAxisRange(0,EventRange*EventBin, "X");   //cut the X range so that reference plot and test run plot match.
	  z_max = hRef->GetMaximum();
	  z_min = hRef->GetMinimum();
	  z_maxT = hTRef->GetMaximum();
	  z_minT = hTRef->GetMinimum(); 
	  c->Clear();
	  c->Draw();
	  c->Divide(1, 2);
	  c->cd(1);
	  gPad->SetLogz(1);    //set log scale on z.
	  if(hEvents->GetMaximum() > z_max) z_max=hEvents->GetMaximum();
	  if(hEvents->GetMinimum() < z_min) z_min=hEvents->GetMinimum();
	  hEvents->SetAxisRange(z_min-1, z_max+1, "Z");
	  hEvents->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);    //set log scale on z.
	  hRef->SetAxisRange(z_min-1, z_max+1, "Z");
	  hRef->SetTitle("reference run");
	  hRef->Draw("colz");
	  c->Update();
	  c->Write("cEvents");
	  c->Print(outPDFName.c_str());
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);    //set log scale on z.
	  if(hTEvents->GetMaximum() > z_maxT) z_maxT = hTEvents->GetMaximum();
	  if(hTEvents->GetMinimum() < z_minT) z_minT = hTEvents->GetMinimum();
	  hTEvents->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTEvents->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);    //set log scale on z.
	  hTRef->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTRef->SetTitle("reference run");
	  hTRef->Draw("colz");
	  c->Update();
	  c->Write("cTEvent");
	  c->Print(outPDFName.c_str());
	  delete hEvents; delete hTEvents;

	  
	  hRef = (TH2D*)refFile->Get("hSlot");
	  hTRef = (TH2D*)refFile->Get("hTSlot");
	  hRef->Scale(MAX_EVENT/ref_event); //scale the reference plots with the same total number of events
	  hTRef->Scale(1/(double)EventBin);
	  hTRef->SetAxisRange(0,EventRange*EventBin, "X");   //cut the X range so that reference plot and test run plot match.
	  z_max = hRef->GetMaximum();
	  z_min = hRef->GetMinimum();
	  z_maxT = hTRef->GetMaximum();
	  z_minT = hTRef->GetMinimum();
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);    //set log scale on z.
	  if(hSlot->GetMaximum() > z_max) z_max=hSlot->GetMaximum();
	  if(hSlot->GetMinimum() < z_min) z_min=hSlot->GetMinimum();
	  hSlot->SetAxisRange(z_min-1, z_max+1, "Z");
	  hSlot->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);    //set log scale on z.
	  hRef->SetAxisRange(z_min-1, z_max+1, "Z");
	  hRef->SetTitle("reference run");
	  hRef->Draw("colz");
	  c->Update();
	  c->Write("cSlot");
	  c->Print(outPDFName.c_str());
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);    //set log scale on z.
	  hTSlot->Scale(1/(double)EventBin);
	  if(hTSlot->GetMaximum() > z_maxT) z_maxT = hTSlot->GetMaximum();
	  if(hTSlot->GetMinimum() < z_minT) z_minT = hTSlot->GetMinimum();
	  hTSlot->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTSlot->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);    //set log scale on z.
	  hTRef->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTRef->SetTitle("reference run");
	  hTRef->Draw("colz");
	  c->Update();
	  c->Write("cTSlot");
	  c->Print(outPDFName.c_str());
	  delete hSlot; delete hTSlot;
	  
	  hRef = (TH2D*)refFile->Get("hOverflow");
	  hTRef = (TH2D*)refFile->Get("hTOverflow");
	  hRef->Scale(MAX_EVENT/ref_event); //scale the reference plots with the same total number of events
	  hTRef->SetAxisRange(0,EventRange*EventBin, "X");   //cut the X range so that reference plot and test run plot match.
	  z_max = hRef->GetMaximum();
	  z_min = hRef->GetMinimum();
	  z_maxT = hTRef->GetMaximum();
	  z_minT = hTRef->GetMinimum();
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);    //set log scale on z.
	  if(hOverflow->GetMaximum() > z_max) z_max=hOverflow->GetMaximum();
	  if(hOverflow->GetMinimum() < z_min) z_min = hOverflow->GetMinimum();
	  hOverflow->SetAxisRange(z_min-1, z_max+1, "Z");
	  hOverflow->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);    //set log scale on z.
	  hRef->SetAxisRange(z_min-1, z_max+1, "Z");
	  hRef->SetTitle("reference run");
	  hRef->Draw("colz");
	  c->Update();
	  c->Write("cOverflow");
	  c->Print(outPDFName.c_str());
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  gPad->SetLogz(0);
	  c->cd(1);
	  if(hTOverflow->GetMaximum() > z_maxT) z_maxT = hTOverflow->GetMaximum();
	  if(hTOverflow->GetMinimum() < z_minT) z_minT = hTOverflow->GetMinimum();
	  hTOverflow->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTOverflow->Draw("colz");
	  c->cd(2);
	  hTRef->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTRef->SetTitle("reference run");
	  hTRef->Draw("colz");
	  c->Update();
	  c->Write("cTOverflow");
	  c->Print(outPDFName.c_str());
	  delete hOverflow; delete hTOverflow;


	  hRef=(TH2D*)refFile->Get("hFull");
	  hTRef=(TH2D*)refFile->Get("hTFull"); 
	  hRef->Scale(MAX_EVENT/ref_event); //scale the reference plots with the same total number of events
	  hTRef->SetAxisRange(0,EventRange*EventBin, "X");   //cut the X range so that reference plot and test run plot match.
	  z_max = hRef->GetMaximum();
	  z_min = hRef->GetMinimum();
	  z_maxT = hTRef->GetMaximum();
	  z_minT = hTRef->GetMinimum(); 
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hFull->GetMaximum() > z_max) z_max=hFull->GetMaximum();
	  if(hFull->GetMinimum() < z_min) z_min = hFull->GetMinimum();
	  hFull->SetAxisRange(z_min-1, z_max+1, "Z");
	  hFull->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hRef->SetAxisRange(z_min-1, z_max+1, "Z");
	  hRef->SetTitle("reference run");
	  hRef->Draw("colz");
	  c->Update();
	  c->Write("cFull");
	  c->Print(outPDFName.c_str());
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  gPad->SetLogz(0);
	  c->cd(1);
	  if(hTFull->GetMaximum() > z_maxT) z_maxT=hTFull->GetMaximum();
	  if(hTFull->GetMinimum() < z_minT) z_minT=hTFull->GetMinimum();
	  hTFull->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTFull->Draw("colz");
	  c->cd(2);
	  hTRef->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTRef->SetTitle("reference run");
	  hTRef->Draw("colz");
	  c->Update();
	  c->Write("cTFull");
	  c->Print(outPDFName.c_str());
	  delete hFull; delete hTFull;


	  hRef=(TH2D*)refFile->Get("hId");
	  hTRef=(TH2D*)refFile->Get("hTId");
	  hRef->Scale(MAX_EVENT/ref_event); //scale the reference plots with the same total number of events
	  hTRef->SetAxisRange(0,EventRange*EventBin, "X");   //cut the X range so that reference plot and test run plot match.
	  z_max = hRef->GetMaximum();
	  z_min = hRef->GetMinimum();
	  z_maxT = hTRef->GetMaximum();
	  z_minT = hTRef->GetMinimum();
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hId->GetMaximum() > z_max) z_max = hId->GetMaximum();
	  if(hId->GetMinimum() < z_min) z_min = hId->GetMinimum();
	  hId->SetAxisRange(z_min-1, z_max+1, "Z");
	  hId->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hRef->SetAxisRange(z_min-1, z_max+1, "Z");
	  hRef->SetTitle("reference run");
	  hRef->Draw("colz");
	  c->Update();
	  c->Write("cId");
	  c->Print(outPDFName.c_str());
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hTId->GetMaximum() > z_maxT) z_maxT=hTId->GetMaximum();
	  if(hTId->GetMinimum() < z_minT) z_minT = hTId->GetMinimum();
	  hTId->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTId->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hTRef->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTRef->SetTitle("reference run");
	  hTRef->Draw("colz");
	  c->Update();
	  c->Write("cTId");
	  c->Print(outPDFName.c_str());
	  delete hId; delete hTId;


	  hRef=(TH2D*)refFile->Get("hnwords");
	  hTRef=(TH2D*)refFile->Get("hTnwords");
	  hRef->Scale(MAX_EVENT/ref_event); //scale the reference plots with the same total number of events
	  hTRef->SetAxisRange(0,EventRange*EventBin, "X");   //cut the X range so that reference plot and test run plot match.
	  z_max = hRef->GetMaximum();
	  z_min = hRef->GetMinimum();
	  z_maxT = hTRef->GetMaximum();
	  z_minT = hTRef->GetMinimum();
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hnwords->GetMaximum() > z_max) z_max=hnwords->GetMaximum();
	  if(hnwords->GetMinimum() < z_min) z_min = hnwords->GetMinimum();  
	  hnwords->SetAxisRange(z_min-1, z_max+1, "Z");
	  hnwords->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hRef->SetAxisRange(z_min-1, z_max+1, "Z");
	  hRef->SetTitle("reference run");
	  hRef->Draw("colz");
	  c->Update();
	  c->Write("cnwords");
	  c->Print(outPDFName.c_str());


	  //nword distribution plots for every slot, with possible nword range shown.
	  for(int i=1; i<= NFEMs; i++){
		c->Clear();
		c->Draw();
		TH1D* h_nwordPerSlot = hnwords->ProjectionY("",i,i);
		h_nwordPerSlot->SetTitle(Form("nword distribution for slot %d; nwords;entries", firstFEM+i-1));
		h_nwordPerSlot->Draw("hist[");
		c->Update();
		TLine* l_max = new TLine(nword_Huff[0], gPad->GetUymin(), nword_Huff[0], gPad->GetUymax());
		TLine* l_min = new TLine(nword_Huff[1], gPad->GetUymin(), nword_Huff[1], gPad->GetUymax());
		l_max->SetLineColor(kRed);
		l_min->SetLineColor(kRed);
		l_max->Draw();
		l_min->Draw();
		c->Update();
		c->Write(Form("cnword_slot%d",firstFEM+i-1));
		c->Print(outPDFName.c_str());
		delete h_nwordPerSlot;
		delete l_max;
		delete l_min;
	  }



	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hTnwords->GetMaximum() > z_maxT) z_maxT=hTnwords->GetMaximum();
	  if(hTnwords->GetMinimum() < z_minT) z_minT=hTnwords->GetMinimum();
	  hTnwords->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTnwords->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hTRef->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTRef->SetAxisRange(0,EventRange*EventBin, "X");
	  hTRef->SetTitle("reference run");
	  hTRef->Draw("colz");
	  c->Update();
	  c->Write("cTnwords");
	  c->Print(outPDFName.c_str());
	  delete hnwords; delete hTnwords;
	 
	  hRef=(TH2D*)refFile->Get("hwordcount");
	  hTRef=(TH2D*)refFile->Get("hTwordcount");
	  hRef->Scale(MAX_EVENT/ref_event); //scale the reference plots with the same total number of events
	  hTRef->SetAxisRange(0,EventRange*EventBin, "X");   //cut the X range so that reference plot and test run plot match.
	  z_max = hRef->GetMaximum();
	  z_min = hRef->GetMinimum();
	  z_maxT = hTRef->GetMaximum();
	  z_minT = hTRef->GetMinimum();
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hwordcount->GetMaximum() > z_max) z_max=hwordcount->GetMaximum();
	  if(hwordcount->GetMinimum() < z_min) z_min=hwordcount->GetMinimum();
	  hwordcount->SetAxisRange(z_min-1, z_max+1, "Z");
	  hwordcount->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hRef->SetAxisRange(z_min-1, z_max+1, "Z");
	  hRef->SetTitle("reference run");
	  hRef->Draw("colz");
	  c->Update();
	  c->Write("cwordcount");
	  c->Print(outPDFName.c_str());
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hTwordcount->GetMaximum() > z_maxT) z_maxT=hTwordcount->GetMaximum();
	  if(hTwordcount->GetMinimum() < z_minT) z_minT=hTwordcount->GetMinimum();
	  hTwordcount->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTwordcount->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hTRef->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTRef->SetAxisRange(0,EventRange*EventBin, "X");
	  hTRef->SetTitle("reference run");
	  hTRef->Draw("colz");
	  c->Update();
	  c->Write("cTwordcount");
	  c->Print(outPDFName.c_str());
	  delete hwordcount; delete hTwordcount;

	  hRef=(TH2D*)refFile->Get("hdiff_nword");
	  hTRef=(TH2D*)refFile->Get("hTdiff_nword");
	  hRef->Scale(MAX_EVENT/ref_event); //scale the reference plots with the same total number of events
	  hTRef->SetAxisRange(0,EventRange*EventBin, "X");   //cut the X range so that reference plot and test run plot match.
	  z_max = hRef->GetMaximum();
	  z_min = hRef->GetMinimum();
	  z_maxT = hTRef->GetMaximum();
	  z_minT = hTRef->GetMinimum();
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hdiff_nword->GetMaximum() > z_max) z_max=hdiff_nword->GetMaximum();
	  if(hdiff_nword->GetMinimum() < z_min) z_min=hdiff_nword->GetMinimum();
	  hdiff_nword->SetAxisRange(z_min-1, z_max+1, "Z");
	  hdiff_nword->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hRef->SetAxisRange(z_min-1, z_max+1, "Z");
	  hRef->SetTitle("reference run");
	  hRef->Draw("colz");
	  c->Update();
	  c->Write("cdiff_nwords");
	  c->Print(outPDFName.c_str());
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  gPad->SetLogz(0);
	  c->cd(1);
	  hTdiff_nword->Draw("colz");
	  c->cd(2);
	  hTRef->SetAxisRange(-10*EventBin, 10*EventBin, "Z");
	  hTRef->SetTitle("reference run");
	  hTRef->Draw("colz");
	  c->Update();
	  c->Write("cTdiff_nwords");
	  c->Print(outPDFName.c_str());
	  delete hdiff_nword; delete hTdiff_nword;


	  hRef=(TH2D*)refFile->Get("hDeltaFrame");
	  hTRef=(TH2D*)refFile->Get("hTDeltaFrame");
	  hRef->Scale(MAX_EVENT/ref_event); //scale the reference plots with the same total number of events
	  hTRef->SetAxisRange(0,EventRange*EventBin, "X");   //cut the X range so that reference plot and test run plot match.
	  z_max = hRef->GetMaximum();
	  z_min = hRef->GetMinimum();
	  z_maxT = hTRef->GetMaximum();
	  z_minT = hTRef->GetMinimum();
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hDeltaFrame->GetMaximum() > z_max) z_max=hDeltaFrame->GetMaximum();
	  if(hDeltaFrame->GetMinimum() < z_min) z_min=hDeltaFrame->GetMinimum();
	  hDeltaFrame->SetAxisRange(z_min-1, z_max+1, "Z");
	  hDeltaFrame->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hRef->SetAxisRange(z_min-1, z_max+1, "Z");
	  hRef->SetTitle("reference run");
	  hRef->Draw("colz");
	  c->Update();
	  c->Write("cDeltaFrame");
	  c->Print(outPDFName.c_str());
	  //frame difference plot between slots within 1 event
	  for(int i =0; i<NFEMs-1; i++){
		c->Clear();
		c->Draw();
		gPad->SetLogy(1);
		hFrameComp[i]->Draw("hist[");
		c->Update();
		c->Write(hFrameComp[i]->GetTitle());
		c->Print(outPDFName.c_str());
		delete hFrameComp[i];
	  }
	  //end of plotting frame difference within 1 event
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hTDeltaFrame->GetMaximum() > z_maxT) z_maxT=hTDeltaFrame->GetMaximum();
	  if(hTDeltaFrame->GetMinimum() > z_minT) z_minT=hTDeltaFrame->GetMinimum();
	  hTDeltaFrame->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTDeltaFrame->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hTRef->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  //hTRef->SetMinimum(1);  //set the z range to be the same as hTDeltaFrame
	  hTRef->SetTitle("reference run");
	  hTRef->Draw("colz");
	  c->Update();
	  c->Write("cTDeltaFrame");
	  c->Print(outPDFName.c_str());
	  delete hDeltaFrame; delete hTDeltaFrame;



	  hRef=(TH2D*)refFile->Get("hDeltaEvent");
	  hTRef=(TH2D*)refFile->Get("hTDeltaEvent");
	  hRef->Scale(MAX_EVENT/ref_event); //scale the reference plots with the same total number of events
	  hTRef->SetAxisRange(0,EventRange*EventBin, "X");   //cut the X range so that reference plot and test run plot match.
	  z_max = hRef->GetMaximum();
	  z_min = hRef->GetMinimum();
	  z_maxT = hTRef->GetMaximum();
	  z_minT = hTRef->GetMinimum();
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hDeltaEvent->GetMaximum() > z_max) z_max=hDeltaEvent->GetMaximum();
	  if(hDeltaEvent->GetMinimum() < z_min) z_min=hDeltaEvent->GetMinimum();
	  hDeltaEvent->SetAxisRange(z_min-1, z_max+1, "Z");
	  hDeltaEvent->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hRef->SetAxisRange(z_min-1, z_max+1, "Z");
	  hRef->SetTitle("reference run");
	  hRef->Draw("colz");
	  c->Update();
	  c->Write("cDeltaEvent");
	  c->Print(outPDFName.c_str());
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hTDeltaEvent->GetMaximum() > z_maxT) z_maxT=hTDeltaEvent->GetMaximum();
	  if(hTDeltaEvent->GetMinimum() < z_minT) z_minT=hTDeltaEvent->GetMinimum();
	  hTDeltaEvent->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTDeltaEvent->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hTRef->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTRef->SetTitle("reference run");
	  hTRef->Draw("colz");
	  c->Update();
	  c->Write("cTDeltaEvent");
	  c->Print(outPDFName.c_str());
	  delete hDeltaEvent; delete hTDeltaEvent;



	 
	  hRef=(TH2D*)refFile->Get("hTTriggerRate");
	  hRef->SetAxisRange(0, ceil(MAX_EVENT/triggerRate), "X");  
	  z_max = hRef->GetMaximum();
	  z_min = hRef->GetMinimum();
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hTTriggerRate->GetMaximum() > z_max) z_max=hTTriggerRate->GetMaximum();
	  if(hTTriggerRate->GetMinimum() < z_min) z_min=hTTriggerRate->GetMinimum();
	  hTTriggerRate->SetAxisRange(z_min, z_max+1, "Z");
	  hTTriggerRate->SetAxisRange(0, ceil(MAX_EVENT/triggerRate), "X");
	  hTTriggerRate->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hRef->SetAxisRange(z_min, z_max+1, "Z");
	  hRef->SetTitle("reference run");
	  hRef->Draw("colz");
	  c->Update();
	  c->Write("cTriggerRate");
	  c->Print(outPDFName.c_str());
	  delete hTTriggerRate;

	  hRef=(TH2D*)refFile->Get("hFrameTrigger");
	  hRef->Scale(MAX_EVENT/ref_event);
	  z_max = hRef->GetMaximum();
	  z_min = hRef->GetMinimum();
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hFrameTrigger->GetMaximum() > z_max) z_max=hFrameTrigger->GetMaximum();
	  if(hFrameTrigger->GetMinimum() < z_min) z_min=hFrameTrigger->GetMinimum();
	  hFrameTrigger->SetAxisRange(z_min-1, z_max+1, "Z");
	  hFrameTrigger->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hRef->SetAxisRange(z_min-1, z_max+1, "Z");
	  hRef->SetTitle("reference run");
	  hRef->Draw("colz");
	  c->Update();
	  c->Write("cFrameTrigger");
	  c->Print(outPDFName.c_str());
	  delete hFrameTrigger;

	  hRef=(TH2D*)refFile->Get("hTriggerFrame");
	  hTRef=(TH2D*)refFile->Get("hTTriggerFrame");
	  hRef->Scale(MAX_EVENT/ref_event); //scale the reference plots with the same total number of events
	  hTRef->SetAxisRange(0,EventRange*EventBin, "X");   //cut the X range so that reference plot and test run plot match.
	  z_max = hRef->GetMaximum();
	  z_min = hRef->GetMinimum();
	  z_maxT = hTRef->GetMaximum();
	  z_minT = hTRef->GetMinimum();
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hTriggerFrame->GetMaximum() > z_max) z_max=hTriggerFrame->GetMaximum();
	  if(hTriggerFrame->GetMinimum() < z_min) z_min=hTriggerFrame->GetMinimum();
	  hTriggerFrame->SetAxisRange(z_min-1, z_max+1, "Z");
	  hTriggerFrame->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hRef->SetAxisRange(z_min-1, z_max+1, "Z");
	  hRef->SetTitle("reference run");
	  hRef->Draw("colz");
	  c->Update();
	  c->Write("cTriggerFrame");
	  c->Print(outPDFName.c_str());
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hTTriggerFrame->GetMaximum() > z_maxT) z_maxT=hTTriggerFrame->GetMaximum();
	  if(hTTriggerFrame->GetMinimum() < z_minT) z_minT=hTTriggerFrame->GetMinimum();
	  hTTriggerFrame->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTTriggerFrame->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hTRef->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTRef->SetAxisRange(0,EventRange*EventBin, "X");
	  hTRef->SetTitle("reference run");
	  hTRef->Draw("colz");
	  c->Update();
	  c->Write("cTTriggerFrame");
	  c->Print(outPDFName.c_str());
	  delete hTriggerFrame; delete hTTriggerFrame;


	  hRef=(TH2D*)refFile->Get("hTriggerSample");
	  hTRef=(TH2D*)refFile->Get("hTTriggerSample");
	  hRef->Scale(MAX_EVENT/ref_event); //scale the reference plots with the same total number of events
	  hTRef->SetAxisRange(0,EventRange*EventBin, "X");   //cut the X range so that reference plot and test run plot match.
	  z_max = hRef->GetMaximum();
	  z_min = hRef->GetMinimum();
	  z_maxT = hTRef->GetMaximum();
	  z_minT = hTRef->GetMinimum();
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hTriggerSample->GetMaximum() > z_max) z_max=hTriggerSample->GetMaximum();
	  if(hTriggerSample->GetMinimum() < z_min) z_min=hTriggerSample->GetMinimum();
	  hTriggerSample->SetAxisRange(z_min-1, z_max+1, "Z");
	  hTriggerSample->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hRef->SetAxisRange(z_min-1, z_max+1, "Z");
	  hRef->SetTitle("reference run");
	  hRef->Draw("colz");
	  c->Update();
	  c->Write("cTriggerSample");
	  c->Print(outPDFName.c_str());

	 //trigger sample distribution for each slot
	  for(int i=1; i<= NFEMs; i++){
		TH1D* htriggersample_Projection = hTriggerSample->ProjectionY("", i,i);
		htriggersample_Projection->SetTitle(Form("trigger sample distribution at slot %d; trigger sample: entries", firstFEM+i-1));
		c->Clear();
		c->Draw();
		htriggersample_Projection->Draw("hist[");
		c->Update();
		c->Write(Form("cTriggerSample_slot%d",  firstFEM+i-1));
		c->Print(outPDFName.c_str());
		delete htriggersample_Projection;
	  }
	  //end of trigger sample distribution for each slot.
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hTTriggerSample->GetMaximum() > z_maxT) z_maxT=hTTriggerSample->GetMaximum();
	  if(hTTriggerSample->GetMinimum() < z_minT) z_minT=hTTriggerSample->GetMinimum();
	  hTTriggerSample->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTTriggerSample->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hTRef->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTRef->SetTitle("reference run");
	  hTRef->Draw("colz");
	  c->Update();
	  c->Write("cTTriggerSample");
	  c->Print(outPDFName.c_str());
	  delete hTriggerSample; delete hTTriggerSample;

	  
	  hRef=(TH2D*)refFile->Get("hTriggerSamComp");
	  hRef->Scale(MAX_EVENT/ref_event); //scale the reference plots with the same total number of events
	  z_max = hRef->GetMaximum();
	  z_min = hRef->GetMinimum();
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hTriggerSamComp->GetMaximum() > z_max) z_max=hTriggerSamComp->GetMaximum();
	  if(hTriggerSamComp->GetMinimum() < z_min) z_min=hTriggerSamComp->GetMinimum();
	  hTriggerSamComp->SetAxisRange(z_min-1, z_max+1, "Z");
	  hTriggerSamComp->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hRef->SetAxisRange(z_min-1, z_max+1, "Z");
	  hRef->SetTitle("reference run");
	  hRef->Draw("colz");
	  c->Update();
	  c->Write("cTriggerSamComp");
	  c->Print(outPDFName.c_str());


	  hRef=(TH2D*)refFile->Get("hChecksum");
	  hTRef=(TH2D*)refFile->Get("hTChecksum");
	  hRef->Scale(MAX_EVENT/ref_event); //scale the reference plots with the same total number of events
	  hTRef->SetAxisRange(0,EventRange*EventBin, "X");   //cut the X range so that reference plot and test run plot match.
	  z_max = hRef->GetMaximum();
	  z_min = hRef->GetMinimum();
	  z_maxT = hTRef->GetMaximum();
	  z_minT = hTRef->GetMinimum();
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hChecksum->GetMaximum() > z_max) z_max=hChecksum->GetMaximum();
	  if(hChecksum->GetMinimum() < z_min) z_min=hChecksum->GetMinimum();
	  hChecksum->SetAxisRange(z_min-1, z_max+1, "Z");
	  hChecksum->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hRef->SetAxisRange(z_min-1, z_max+1, "Z");
	  hRef->SetTitle("reference run");
	  hRef->Draw("colz");
	  c->Update();
	  c->Write("cChecksum");
	  c->Print(outPDFName.c_str());
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hTChecksum->GetMaximum() > z_maxT) z_maxT=hTChecksum->GetMaximum();
	  if(hTChecksum->GetMinimum() < z_minT) z_minT=hTChecksum->GetMinimum();
	  hTChecksum->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTChecksum->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hTRef->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTRef->SetAxisRange(0,EventRange*EventBin, "X");
	  hTRef->SetTitle("reference run");
	  hTRef->Draw("colz");
	  c->Update();
	  c->Write("cTChecksum");
	  c->Print(outPDFName.c_str());
	  delete hChecksum; delete hTChecksum;


	  hRef=(TH2D*)refFile->Get("hMychecksum");
	  hTRef=(TH2D*)refFile->Get("hTMychecksum");
	  hRef->Scale(MAX_EVENT/ref_event); //scale the reference plots with the same total number of events
	  hTRef->SetAxisRange(0,EventRange*EventBin, "X");   //cut the X range so that reference plot and test run plot match.
	  z_max = hRef->GetMaximum();
	  z_min = hRef->GetMinimum();
	  z_maxT = hTRef->GetMaximum();
	  z_minT = hTRef->GetMinimum();
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hMychecksum->GetMaximum() > z_max) z_max=hMychecksum->GetMaximum();
	  if(hMychecksum->GetMinimum() < z_min) z_min=hMychecksum->GetMinimum();
	  hMychecksum->SetAxisRange(z_min-1, z_max+1, "Z");
	  hMychecksum->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hRef->SetAxisRange(z_min-1, z_max+1, "Z");
	  hRef->SetTitle("reference run");
	  hRef->Draw("colz");
	  c->Update();
	  c->Write("cMychecksum");
	  c->Print(outPDFName.c_str());
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hTMychecksum->GetMaximum() > z_maxT) z_maxT=hTMychecksum->GetMaximum();
	  if(hTMychecksum->GetMinimum() < z_minT) z_minT=hTMychecksum->GetMinimum();
	  hTMychecksum->SetAxisRange(z_minT-1,z_maxT+1,"Z");
	  hTMychecksum->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hTRef->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTRef->SetTitle("reference run");
	  hTRef->Draw("colz");
	  c->Update();
	  c->Write("cTMychecksum");
	  c->Print(outPDFName.c_str());
	  delete hMychecksum; delete hTMychecksum;


	  hRef=(TH2D*)refFile->Get("hdiff_checksum");
	  hTRef=(TH2D*)refFile->Get("hTdiff_checksum");
	  hRef->Scale(MAX_EVENT/ref_event); //scale the reference plots with the same total number of events
	  hTRef->SetAxisRange(0,EventRange*EventBin, "X");   //cut the X range so that reference plot and test run plot match.
	  z_max = hRef->GetMaximum();
	  z_min = hRef->GetMinimum();
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hdiff_checksum->GetMaximum() > z_max ) z_max=hdiff_checksum->GetMaximum();
	  if(hdiff_checksum->GetMinimum() < z_min) z_min=hdiff_checksum->GetMinimum();
	  hdiff_checksum->SetAxisRange(z_min-1, z_max+1, "Z");
	  hdiff_checksum->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hRef->SetAxisRange(z_min-1, z_max+1, "Z");
	  hRef->SetTitle("reference run");
	  hRef->Draw("colz");
	  c->Update();
	  c->Write("cdiff_checksum");
	  c->Print(outPDFName.c_str());
	  //zomming in version of diff_checksum
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  hdiff_checksum->SetAxisRange(-100,100, "Y");
	  hdiff_checksum->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hRef->SetAxisRange(-100,100, "Y");
	  hRef->Draw("colz");
	  c->Update();
	  c->Write("cdiff_checksum_zoom");
	  c->Print(outPDFName.c_str());
	  //end of zooming in plot
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  gPad->SetLogz(0);
	  c->cd(1);
	  hTdiff_checksum->Draw("colz");
	  c->cd(2);
	  hTRef->SetAxisRange(-4100*EventBin, 4100*EventBin, "Z");
	  hTRef->SetTitle("reference run");
	  hTRef->Draw("colz");
	  c->Update();
	  c->Write("cTdiff_checksum");
	  c->Print(outPDFName.c_str());
	  delete hdiff_checksum; delete hTdiff_checksum;


	  hRef=(TH2D*)refFile->Get("hRMS");
	  hTRef=(TH2D*)refFile->Get("hTRMS");
	  hRef->Scale(MAX_EVENT/ref_event); //scale the reference plots with the same total number of events
	  hTRef->SetAxisRange(0,EventRange*EventBin, "X");   //cut the X range so that reference plot and test run plot match.
	  z_max = hRef->GetMaximum();
	  z_min = hRef->GetMinimum();
	  z_maxT = hTRef->GetMaximum();
	  z_minT = hTRef->GetMinimum();
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hRMS->GetMaximum() > z_max ) z_max=hRMS->GetMaximum();
	  if(hRMS->GetMinimum() < z_min ) z_min=hRMS->GetMinimum();
	  hRMS->SetAxisRange(z_min-1, z_max+1, "Z");
	  hRMS->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hRef->SetAxisRange(z_min-1, z_max+1, "Z");
	  hRef->SetTitle("reference run");
	  hRef->Draw("colz");
	  c->Update();
	  c->Write("cRMS");
	  c->Print(outPDFName.c_str());
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hTRMS->GetMaximum() > z_maxT) z_maxT=hTRMS->GetMaximum();
	  if(hTRMS->GetMinimum() < z_minT) z_minT=hTRMS->GetMinimum();
	  hTRMS->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTRMS->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hTRef->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTRef->SetTitle("reference run");
	  hTRef->Draw("colz");
	  c->Update();
	  c->Write("cTRMS");
	  c->Print(outPDFName.c_str());
	  delete hRMS; delete hTRMS;


	  hRef=(TH2D*)refFile->Get("hMax");
	  hTRef=(TH2D*)refFile->Get("hTMax");
	  hRef->Scale(MAX_EVENT/ref_event); //scale the reference plots with the same total number of events
	  hTRef->SetAxisRange(0,EventRange*EventBin, "X");   //cut the X range so that reference plot and test run plot match.
	  z_max = hRef->GetMaximum();
	  z_min = hRef->GetMinimum();
	  z_maxT = hTRef->GetMaximum();
	  z_minT = hTRef->GetMinimum();
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hMax->GetMaximum() > z_max) z_max=hMax->GetMaximum();
	  if(hMax->GetMinimum() < z_min) z_min=hMax->GetMinimum();
	  hMax->SetAxisRange(z_min-1, z_max+1,"Z");
	  hMax->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hRef->SetAxisRange(z_min-1, z_max+1, "Z");
	  hRef->SetTitle("reference run");
	  hRef->Draw("colz");
	  c->Update();
	  c->Write("cMax");
	  c->Print(outPDFName.c_str());
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hTMax->GetMaximum() > z_maxT) z_maxT=hTMax->GetMaximum();
	  if(hTMax->GetMinimum() < z_minT) z_minT=hTMax->GetMinimum();
	  hTMax->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTMax->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hTRef->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTRef->SetTitle("reference run");
	  hTRef->Draw("colz");
	  c->Update();
	  c->Write("cTMax");
	  c->Print(outPDFName.c_str());
	  delete hMax; delete hTMax;



	  hRef=(TH2D*)refFile->Get("hMin");
	  hTRef=(TH2D*)refFile->Get("hTMin");
	  hRef->Scale(MAX_EVENT/ref_event); //scale the reference plots with the same total number of events
	  hTRef->SetAxisRange(0,EventRange*EventBin, "X");   //cut the X range so that reference plot and test run plot match.
	  z_max = hRef->GetMaximum();
	  z_min = hRef->GetMinimum();
	  z_maxT = hTRef->GetMaximum();
	  z_minT = hTRef->GetMinimum();
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hMin->GetMaximum() > z_max) z_max=hMin->GetMaximum();
	  if(hMin->GetMinimum() < z_min) z_min=hMin->GetMinimum();
	  hMin->SetAxisRange(z_min-1, z_max+1,"Z");
	  hMin->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hRef->SetAxisRange(z_min-1, z_max+1, "Z");
	  hRef->SetTitle("reference run");
	  hRef->Draw("colz");
	  c->Update();
	  c->Write("cMin");
	  c->Print(outPDFName.c_str());
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hTMin->GetMaximum() > z_maxT) z_maxT=hTMin->GetMaximum();
	  if(hTMin->GetMinimum() < z_minT) z_minT=hTMin->GetMinimum();
	  hTMin->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTMin->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hTRef->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTRef->SetTitle("reference run");
	  hTRef->Draw("colz");
	  c->Update();
	  c->Write("cTMin");
	  c->Print(outPDFName.c_str());
	  delete hMin; delete hTMin;


	  hRef=(TH2D*)refFile->Get("hMean");
	  hTRef=(TH2D*)refFile->Get("hTMean");
	  hRef->Scale(MAX_EVENT/ref_event); //scale the reference plots with the same total number of events
	  hTRef->SetAxisRange(0,EventRange*EventBin, "X");   //cut the X range so that reference plot and test run plot match.
	  z_max = hRef->GetMaximum();
	  z_min = hRef->GetMinimum();
	  z_maxT = hTRef->GetMaximum();
	  z_minT = hTRef->GetMinimum();
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hMean->GetMaximum() > z_max) z_max=hMean->GetMaximum();
	  if(hMean->GetMinimum() < z_min) z_min=hMean->GetMinimum();
	  hMean->SetAxisRange(z_min-1, z_max+1,"Z");
	  hMean->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hRef->SetAxisRange(z_min-1, z_max+1, "Z");
	  hRef->SetTitle("reference run");
	  hRef->Draw("colz");
	  c->Update();
	  c->Write("cMean");
	  c->Print(outPDFName.c_str());
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hTMean->GetMaximum() > z_maxT) z_maxT=hTMean->GetMaximum();
	  if(hTMean->GetMinimum() < z_minT) z_minT=hTMean->GetMinimum();
	  hTMean->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTMean->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hTRef->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTRef->SetTitle("reference run");
	  hTRef->Draw("colz");
	  c->Update();
	  c->Write("cTMean");
	  c->Print(outPDFName.c_str());
	  delete hMean; delete hTMean;


	  hRef=(TH2D*)refFile->Get("hNumSample");
	  hTRef=(TH2D*)refFile->Get("hTNumSample");
	  hRef->Scale(MAX_EVENT/ref_event); //scale the reference plots with the same total number of events
	  hTRef->SetAxisRange(0,EventRange*EventBin, "X");   //cut the X range so that reference plot and test run plot match.
	  z_max = hRef->GetMaximum();
	  z_min = hRef->GetMinimum();
	  z_maxT = hTRef->GetMaximum();
	  z_minT = hTRef->GetMinimum();
	  c->Clear();  
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hNumSample->GetMaximum() > z_max) z_max=hNumSample->GetMaximum();
	  if(hNumSample->GetMinimum() < z_min) z_min=hNumSample->GetMinimum();
	  hNumSample->SetAxisRange(z_min-1, z_max+1,"Z");
	  hNumSample->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hRef->SetAxisRange(z_min-1, z_max+1, "Z");
	  hRef->SetTitle("reference run");
	  hRef->Draw("colz");
	  c->Update();
	  c->Write("cNumSample");
	  c->Print(outPDFName.c_str());
	  c->Clear();
	  c->Draw();
	  c->Divide(1,2);
	  c->cd(1);
	  gPad->SetLogz(1);
	  if(hTNumSample->GetMaximum() > z_maxT) z_maxT=hTNumSample->GetMaximum();
	  if(hTNumSample->GetMinimum() < z_minT) z_minT=hTNumSample->GetMinimum();
	  hTNumSample->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTNumSample->Draw("colz");
	  c->cd(2);
	  gPad->SetLogz(1);
	  hTRef->SetAxisRange(z_minT-1, z_maxT+1, "Z");
	  hTRef->SetTitle("reference run");
	  hTRef->Draw("colz");
	  c->Update();
	  c->Write("cTNumSample");
	  c->Print(outPDFName.c_str());
	  delete hNumSample; delete hTNumSample;

	  c->Print((outPDFName+"]").c_str());
	  


	  delete hRef;
	  delete hTRef;
	  delete vRef;
	  delete ref_Huff;
 	  delete c;
	  rootFile.Close();
  }
  


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
