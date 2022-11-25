/***********************************************
* SIDDHARTA-2 Experiment
* Author        : Aleksander Khreptak
* Author e-mail : alex.nuclearboy@gmail.com
* Created on Nov 25, 2022
* Licensed under the Apache License, Version 2.0
***********************************************/

//// Peak detection routine

#include <TROOT.h>
#include <TFile.h>
#include <Riostream.h>
#include <stdlib.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TMinuit.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TSystemDirectory.h>
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include <TSystemFile.h>
#include <TMath.h>
#include <TF1.h>
#include "peakfinder.h"

using namespace std;

void peakfinder() {

  // Get the root input file:
  static const Int_t MAXNFILES = 10;
  TFile *filein[MAXNFILES];
  TFile *fileout[MAXNFILES];

  Int_t nfile = 0;

  TString fname[MAXNFILES];

  TString path = "/home/alex/SIDDHARTA2/SIDDHARTA2_calibration/rootdata/SIDDHARTA2_xrays/input/20220602/";

  TSystemDirectory dir(path, path);
  TList *files = dir.GetListOfFiles();
  if (files) {
    TSystemFile *file; TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname[nfile] = file->GetName();
      if (!file->IsDirectory() && fname[nfile].EndsWith(".root")) {
        filein[nfile] = new TFile(path+fname[nfile].Data(),"read");
        cout<<"Input file #"<<nfile+1<<": "<<fname[nfile].Data()<<endl;
        // Define the output file:
        fileout[nfile] = new TFile("files/histos_" + fname[nfile],"RECREATE");
        nfile++;
      }
    }
  }

  cout<<"Number of root files in the path: "<<nfile<<endl;

  for(Int_t i=0; i<nfile; i++) {  //start the loop on the root input file

    cout<<"Reading input file #"<<i+1<<endl;

    Int_t column, row;
    Double_t runtime=0.;

    Bool_t ISGOODtime;

    nbinsadc = nbinsadc/rebinFactor;

    for(Int_t iBUS=1;iBUS<nBUS;iBUS++) {
      for(Int_t iSDD=1;iSDD<nSDD;iSDD++) {
        hADC[iBUS][iSDD] = new TH1F(Form("hADC_bus%i_sdd%i",iBUS,iSDD),Form("hADC_bus%i_sdd%i",iBUS,iSDD),nbinsadc,minadc,maxadc);
        hADC[iBUS][iSDD]->GetXaxis()->SetTitle("ADC [channel]");
        hADC[iBUS][iSDD]->GetYaxis()->SetTitle("counts / channel");
        hADCall[iBUS][iSDD] = new TH1F(Form("hADCall_bus%i_sdd%i",iBUS,iSDD),Form("hADCall_bus%i_sdd%i",iBUS,iSDD),nbinsadc,minadc,maxadc);
        hADCall[iBUS][iSDD]->GetXaxis()->SetTitle("ADC [channel]");
        hADCall[iBUS][iSDD]->GetYaxis()->SetTitle("counts / channel");
        hCrosstalk[iBUS][iSDD] = new TH1F(Form("hCrosstalk_bus%i_sdd%i",iBUS,iSDD),Form("hCrosstalk_bus%i_sdd%i",iBUS,iSDD),nbinsadc,minadc,maxadc);
        hCrosstalk[iBUS][iSDD]->GetXaxis()->SetTitle("ADC [channel]");
        hCrosstalk[iBUS][iSDD]->GetYaxis()->SetTitle("counts / channel");
        empty[iBUS][iSDD]=true;
      }
    }

    TH2F *hSDDmap = new TH2F ("hSDDmap","hSDDmap",48,1,49,8,1,9);                       // map of the SDDs position around the target
    TH2D *hSDDrate = new TH2D ("hSDDrate","hSDDrate",48,1,49,8,1,9);                    // map of the SDDs rate around the target
    TH2D *hSDDrate_sig = new TH2D ("hSDDrate_sig","hSDDrate_sig",48,1,49,8,1,9);        // map of the SDDs signal rate around the target
    TH2D *hSDDrate_noise = new TH2D ("hSDDrate_noise","hSDD_rate_noise",48,1,49,8,1,9); // map of the SDDs noise rate around the target

    filein[i]->cd();
    tree = (TTree*)filein[i]->Get("ft10");
    nentries = tree->GetEntries();

    tree->GetEntry(0);
    Int_t timestart = (Int_t)tree->GetBranch("date")->GetLeaf("date")->GetValue(0);
    tree->GetEntry(nentries-1);
    Int_t timeend = (Int_t)tree->GetBranch("date")->GetLeaf("date")->GetValue(0);
    Int_t DeltaTime = timeend-timestart;
    runtime = runtime+DeltaTime;
    //cout<<"timestart = "<<timestart<<endl;
    //cout<<"timeend = "<<timeend<<endl;
    cout<<"Runtime: "<<runtime<<endl;

    // Set branch addresses and branch pointers
    tree->SetBranchAddress("buf",&buf);
    tree->SetBranchAddress("kt",&kt);
    tree->SetBranchAddress("nhits",&nhits);
    tree->SetBranchAddress("bus",&bus);
    tree->SetBranchAddress("evnr",&evnr);
    tree->SetBranchAddress("ht",&ht);
    tree->SetBranchAddress("trigg",&trigg);
    tree->SetBranchAddress("sdd",&sdd);
    tree->SetBranchAddress("adc",&adc);
    tree->SetBranchAddress("drift",&drift);
    tree->SetBranchAddress("date",&date);
    tree->SetBranchAddress("ie",&ie);
    tree->SetBranchAddress("ip",&ip);
    tree->SetBranchAddress("dum",&dum);
    tree->SetBranchAddress("edges_veto",&edges_veto);
    tree->SetBranchAddress("v_edge",&v_edge);
    tree->SetBranchAddress("v_val",&v_val);
    tree->SetBranchAddress("v_ch",&v_ch);
    tree->SetBranchAddress("v_found",&v_found);
    tree->SetBranchAddress("v_ch_r",&v_ch_r);
    tree->SetBranchAddress("v_tdc",&v_tdc);
    tree->SetBranchAddress("v_tot",&v_tot);

    for(Int_t ient=0; ient<nentries; ient++) {
      if(ient%100==0) {cout<<"Reading file "<< i+1 << " entry "<< ient <<"/"<< nentries << endl;}
      tree->GetEntry(ient);

      for(Int_t ihit=0; ihit<nhits; ihit++) {
        ISGOODtime=true;
        if(ihit>0){ISGOODtime = CrossTalkTiming(drift[ihit],drift[ihit-1]);}
        SDDhitmap(sdd[ihit], bus[ihit], column, row);
        hSDDmap->SetBinContent(column,row,sdd[ihit]);
        hSDDrate->Fill(column,row);
        if(adc[ihit]>1200){hSDDrate_sig->Fill(column,row);}
        if(adc[ihit]<=1200){hSDDrate_noise->Fill(column,row);}
        hADCall[bus[ihit]][sdd[ihit]]->Fill(adc[ihit]);
        if(ISGOODtime){hADC[bus[ihit]][sdd[ihit]]->Fill(adc[ihit]); empty[bus[ihit]][sdd[ihit]]=false;}
        if(!ISGOODtime){hCrosstalk[bus[ihit]][sdd[ihit]]->Fill(adc[ihit]);}
      }

    }

    hSDDrate->Scale(1/runtime);
    hSDDrate_sig->Scale(1/runtime);
    hSDDrate_noise->Scale(1/runtime);

    for(int i=1; i<=49; i++) {
      for(int j=1;j<9;j++) {
        if(hSDDrate->GetBinContent(i,j)<0.001){
          hSDDrate->SetBinContent(i,j,0);
          hSDDrate_sig->SetBinContent(i,j,0);
          hSDDrate_noise->SetBinContent(i,j,0);
        }
      }
    }

    // Write the histograms
    fileout[i]->cd();
    for(Int_t iBUS=1; iBUS<nBUS; iBUS++) {
      for(Int_t iSDD=1; iSDD<nSDD; iSDD++) {
        if(hADC[iBUS][iSDD]->GetEntries()>0) {
          hADCall[iBUS][iSDD]->Write();
          hADC[iBUS][iSDD]->Write();
          hCrosstalk[iBUS][iSDD]->Write();
        }
      }
    }
    hSDDmap->Write();
    hSDDrate->Write();
    hSDDrate_sig->Write();
    hSDDrate_noise->Write();

    fileout[i]->Close();

    cout<<"Output file #"<<i+1<<": "<<fileout[i]->GetName()<<endl;

    // Creating histgrams
    // Settings global style of histograms
    gStyle->SetOptStat(kFALSE);  //kFALSE or kTRUE
    //gStyle->SetPalette(1,0);
    gStyle->SetPadGridX(kTRUE);
    gStyle->SetPadGridY(kTRUE);
    // Canvas setting
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadRightMargin(0.1);
    gStyle->SetPadBottomMargin(0.12);
    //gStyle->SetPadTopMargin(0.1);
    //gStyle->SetHistLineStyle(1.);
    //gStyle->SetHistLineColor(1);
    //gStyle->SetHistLineWidth(2);

    gStyle->SetTitleSize(0.05,"X");
    gStyle->SetTitleSize(0.05,"Y");
    gStyle->SetTitleOffset(0.9,"X");
    gStyle->SetTitleOffset(1.1,"Y");
    gStyle->SetLabelSize(0.045,"X");
    gStyle->SetLabelSize(0.045,"Y");

    // Text style setting
    gStyle->SetTitleFont(42,"XYZ");
    gStyle->SetLabelFont(42,"XYZ");
    gStyle->SetTextFont(42);

    fname[i].Remove(23,fname[i].Length());

    //////
    TPad *pad[nBUS][nSDD];
    for(Int_t j=1;j<nBUS;j++) {
      for(Int_t k=0;k<8;k++) {
        myCanvas[j][k] = new TCanvas(fname[i]+Form("_hADC_bus%d_sdd%d_%d",j,k*8+1,8*(1+k)),Form("hADC_bus%d_sdd%d_%d",j,k*8+1,8*(k+1)),0,0,1240,1748);
        myCanvas[j][k]->Divide(2,4);  //column, line
        for(Int_t l=1;l<9;l++) {
          pad[j][k] = (TPad*)myCanvas[j][k]->cd(l);
          pad[j][k]->SetLogy();
          pad[j][k]->SetGrid();
          if(!(hADCall[j][l+8*k]->GetEntries()>0)) continue;
          hADCall[j][l+8*k]->SetTitle(Form("Fluorescence spectrum (BUS: %d, SDD: %d)",j,l+8*k));
          hADCall[j][l+8*k]->SetAxisRange(1500,4500,"X");
          hADCall[j][l+8*k]->UseCurrentStyle();
          hADCall[j][l+8*k]->SetLineWidth(2);       
          hADCall[j][l+8*k]->DrawCopy();
          hADC[j][l+8*k]->SetLineColor(2);
          hADC[j][l+8*k]->DrawCopy("same");
          hCrosstalk[j][l+8*k]->SetLineColor(1);
          hCrosstalk[j][l+8*k]->DrawCopy("same");
        }
        myCanvas[j][k]->Print("plots/hADC_spectra_"+fname[i]+Form("_bus%d_sdd%d_%d.pdf",j,k*8+1,8*(k+1)), "pdf");
      }
    }

  }

}

void SDDhitmap(Int_t sddnumber, Int_t busnumber, Int_t &column, Int_t &row) {

  row=0; column=0;
  Int_t SFERA=0;

  //Back side view
  if(sddnumber<=16){sddnumber=sddnumber;SFERA=0;}
  if(sddnumber>=17 && sddnumber<=32){sddnumber=sddnumber-16;SFERA=2;}
  if(sddnumber>=33 && sddnumber<=48){sddnumber=sddnumber-(16*2);SFERA=4;}
  if(sddnumber>=49 && sddnumber<=64){sddnumber=sddnumber-(16*3);SFERA=6;}

  if(sddnumber==1){row=1; column=1+SFERA+((busnumber-1)*8);}     if(sddnumber== 9){row=5; column=1+SFERA+((busnumber-1)*8);}
  if(sddnumber==2){row=4; column=1+SFERA+((busnumber-1)*8);}     if(sddnumber==10){row=8; column=1+SFERA+((busnumber-1)*8);}
  if(sddnumber==3){row=2; column=1+SFERA+((busnumber-1)*8);}     if(sddnumber==11){row=7; column=1+SFERA+((busnumber-1)*8);}
  if(sddnumber==4){row=3; column=1+SFERA+((busnumber-1)*8);}     if(sddnumber==12){row=6; column=1+SFERA+((busnumber-1)*8);}
  if(sddnumber==5){row=3; column=2+SFERA+((busnumber-1)*8);}     if(sddnumber==13){row=6; column=2+SFERA+((busnumber-1)*8);}
  if(sddnumber==6){row=2; column=2+SFERA+((busnumber-1)*8);}     if(sddnumber==14){row=7; column=2+SFERA+((busnumber-1)*8);}
  if(sddnumber==7){row=1; column=2+SFERA+((busnumber-1)*8);}     if(sddnumber==15){row=5; column=2+SFERA+((busnumber-1)*8);}
  if(sddnumber==8){row=4; column=2+SFERA+((busnumber-1)*8);}     if(sddnumber==16){row=8; column=2+SFERA+((busnumber-1)*8);}

}

bool CrossTalkTiming(Short_t drift, Short_t drift_pre) {

  Bool_t ISGOODtime = true;
  Double_t t1=0, t2=0, timediff=0;
  t1=drift+32768;
  t2=drift_pre+32768;
  if(t1>t2) {timediff = t1-t2;}
  if(t2>t1) {timediff = t1+(32768.*2)-t2;}
  if((timediff>0. && timediff<625)) {ISGOODtime = false;}   //625 = 5 microseconds
  return ISGOODtime;

}

Int_t SFERAnumber(Int_t sdd) {

  Int_t SFERA=0;
  if(sdd<=16){SFERA=1;}
  if(sdd>=17 && sdd<=32){SFERA=2;}
  if(sdd>=33 && sdd<=48){SFERA=3;}
  if(sdd>=49 && sdd<=64){SFERA=4;}

  return SFERA;

}
