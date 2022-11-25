#ifndef ___PEAKFINDER_H___
#define ___PEAKFINDER_H___

#include  <TGraph.h>
#include  <TGraphErrors.h>
#include  <TClonesArray.h>
#include  <TF1.h>
#include  <TH1.h>
#include  <TH2.h>
#include  <TH1F.h>
#include  <TH2F.h>
#include  <TH1D.h>
#include  <TH2D.h>
#include  <TCanvas.h>
#include  <TLegend.h>
#include  <TPaveText.h>
#include  <TStyle.h>
#include  <TString.h>

static const Int_t nBUS = 7;   //last BUS #6
static const Int_t nSDD = 65;  //last SDD #64

// Histograms booking:
Int_t nbinsadc = 10000;
Int_t minadc = 0.;
Int_t maxadc = 10000.;
Int_t rebinFactor = 4;

TH1F *hADC[nBUS][nSDD],*hADCall[nBUS][nSDD],*hADCPF[nBUS][nSDD];
TH1F *hCrosstalk[nBUS][nSDD];

Bool_t empty[nBUS][nSDD];

TTree *tree;

// Fixed dimensions of arrays stored in the TTree, if any
// Declaration of leaf types:
static const Int_t MAXHITS = 40000;
Int_t   buf   = 0;
Short_t kt[4] = {0};
Int_t   nhits = 0 ;
Int_t   bus[MAXHITS]   = {0};
Short_t evnr[MAXHITS]  = {0};
Short_t ht[MAXHITS]    = {0};
Short_t trigg[MAXHITS] = {0};
Short_t sdd[MAXHITS]   = {0};
Short_t adc[MAXHITS]   = {0};
Short_t drift[MAXHITS] = {0};
Int_t   date = 0;
Short_t ie   = 0;
Short_t ip   = 0;
Short_t dum  = 0;

Int_t   edges_veto  = 0;
Short_t v_edge[280] = {0};
Short_t v_val[280]  = {0};
Short_t v_ch[280]   = {0};
Int_t   v_found     = 0;
Short_t v_ch_r[120] = {0};
Short_t v_tdc[120]  = {0};
Short_t v_tot[120]  = {0};

Long64_t nentries = 0;

void    SDDhitmap(Int_t sddnumber, Int_t busnumber, Int_t &column, Int_t &row);
Bool_t  CrossTalkTiming(Short_t drift, Short_t drift_pre);
Int_t   SFERAnumber(Int_t sdd);

TCanvas* myCanvas[nBUS][nSDD];

#endif
