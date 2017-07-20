#include <TAxis.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TVector2.h>

// Plotting macro for trigger efficiency plots
// Use triggerEff.C and the batch submission scripts to get the input needed for this macro

// Hardcoded bin settings
  double pi=TMath::Pi();
  //const int nBinsMetOnly=59; double metOnlyBins[nBinsMetOnly+1]={0,5,10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 500, 600, 700, 800, 1000, 1200, 1500};
  const int nBinsMetOnly=44; double metOnlyBins[nBinsMetOnly+1]={0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 220, 240, 260, 280, 300, 350, 400, 450, 500, 600, 700, 800, 1000, 1200};
  //const int nBinsMetOnly=36; double metOnlyBins[nBinsMetOnly+1]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 250, 300, 350, 400, 600, 1000, 1500};
  //const int nBinsMet=7; double metBins[nBinsMet+1]={100,150,200,250,300,400,800,1500};
  const int nBinsMet=7; double metBins[nBinsMet+1]={100,125,150,175,200,250,300,400};
  const int nBinsPhi=6; double phiBins[nBinsPhi+1]={0.,pi/3.,2.*pi/3.,pi,4.*pi/3., 5.*pi/3., 2.*pi};
  const int nBinsLeptonPt=6; double leptonPtBins[nBinsLeptonPt+1]={45,50,60,80,100,150,300};
  const int nBinsJetPt=8; double jetPtBins[nBinsJetPt+1]={80,100,120,140,160,200,250,500,1000};

  
  // fine eta binning that is used for the JECs
  //const int nBinsJetEta=36; double jetEtaBins[nBinsJetEta+1]={-5.191, -3.839, -3.489, -3.139, -2.964, -2.853, -2.650, -2.500, -2.322, -2.172, -1.930, -1.653, -1.479, -1.305, -1.044, -0.783, -0.522, -0.261,  0.000,  0.261,  0.522,  0.783,  1.044,  1.305,  1.479,  1.653,  1.930,  2.172,  2.322,  2.500,  2.650,  2.853,  2.964,  3.139,  3.489,  3.839,  5.191};
  
  // coarser eta binning for low statistics
  const int nBinsJetEta=9; double jetEtaBins[nBinsJetEta+1]={-5.191, -2.964, -2.5, -1.479, -0.522,  0.522,  1.479,  2.5,  2.964,  5.191};
  
  //const int nBinsMjj=4; double mjjBins[nBinsMjj+1]={0,500,1000,2000,4000};
  const int nBinsMjj=4; double mjjBins[nBinsMjj+1]={0,200,500,1000,2000};

void plotEfficiency(
  TString inputFilePath, //full path to input tree
  TString flavor="electrons", // (electrons|muons)
  TString mode="monojet", // (monojet|vbf)
  bool type1=false,
  TString extra_string="nice",
  int fitFunction=1 // 1=tanh, 2=arctan, 3=erf, 4=generalized logistic
) {
  TString fitFunctionName;
  if     (fitFunction==1) fitFunctionName="tanh"; 
  else if(fitFunction==2) fitFunctionName="arctan"; 
  else if(fitFunction==3) fitFunctionName="erf"; 
  else if(fitFunction==4) fitFunctionName="glf"; 
  TFile *inputFile = TFile::Open(inputFilePath, "READ");
  // set up input tree
  TTree *effTree = (TTree*) inputFile->Get("effTree");
  bool passMetTriggers, passJetMetCleaning, jet1ExtraClean, jet2ExtraClean, jet3ExtraClean, jet4ExtraClean;
  float leptonPt, leptonEta, leptonPhi, pfMetNoMu, pfMetNoMuPhi, jet1Pt, jet1Eta, jet1Phi, jet1M, jet1PtUp, jet1PtDown, jet2Pt, jet2Eta, jet2Phi, jet2M, jet2PtUp, jet2PtDown,
        jet3Pt, jet3Eta, jet3Phi, jet3M, jet3PtUp, jet3PtDown, jet4Pt, jet4Eta, jet4Phi, jet4M, jet4PtUp, jet4PtDown;
  unsigned nTau;
  effTree->SetBranchAddress("passMetTriggers"    , &passMetTriggers   );
  effTree->SetBranchAddress("passJetMetCleaning" , &passJetMetCleaning);
  effTree->SetBranchAddress("leptonPt"           , &leptonPt          );
  effTree->SetBranchAddress("leptonEta"          , &leptonEta         );
  effTree->SetBranchAddress("leptonPhi"          , &leptonPhi         );
  effTree->SetBranchAddress("pfMetNoMu"          , &pfMetNoMu         );
  effTree->SetBranchAddress("pfMetNoMuPhi"       , &pfMetNoMuPhi      );
  effTree->SetBranchAddress("jet1Pt"             , &jet1Pt            );
  effTree->SetBranchAddress("jet1Eta"            , &jet1Eta           );
  effTree->SetBranchAddress("jet1Phi"            , &jet1Phi           );
  effTree->SetBranchAddress("jet1M"              , &jet1M             );
  //effTree->SetBranchAddress("jet1PtUp"           , &jet1PtUp          );
  //effTree->SetBranchAddress("jet1PtDown"         , &jet1PtDown        );
  effTree->SetBranchAddress("jet1ExtraClean"     , &jet1ExtraClean    );
  effTree->SetBranchAddress("jet2Pt"             , &jet2Pt            );
  effTree->SetBranchAddress("jet2Eta"            , &jet2Eta           );
  effTree->SetBranchAddress("jet2Phi"            , &jet2Phi           );
  effTree->SetBranchAddress("jet2M"              , &jet2M             );
  //effTree->SetBranchAddress("jet2PtUp"           , &jet2PtUp          );
  //effTree->SetBranchAddress("jet2PtDown"         , &jet2PtDown        );
  effTree->SetBranchAddress("jet2ExtraClean"     , &jet2ExtraClean    );
  effTree->SetBranchAddress("jet3Pt"             , &jet3Pt            );
  effTree->SetBranchAddress("jet3Eta"            , &jet3Eta           );
  effTree->SetBranchAddress("jet3Phi"            , &jet3Phi           );
  effTree->SetBranchAddress("jet3M"              , &jet3M             );
  //effTree->SetBranchAddress("jet3PtUp"           , &jet3PtUp          );
  //effTree->SetBranchAddress("jet3PtDown"         , &jet3PtDown        );
  effTree->SetBranchAddress("jet3ExtraClean"     , &jet3ExtraClean    );
  effTree->SetBranchAddress("jet4Pt"             , &jet4Pt            );
  effTree->SetBranchAddress("jet4Eta"            , &jet4Eta           );
  effTree->SetBranchAddress("jet4Phi"            , &jet4Phi           );
  effTree->SetBranchAddress("jet4M"              , &jet4M             );
  //effTree->SetBranchAddress("jet4PtUp"           , &jet4PtUp          );
  //effTree->SetBranchAddress("jet4PtDown"         , &jet4PtDown        );
  effTree->SetBranchAddress("jet4ExtraClean"     , &jet4ExtraClean    );
  effTree->SetBranchAddress("nTau"               , &nTau              );

  // set up histograms
  TH2F *hPass_metPhi_met = new TH2F("hPass_metPhi_met", "hPass_metPhi_met", nBinsMet, metBins, nBinsPhi, phiBins); hPass_metPhi_met->SetDirectory(0); hPass_metPhi_met->Sumw2();
  TH2F *hPass_jet1Pt_met = new TH2F("hPass_jet1Pt_met", "hPass_jet1Pt_met", nBinsMet, metBins, nBinsJetPt, jetPtBins); hPass_jet1Pt_met->SetDirectory(0); hPass_jet1Pt_met->Sumw2();
  TH2F *hPass_jet1Eta_met = new TH2F("hPass_jet1Eta_met", "hPass_jet1Eta_met", nBinsMet, metBins, nBinsJetEta, jetEtaBins); hPass_jet1Eta_met->SetDirectory(0); hPass_jet1Eta_met->Sumw2();
  TH2F *hPass_jet2Eta_met = new TH2F("hPass_jet2Eta_met", "hPass_jet2Eta_met", nBinsMet, metBins, nBinsJetEta, jetEtaBins); hPass_jet2Eta_met->SetDirectory(0); hPass_jet2Eta_met->Sumw2();
  TH2F *hPass_mjj_met = new TH2F("hPass_mjj_met", "hPass_mjj_met", nBinsMet, metBins, nBinsMjj, mjjBins); hPass_mjj_met->SetDirectory(0); hPass_mjj_met->Sumw2();
  TH2F *hPass_leptonPt_met = new TH2F("hPass_leptonPt_met", "hPass_leptonPt_met", nBinsMet, metBins, nBinsLeptonPt, leptonPtBins); hPass_mjj_met->SetDirectory(0); hPass_leptonPt_met->Sumw2();
  TH2F *hFail_metPhi_met = new TH2F("hFail_metPhi_met", "hFail_metPhi_met", nBinsMet, metBins, nBinsPhi, phiBins); hFail_metPhi_met->SetDirectory(0); hFail_metPhi_met->Sumw2();
  TH2F *hFail_jet1Pt_met = new TH2F("hFail_jet1Pt_met", "hFail_jet1Pt_met", nBinsMet, metBins, nBinsJetPt, jetPtBins); hFail_jet1Pt_met->SetDirectory(0); hFail_jet1Pt_met->Sumw2();
  TH2F *hFail_jet1Eta_met = new TH2F("hFail_jet1Eta_met", "hFail_jet1Eta_met", nBinsMet, metBins, nBinsJetEta, jetEtaBins); hFail_jet1Eta_met->SetDirectory(0); hFail_jet1Eta_met->Sumw2();
  TH2F *hFail_jet2Eta_met = new TH2F("hFail_jet2Eta_met", "hFail_jet2Eta_met", nBinsMet, metBins, nBinsJetEta, jetEtaBins); hFail_jet2Eta_met->SetDirectory(0); hFail_jet2Eta_met->Sumw2();
  TH2F *hFail_mjj_met = new TH2F("hFail_mjj_met", "hFail_mjj_met", nBinsMet, metBins, nBinsMjj, mjjBins); hFail_mjj_met->SetDirectory(0); hFail_mjj_met->Sumw2();
  TH2F *hFail_leptonPt_met = new TH2F("hFail_leptonPt_met", "hFail_leptonPt_met", nBinsMet, metBins, nBinsLeptonPt, leptonPtBins); hFail_mjj_met->SetDirectory(0); hFail_leptonPt_met->Sumw2();
  
  TH1F *hPass_met = new TH1F("hPass_met", "hPass_met", nBinsMetOnly, metOnlyBins);
  TH1F *hFail_met = new TH1F("hFail_met", "hFail_met", nBinsMetOnly, metOnlyBins);
  Long64_t nEntries=effTree->GetEntries();
  for(Long64_t iEntry=0; iEntry<nEntries; iEntry++) {
    effTree->GetEntry(iEntry);

    double mjj=-1;
    TVector2 pfMetP2, pfMetNoMuP2, leptonP2;
    TLorentzVector selectedJet1, selectedJet2, selectedDijetSystem;
    if(mode=="monojet") {
      if     (jet1Pt>=100 && fabs(jet1Eta)<2.5 && jet1ExtraClean) selectedJet1.SetPtEtaPhiM(jet1Pt, jet1Eta, jet1Phi, jet1M); 
      else if(jet2Pt>=100 && fabs(jet2Eta)<2.5 && jet2ExtraClean) selectedJet1.SetPtEtaPhiM(jet2Pt, jet2Eta, jet2Phi, jet2M); 
      else if(jet3Pt>=100 && fabs(jet3Eta)<2.5 && jet3ExtraClean) selectedJet1.SetPtEtaPhiM(jet3Pt, jet3Eta, jet3Phi, jet3M); 
      else if(jet4Pt>=100 && fabs(jet4Eta)<2.5 && jet4ExtraClean) selectedJet1.SetPtEtaPhiM(jet4Pt, jet4Eta, jet4Phi, jet4M); 
      else continue;
      // if(!passJetMetCleaning) continue; 
      if(flavor=="electrons") {
        if(leptonPt<45) continue;
        //if(fabs(leptonEta) > 1.4442) continue;
      } else if (flavor=="muons") {
        if(leptonPt<25) continue;
      }
    } else if(mode=="vbf") {
      unsigned char nSelectedJet1, nSelectedJet2;
      // Find a jet with pT at least 80 GeV that is either forward or passes the cleaning cuts
      if     (jet1Pt>=80 && (fabs(jet1Eta)>=2.5 || jet1ExtraClean)) { selectedJet1.SetPtEtaPhiM(jet1Pt, jet1Eta, jet1Phi, jet1M); nSelectedJet1=1;}
      //else if(jet2Pt>=80 && (fabs(jet2Eta)>=2.5 || jet2ExtraClean)) { selectedJet1.SetPtEtaPhiM(jet2Pt, jet2Eta, jet2Phi, jet2M); nSelectedJet1=2;}
      //else if(jet3Pt>=80 && (fabs(jet3Eta)>=2.5 || jet3ExtraClean)) { selectedJet1.SetPtEtaPhiM(jet3Pt, jet3Eta, jet3Phi, jet3M); nSelectedJet1=3;}
      //else if(jet4Pt>=80 && (fabs(jet4Eta)>=2.5 || jet4ExtraClean)) { selectedJet1.SetPtEtaPhiM(jet4Pt, jet4Eta, jet4Phi, jet4M); nSelectedJet1=4;}
      else continue;
      // Find another jet with pT at least 40 GeV
      if     (nSelectedJet1!=1 && jet1Pt>=40) { selectedJet2.SetPtEtaPhiM(jet1Pt, jet1Eta, jet1Phi, jet1M); nSelectedJet2=1; }
      else if(nSelectedJet1!=2 && jet2Pt>=40) { selectedJet2.SetPtEtaPhiM(jet2Pt, jet2Eta, jet2Phi, jet2M); nSelectedJet2=2; }
      else if(nSelectedJet1!=3 && jet3Pt>=40) { selectedJet2.SetPtEtaPhiM(jet3Pt, jet3Eta, jet3Phi, jet3M); nSelectedJet2=3; }
      else if(nSelectedJet1!=4 && jet4Pt>=40) { selectedJet2.SetPtEtaPhiM(jet4Pt, jet4Eta, jet4Phi, jet4M); nSelectedJet2=4; }
      else continue;
      
      selectedDijetSystem = selectedJet1 + selectedJet2;
      mjj = selectedDijetSystem.M();

      // if(!passJetMetCleaning) continue;
      if(flavor=="electrons") {
        if(leptonPt<45) continue;
      } else if (flavor=="muons") {
        if(leptonPt<25) continue;
      }
    }  
    //leptonP2.SetMagPhi(leptonPt, leptonPhi);
    //pfMetNoMuP2.SetMagPhi(pfMetNoMu, pfMetNoMuPhi);
    //pfMetP2 = (flavor=="electrons") ? pfMetNoMuP2 : pfMetNoMuP2 - leptonP2;
    //bool passMTwindow=false;
    //double mT = TMath::Sqrt( 2.*leptonPt*pfMetP2.Mod() * (1. - TMath::Cos(pfMetP2.DeltaPhi(leptonP2) )));
    //if( TMath::Abs(mT-80.)>10.) continue;
    //if( mT < 80.) continue;

    if(passMetTriggers) {
      hPass_metPhi_met->Fill(TMath::Min(double(pfMetNoMu), metBins[nBinsMet]-0.001), pfMetNoMuPhi);
      hPass_jet1Pt_met->Fill(TMath::Min(double(pfMetNoMu), metBins[nBinsMet]-0.001), TMath::Min(double(selectedJet1.Pt()), jetPtBins[nBinsJetPt]-0.001));
      hPass_jet1Eta_met->Fill(TMath::Min(double(pfMetNoMu), metBins[nBinsMet]-0.001), selectedJet1.Eta());
      hPass_jet2Eta_met->Fill(TMath::Min(double(pfMetNoMu), metBins[nBinsMet]-0.001), selectedJet2.Eta());
      hPass_mjj_met->Fill(TMath::Min(double(pfMetNoMu), metBins[nBinsMet]-0.001), TMath::Min(mjj, mjjBins[nBinsMjj]));
      hPass_leptonPt_met->Fill(TMath::Min(double(pfMetNoMu), metBins[nBinsMet]-0.001), TMath::Min(double(leptonPt), leptonPtBins[nBinsLeptonPt]));
      hPass_met->Fill(TMath::Min(double(pfMetNoMu), metOnlyBins[nBinsMetOnly]-0.001));
    } else {
      hFail_metPhi_met->Fill(TMath::Min(double(pfMetNoMu), metBins[nBinsMet]-0.001), pfMetNoMuPhi);
      hFail_jet1Pt_met->Fill(TMath::Min(double(pfMetNoMu), metBins[nBinsMet]-0.001), TMath::Min(double(selectedJet1.Pt()), jetPtBins[nBinsJetPt]-0.001));
      hFail_jet1Eta_met->Fill(TMath::Min(double(pfMetNoMu), metBins[nBinsMet]-0.001), selectedJet1.Eta());
      hFail_jet2Eta_met->Fill(TMath::Min(double(pfMetNoMu), metBins[nBinsMet]-0.001), selectedJet2.Eta());
      hFail_mjj_met->Fill(TMath::Min(double(pfMetNoMu), metBins[nBinsMet]-0.001), TMath::Min(mjj, mjjBins[nBinsMjj]));
      hFail_leptonPt_met->Fill(TMath::Min(double(pfMetNoMu), metBins[nBinsMet]-0.001), TMath::Min(double(leptonPt), leptonPtBins[nBinsLeptonPt]-0.001));
      hFail_met->Fill(TMath::Min(double(pfMetNoMu), metOnlyBins[nBinsMetOnly]-0.001));
      //if(pfMetNoMu>400) printf("event dump: run %d event %d
    }

  }
  // Clone histograms for the total and efficiency from passing histos
  TH2F* hTotal_metPhi_met=(TH2F*)hPass_metPhi_met->Clone("hTotal_metPhi_met"); hTotal_metPhi_met->SetDirectory(0);
  TH2F* hTotal_jet1Pt_met=(TH2F*)hPass_jet1Pt_met->Clone("hTotal_jet1Pt_met"); hTotal_jet1Pt_met->SetDirectory(0);
  TH2F* hTotal_jet1Eta_met=(TH2F*)hPass_jet1Eta_met->Clone("hTotal_jet1Eta_met"); hTotal_jet1Eta_met->SetDirectory(0);
  TH2F* hTotal_jet2Eta_met=(TH2F*)hPass_jet2Eta_met->Clone("hTotal_jet2Eta_met"); hTotal_jet2Eta_met->SetDirectory(0);
  TH2F* hTotal_mjj_met=(TH2F*)hPass_mjj_met->Clone("hTotal_mjj_met"); hTotal_mjj_met->SetDirectory(0);
  TH2F* hTotal_leptonPt_met=(TH2F*)hPass_leptonPt_met->Clone("hTotal_leptonPt_met"); hTotal_leptonPt_met->SetDirectory(0);
  TH2F* hEff_metPhi_met=(TH2F*)hPass_metPhi_met->Clone("hEff_metPhi_met"); hEff_metPhi_met->SetDirectory(0);
  TH2F* hEff_jet1Pt_met=(TH2F*)hPass_jet1Pt_met->Clone("hEff_jet1Pt_met"); hEff_jet1Pt_met->SetDirectory(0);
  TH2F* hEff_jet1Eta_met=(TH2F*)hPass_jet1Eta_met->Clone("hEff_jet1Eta_met"); hEff_jet1Eta_met->SetDirectory(0);
  TH2F* hEff_jet2Eta_met=(TH2F*)hPass_jet2Eta_met->Clone("hEff_jet2Eta_met"); hEff_jet2Eta_met->SetDirectory(0);
  TH2F* hEff_mjj_met=(TH2F*)hPass_mjj_met->Clone("hEff_mjj_met"); hEff_mjj_met->SetDirectory(0);
  TH2F* hEff_leptonPt_met=(TH2F*)hPass_leptonPt_met->Clone("hEff_leptonPt_met"); hEff_leptonPt_met->SetDirectory(0);
  
  TH1F* hTotal_met=(TH1F*)hPass_met->Clone("hTotal_met"); hTotal_met->SetDirectory(0);
  TGraphAsymmErrors *hEff_met=new TGraphAsymmErrors;
  
  // Add failing histo to get the total
  hTotal_metPhi_met->Add(hFail_metPhi_met);
  hTotal_jet1Pt_met->Add(hFail_jet1Pt_met);
  hTotal_jet1Eta_met->Add(hFail_jet1Eta_met);
  hTotal_jet2Eta_met->Add(hFail_jet2Eta_met);
  hTotal_mjj_met->Add(hFail_mjj_met);
  hTotal_leptonPt_met->Add(hFail_leptonPt_met);
  hTotal_met->Add(hFail_met);
  
  // Divide histograms (2D plots don't have asymmetric errors plotted for now) 
  hEff_metPhi_met->Divide(hTotal_metPhi_met);
  hEff_jet1Pt_met->Divide(hTotal_jet1Pt_met);
  hEff_jet1Eta_met->Divide(hTotal_jet1Eta_met);
  hEff_jet2Eta_met->Divide(hTotal_jet2Eta_met);
  hEff_mjj_met->Divide(hTotal_mjj_met);
  hEff_leptonPt_met->Divide(hTotal_leptonPt_met);
  hEff_met->Divide(hPass_met, hTotal_met);
  
  //for(int i=0; i<=hEff_met->GetN(); i++) {
  //  if(i>=hPass_met->GetNbinsX()) break;
  //  double x,y;
  //  hEff_met->GetPoint(i,x,y);
  //  double error = TMath::Sqrt(x);
  //  hEff_met->SetPointEXlow(i,error);
  //  hEff_met->SetPointEXhigh(i,error);
  //}
  
  // BEGIN DRAWING

  // Set some options
  char title[400]; 
  if(flavor=="electrons") {
    if(type1) sprintf(title,"Efficiency of HLT_PFMETTypeOne*_PFMHT*_IDTight in SingleElectron Run2017 A/B");
    else      sprintf(title,"Efficiency of HLT_PFMETNoMu*_PFMHTNoMu*_IDTight in SingleElectron Run2017 A/B");
  } else sprintf(title, "Efficiency of HLT_PFMETNoMu*_PFMHTNoMu*_IDTight in SingleMuon Run2017 A/B");
  gStyle->SetPalette(kLightTemperature);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetPaintTextFormat("3.2f");
  // Draw the 2D plots
  TPaveText *selectionPave2D = new TPaveText(0.7,0.9,1,1, "NDC");
  selectionPave2D->SetFillColorAlpha(0,0.4);
  selectionPave2D->SetFillStyle(0);
  selectionPave2D->SetTextColor(14);
  selectionPave2D->SetBorderSize(0);
  if     (mode=="monojet") selectionPave2D->AddText("#bf{#it{#splitline{\"Monojet\" topology}{Central jet, p_{T}>100, r_{CH}>10%, r_{NH}<80%}}}");
  else if(mode=="vbf"    ) selectionPave2D->AddText("#bf{#it{#splitline{VBF topology}{Two jets p_{T}>80(40), leading jet r_{CH}>10%, r_{NH}<80% if central}}}");
  
  TCanvas *cEff_metPhi_met=new TCanvas("cEff_metPhi_met", "cEff_metPhi_met",800,600);
  cEff_metPhi_met->SetLogx();
  hEff_metPhi_met->GetXaxis()->SetTitle(flavor=="electrons"?"E_{T}^{miss} [GeV]":"E_{T}^{miss} + p_{T}^{#mu} [GeV]");
  hEff_metPhi_met->GetXaxis()->SetMoreLogLabels();
  hEff_metPhi_met->GetYaxis()->SetTitle(flavor=="electrons"?"#phi(E_{T}^{miss}) [rad]":"#phi(E_{T}^{miss} + p_{T}^{#mu}) [rad]");
  hEff_metPhi_met->SetMarkerColor(kBlack);
  hEff_metPhi_met->SetMarkerSize(1.2);
  hEff_metPhi_met->SetTitle(title);
  hEff_metPhi_met->Draw("COLZ TEXTE");
  cEff_metPhi_met->Modified(); cEff_metPhi_met->Update();
  ((TPaveText*)cEff_metPhi_met->GetPrimitive("title"))->SetX1NDC(0);
  ((TPaveText*)cEff_metPhi_met->GetPrimitive("title"))->SetX2NDC(0.7);
  ((TPaveText*)cEff_metPhi_met->GetPrimitive("title"))->SetY1NDC(0.9);
  ((TPaveText*)cEff_metPhi_met->GetPrimitive("title"))->SetY2NDC(1);
  cEff_metPhi_met->Modified(); cEff_metPhi_met->Update();
  selectionPave2D->Draw("SAME");
  
  TCanvas *cEff_jet1Pt_met=new TCanvas("cEff_jet1Pt_met", "cEff_jet1Pt_met",800,600);
  cEff_jet1Pt_met->SetLogx();
  cEff_jet1Pt_met->SetLogy();
  hEff_jet1Pt_met->GetXaxis()->SetTitle(flavor=="electrons"?"E_{T}^{miss} [GeV]":"E_{T}^{miss} + p_{T}^{#mu} [GeV]");
  hEff_jet1Pt_met->GetXaxis()->SetMoreLogLabels();
  hEff_jet1Pt_met->GetYaxis()->SetTitle("Jet 1 p_{T} [GeV]");
  hEff_jet1Pt_met->GetYaxis()->SetMoreLogLabels();
  hEff_jet1Pt_met->GetYaxis()->SetTitleOffset(1.4);
  hEff_jet1Pt_met->SetMarkerColor(kBlack);
  hEff_jet1Pt_met->SetMarkerSize(1.2);
  hEff_jet1Pt_met->SetTitle(title);
  hEff_jet1Pt_met->Draw("COLZ TEXTE");
  if(mode=="monojet") hEff_jet1Pt_met->GetYaxis()->SetRangeUser(100., leptonPtBins[nBinsLeptonPt]);
  cEff_jet1Pt_met->Modified(); cEff_jet1Pt_met->Update();
  ((TPaveText*)cEff_jet1Pt_met->GetPrimitive("title"))->SetX1NDC(0);
  ((TPaveText*)cEff_jet1Pt_met->GetPrimitive("title"))->SetX2NDC(0.7);
  ((TPaveText*)cEff_jet1Pt_met->GetPrimitive("title"))->SetY1NDC(0.9);
  ((TPaveText*)cEff_jet1Pt_met->GetPrimitive("title"))->SetY2NDC(1);
  cEff_jet1Pt_met->Modified(); cEff_jet1Pt_met->Update();
  selectionPave2D->Draw("SAME");
  
  TCanvas *cEff_jet1Eta_met=new TCanvas("cEff_jet1Eta_met", "cEff_jet1Eta_met",800,600);
  cEff_jet1Eta_met->SetLogx();
  hEff_jet1Eta_met->GetXaxis()->SetTitle(flavor=="electrons"?"E_{T}^{miss} [GeV]":"E_{T}^{miss} + p_{T}^{#mu} [GeV]");
  hEff_jet1Eta_met->GetXaxis()->SetMoreLogLabels();
  hEff_jet1Eta_met->GetYaxis()->SetTitle("Jet 1 pseudorapidity");
  hEff_jet1Eta_met->SetMarkerColor(kBlack);
  hEff_jet1Eta_met->SetMarkerSize(1.2);
  hEff_jet1Eta_met->SetTitle(title);
  hEff_jet1Eta_met->Draw("COLZ TEXTE");
  if(mode=="monojet") hEff_jet1Eta_met->GetYaxis()->SetRangeUser(-2.499,2.499);
  cEff_jet1Eta_met->Modified(); cEff_jet1Eta_met->Update();
  ((TPaveText*)cEff_jet1Eta_met->GetPrimitive("title"))->SetX1NDC(0);
  ((TPaveText*)cEff_jet1Eta_met->GetPrimitive("title"))->SetX2NDC(0.7);
  ((TPaveText*)cEff_jet1Eta_met->GetPrimitive("title"))->SetY1NDC(0.9);
  ((TPaveText*)cEff_jet1Eta_met->GetPrimitive("title"))->SetY2NDC(1);
  cEff_jet1Eta_met->Modified(); cEff_jet1Eta_met->Update();
  selectionPave2D->Draw("SAME");
  
  TCanvas *cEff_jet2Eta_met, *cEff_mjj_met;
  if(mode=="vbf") {
    cEff_jet2Eta_met=new TCanvas("cEff_jet2Eta_met", "cEff_jet2Eta_met",800,600);
    cEff_jet2Eta_met->SetLogx();
    hEff_jet2Eta_met->GetXaxis()->SetTitle(flavor=="electrons"?"E_{T}^{miss} [GeV]":"E_{T}^{miss} + p_{T}^{#mu} [GeV]");
    hEff_jet2Eta_met->GetXaxis()->SetMoreLogLabels();
    hEff_jet2Eta_met->GetYaxis()->SetTitle("Jet 2 pseudorapidity");
    hEff_jet2Eta_met->SetMarkerColor(kBlack);
    hEff_jet2Eta_met->SetMarkerSize(1.2);
    hEff_jet2Eta_met->SetTitle(title);
    hEff_jet2Eta_met->Draw("COLZ TEXTE");
    cEff_jet2Eta_met->Modified(); cEff_jet2Eta_met->Update();
    ((TPaveText*)cEff_jet2Eta_met->GetPrimitive("title"))->SetX1NDC(0);
    ((TPaveText*)cEff_jet2Eta_met->GetPrimitive("title"))->SetX2NDC(0.7);
    ((TPaveText*)cEff_jet2Eta_met->GetPrimitive("title"))->SetY1NDC(0.9);
    ((TPaveText*)cEff_jet2Eta_met->GetPrimitive("title"))->SetY2NDC(1);
    cEff_jet2Eta_met->Modified(); cEff_jet2Eta_met->Update();
    selectionPave2D->Draw("SAME");
    
    cEff_mjj_met=new TCanvas("cEff_mjj_met", "cEff_mjj_met",800,600);
    cEff_mjj_met->SetLogx();
    hEff_mjj_met->GetXaxis()->SetTitle(flavor=="electrons"?"E_{T}^{miss} [GeV]":"E_{T}^{miss} + p_{T}^{#mu} [GeV]");
    hEff_mjj_met->GetXaxis()->SetMoreLogLabels();
    hEff_mjj_met->GetYaxis()->SetTitle("Dijet mass [GeV]");
    hEff_mjj_met->GetYaxis()->SetTitleOffset(1.4);
    hEff_mjj_met->SetMarkerColor(kBlack);
    hEff_mjj_met->SetMarkerSize(1.2);
    hEff_mjj_met->SetTitle(title);
    hEff_mjj_met->Draw("COLZ TEXTE");
    cEff_mjj_met->Modified(); cEff_mjj_met->Update();
    ((TPaveText*)cEff_mjj_met->GetPrimitive("title"))->SetX1NDC(0);
    ((TPaveText*)cEff_mjj_met->GetPrimitive("title"))->SetX2NDC(0.7);
    ((TPaveText*)cEff_mjj_met->GetPrimitive("title"))->SetY1NDC(0.9);
    ((TPaveText*)cEff_mjj_met->GetPrimitive("title"))->SetY2NDC(1);
    cEff_mjj_met->Modified(); cEff_mjj_met->Update();
    selectionPave2D->Draw("SAME");
  }
  TCanvas *cEff_leptonPt_met=new TCanvas("cEff_leptonPt_met", "cEff_leptonPt_met",800,600);
  cEff_leptonPt_met->SetLogx();
  cEff_leptonPt_met->SetLogy();
  hEff_leptonPt_met->GetXaxis()->SetTitle(flavor=="electrons"?"E_{T}^{miss} [GeV]":"E_{T}^{miss} + p_{T}^{#mu} [GeV]");
  hEff_leptonPt_met->GetXaxis()->SetMoreLogLabels();
  hEff_leptonPt_met->GetYaxis()->SetTitle(flavor=="electrons"?"Electron p_{T} [GeV]":"Muon p_{T} [GeV]");
  hEff_leptonPt_met->GetYaxis()->SetMoreLogLabels();
  hEff_leptonPt_met->GetYaxis()->SetTitleOffset(1.4);
  hEff_leptonPt_met->SetMarkerColor(kBlack);
  hEff_leptonPt_met->SetMarkerSize(1.2);
  hEff_leptonPt_met->SetTitle(title);
  hEff_leptonPt_met->Draw("COLZ TEXTE");
  cEff_leptonPt_met->Modified(); cEff_leptonPt_met->Update();
  ((TPaveText*)cEff_leptonPt_met->GetPrimitive("title"))->SetX1NDC(0);
  ((TPaveText*)cEff_leptonPt_met->GetPrimitive("title"))->SetX2NDC(0.7);
  ((TPaveText*)cEff_leptonPt_met->GetPrimitive("title"))->SetY1NDC(0.9);
  ((TPaveText*)cEff_leptonPt_met->GetPrimitive("title"))->SetY2NDC(1);
  cEff_leptonPt_met->Modified(); cEff_leptonPt_met->Update();
  selectionPave2D->Draw("SAME");
  cEff_leptonPt_met->RedrawAxis();

  // Draw the 1D plot with the fit
  TCanvas *cEff_met=new TCanvas("cEff_met", "cEff_met",1024,640);
  cEff_met->SetRightMargin(0.3);
  hEff_met->GetXaxis()->SetTitle(flavor=="electrons"?"E_{T}^{miss} [GeV]":"E_{T}^{miss} + p_{T}^{#mu} [GeV]");
  hEff_met->GetXaxis()->SetMoreLogLabels();
  hEff_met->GetYaxis()->SetTitle("Efficiency");
  hEff_met->SetMarkerStyle(20);
  hEff_met->SetMarkerSize(0.8);
  hEff_met->Draw("AP E1");
  hEff_met->SetTitle(title);
  hEff_met->GetYaxis()->SetRangeUser(0,1.2);
  //hEff_met->GetXaxis()->SetRangeUser(hTotal_met->GetXaxis()->GetBinLowEdge(1), hTotal_met->GetXaxis()->GetBinUpEdge(hTotal_met->GetNbinsX())-.001);
  hEff_met->GetXaxis()->SetRangeUser(0,1199.9);
  //TLine *oneline = new TLine(hTotal_met->GetXaxis()->GetBinLowEdge(1), 1, hTotal_met->GetXaxis()->GetBinUpEdge(hTotal_met->GetNbinsX()), 1);
  TLine *oneline = new TLine(0, 1, 1199.9, 1);
  oneline->SetLineStyle(2);
  oneline->Draw("SAME");
  TF1 *f_sigmoid;
  if     (fitFunction==1) { f_sigmoid = new TF1("sigmoidFit", "[2]*0.5*(1 + TMath::TanH( [0]*( x-[1] )))", 0, 1500); f_sigmoid->SetLineColor(kGreen); } 
  else if(fitFunction==2) { f_sigmoid = new TF1("sigmoidFit", "[2]*(0.5 + 1./TMath::Pi()*TMath::ATan( TMath::Pi()*0.5*[0]*( x-[1] )))", 0, 1500); f_sigmoid->SetLineColor(kRed); }
  else if(fitFunction==3) { f_sigmoid = new TF1("sigmoidFit", "[2]*0.5*(1+TMath::Erf( TMath::Sqrt(TMath::Pi())*0.5*[0]*( x-[1] )))", 0, 1500); f_sigmoid->SetLineColor(kOrange+4); }
  else if(fitFunction==4) { f_sigmoid = new TF1("sigmoidFit", "[2]*pow(1+[3]*TMath::Exp(-[0]*(x-[1])),-1./[4])", 0, 1500); f_sigmoid->SetLineColor(kMagenta); }
  if(fitFunction==1 || fitFunction==2 || fitFunction==3) {
    f_sigmoid->SetParNames("k", "E_{c}","#varepsilon_{obs}");
    f_sigmoid->SetParameter(0,0.03); f_sigmoid->SetParLimits(0,0,100);
    f_sigmoid->SetParameter(1,160);   f_sigmoid->SetParLimits(1,100,300);
    f_sigmoid->SetParameter(2,.95);  f_sigmoid->SetParLimits(2,0,1); 
    //f_sigmoid->SetParameter(2,1);  f_sigmoid->SetParLimits(2,1,1); 
  } else if(fitFunction==4) {
    f_sigmoid->SetParNames("k", "E_{c}","#varepsilon_{obs}","Q","#nu");
    f_sigmoid->SetParameter(0,.03);  f_sigmoid->SetParLimits(0,0,100);
    f_sigmoid->SetParameter(1,160);  f_sigmoid->SetParLimits(1,100,300);
    f_sigmoid->SetParameter(2,.95);  f_sigmoid->SetParLimits(2,0,1); 
    f_sigmoid->SetParameter(3,1);    f_sigmoid->SetParLimits(3,0,10); 
    f_sigmoid->SetParameter(4,1);    f_sigmoid->SetParLimits(4,0,100); 
  } else {
    printf("bad fit function\n");
    return;
  }
  hEff_met->Fit(f_sigmoid, "ME","", 0,1500);
  TPaveText *formula = new TPaveText(0.72,0.8,0.96,0.9, "NDC");
  formula->SetFillColorAlpha(0,0.4);
  formula->SetFillStyle(1001);
  formula->SetLineColor(kBlack);
  formula->SetBorderSize(1);
  if     (fitFunction==1) formula->AddText("f(E_{T}^{miss}) = #varepsilon_{obs}#times #frac{1}{2} #(){1 + tanh(k[E_{T}^{miss}- E_{c}])} ");
  else if(fitFunction==2) formula->AddText("f(E_{T}^{miss}) = #varepsilon_{obs}#times #(){#frac{1}{2} + #frac{1}{#pi} arctan(#frac{k#pi}{2}[E_{T}^{miss}- E_{c}])} ");
  else if(fitFunction==3) formula->AddText("f(E_{T}^{miss}) = #varepsilon_{obs}#times #frac{1}{2} #(){1 + erf(#frac{k#sqrt{#pi}}{2}[E_{T}^{miss}- E_{c}])} ");
  else if(fitFunction==4) formula->AddText("f(E_{T}^{miss}) = #varepsilon_{obs}#times #[]{1 + Q exp(-k[E_{T}^{miss}- E_{c}])}^{-1/#nu}");
  formula->Draw("SAME");
  TPaveText *selectionPave = new TPaveText(0.12,0.77,.4,0.88, "NDC");
  selectionPave->SetFillColorAlpha(0,0.4);
  selectionPave->SetFillStyle(1001);
  selectionPave->SetTextColor(14);
  selectionPave->SetBorderSize(0);
  if     (mode=="monojet") selectionPave->AddText("#bf{#it{#splitline{\"Monojet\" topology}{Central jet, p_{T}>100, r_{CH}>10%, r_{NH}<80%}}}");
  else if(mode=="vbf"    ) selectionPave->AddText("#bf{#it{#splitline{VBF topology}{Two jets p_{T}>80(40), leading jet r_{CH}>10%, r_{NH}<80% if central}}}");
  selectionPave->Draw("SAME");
  cEff_met->Modified(); cEff_met->Update(); 
  TPaveStats *ps;
  ps  = (TPaveStats *)cEff_met->GetPrimitive("stats");
  if(!ps) { printf("bad pointer\n"); return; }
  ps->SetX1NDC(0.72);
  ps->SetX2NDC(0.96);
  ps->SetY1NDC(0.4);
  ps->SetY2NDC(0.78);
  cEff_met->Modified(); cEff_met->Update(); 
  
  system("mkdir -p MitVBFAnalysis/plots");
  cEff_metPhi_met   ->SaveAs(Form("MitVBFAnalysis/plots/hEff_metPhi_met_%s.pdf"   ,extra_string.Data()));     
  cEff_jet1Pt_met   ->SaveAs(Form("MitVBFAnalysis/plots/hEff_jet1Pt_met_%s.pdf"   ,extra_string.Data()));     
  cEff_jet1Eta_met  ->SaveAs(Form("MitVBFAnalysis/plots/hEff_jet1Eta_met_%s.pdf"  ,extra_string.Data()));      
  cEff_leptonPt_met ->SaveAs(Form("MitVBFAnalysis/plots/hEff_leptonPt_met_%s.pdf" ,extra_string.Data()));       
  cEff_met          ->SaveAs(Form("MitVBFAnalysis/plots/hEff_met_%s_%s.pdf"       ,extra_string.Data(), fitFunctionName.Data())); 
  if(mode=="vbf") cEff_jet2Eta_met  ->SaveAs(Form("MitVBFAnalysis/plots/hEff_jet2Eta_met_%s.pdf"  ,extra_string.Data()));      
  if(mode=="vbf") cEff_mjj_met      ->SaveAs(Form("MitVBFAnalysis/plots/hEff_mjj_met_%s.pdf"      ,extra_string.Data()));  
  
  cEff_metPhi_met   ->SaveAs(Form("MitVBFAnalysis/plots/hEff_metPhi_met_%s.png"   ,extra_string.Data()));     
  cEff_jet1Pt_met   ->SaveAs(Form("MitVBFAnalysis/plots/hEff_jet1Pt_met_%s.png"   ,extra_string.Data()));     
  cEff_jet1Eta_met  ->SaveAs(Form("MitVBFAnalysis/plots/hEff_jet1Eta_met_%s.png"  ,extra_string.Data()));      
  cEff_leptonPt_met ->SaveAs(Form("MitVBFAnalysis/plots/hEff_leptonPt_met_%s.png" ,extra_string.Data()));       
  cEff_met          ->SaveAs(Form("MitVBFAnalysis/plots/hEff_met_%s_%s.png"       ,extra_string.Data(), fitFunctionName.Data())); 
  if(mode=="vbf") cEff_jet2Eta_met  ->SaveAs(Form("MitVBFAnalysis/plots/hEff_jet2Eta_met_%s.png"  ,extra_string.Data()));      
  if(mode=="vbf") cEff_mjj_met      ->SaveAs(Form("MitVBFAnalysis/plots/hEff_mjj_met_%s.png"      ,extra_string.Data()));  

}

/*
void mitPalette()
{
  static Int_t  colors[100];
  static Bool_t initialized = kFALSE;
  Double_t Red[3]    = { 1, 138./255., 163/255.};
  Double_t Green[3]  = { 1, 139./255., 31/255.};
  Double_t Blue[3]   = { 1, 140./255., 52/255.};
  Double_t Length[3] = { 0.00, 0.35, 1.00 };
  if(!initialized){
    Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,100);
    for (int i=0; i<100; i++) colors[i] = FI+i;
    initialized = kTRUE;
    return;
  }
  gStyle->SetPalette(100,colors);

}
void mitPalette2()
{
  static Int_t  colors[100];
  static Bool_t initialized = kFALSE;
  Double_t Red[3]    = { 138/255., 1., 163/255.};
  Double_t Green[3]  = { 139/255., 1., 31/255.};
  Double_t Blue[3]   = { 140/255., 1., 52/255.};
  Double_t Length[3] = { 0.00, 0.5, 1.00 };
  if(!initialized){
    Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,100);
    for (int i=0; i<100; i++) colors[i] = FI+i;
    initialized = kTRUE;
    return;
  }
  gStyle->SetPalette(100,colors);

}
void mitPalette3()
{
  static Int_t  colors[100];
  static Bool_t initialized = kFALSE;
  Double_t Red[3]    = { 163/255., 138./255., 1   };
  Double_t Green[3]  = { 31/255. , 139./255., 1   };
  Double_t Blue[3]   = { 52/255. , 140./255., 1   };
  Double_t Length[3] = { 0.00    , 0.65,      1.00};
  if(!initialized){
    Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,100);
    for (int i=0; i<100; i++) colors[i] = FI+i;
    initialized = kTRUE;
    return;
  }
  gStyle->SetPalette(100,colors);

}
*/
