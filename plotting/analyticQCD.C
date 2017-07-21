void analyticQCD(int bin=1) { 
  const int nBinsMjj=50;
  gStyle->SetOptFit(1111);

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls( 1e8);
  TFile *file = TFile::Open("/data/t3home000/snarayan/store/panda/vbf_004_qcd/forD.root", "READ");
  TH2F *histo_met_mjj = ((TH2F*) file->Get("h")->Clone("histo_met_mjj")); histo_met_mjj->SetDirectory(0);
  histo_met_mjj->Rebin2D(1,1); histo_met_mjj->GetYaxis()->SetRangeUser(48,400);
  //TF1 *doubleExp = new TF1("doubleExp", "[0]*TMath::Exp(-[1]*(x-[2]))+[3]*TMath::Exp(-[4]*(x-[5]))",50,400);
  TF1 *doubleExp = new TF1("doubleExp", "[0]*TMath::Exp(-[1]*x) + [2]*TMath::Exp(-[3]*x)",48,400);
  //doubleExp->SetParNames("A_{1}", "k_{1}", "E_{1}", "A_{2}", "k_{2}", "E_{2}");
  doubleExp->SetParNames("A_{1}", "k_{1}", "A_{2}", "k_{2}");
  TCanvas *c1=new TCanvas ("c1","c1");
  c1->SetLogy();
  
  //TH1D *projections[nBinsMjj];
  //TCanvas *canvases[nBinsMjj];
  //for(int nBinMjj=1; nBinMjj<nBinsMjj; nBinMjj++) {
  TH1D *projection=histo_met_mjj->ProjectionY("projection",bin);
  
  // Assign a 10% "theory" uncertainty to each bin (very bad hack)
  for(unsigned int i=1; i<=projection->GetNbinsX(); i++) {
    projection->SetBinError(i, sqrt(pow(projection->GetBinError(i),2) + pow(0.10*projection->GetBinContent(i),2)));
  }
  double integral=projection->Integral();
  doubleExp->SetParameter(0, integral);
  doubleExp->SetParameter(1, 0.05);
  doubleExp->SetParameter(2, 0.001*integral);
  doubleExp->SetParameter(3, 0.02);
  doubleExp->SetParLimits(0, integral, 100.*integral);
  doubleExp->SetParLimits(1, 0.001, 0.2);
  doubleExp->SetParLimits(2, 0, 0.95*integral);
  doubleExp->SetParLimits(3, 0.0001, 0.2);
  projection->Draw();
  projection->Fit(doubleExp,"ME","",48,400);

  histo_met_mjj->FitSlicesY(doubleExp, 0,-1, 0, "N");
  //TObjArray aSlices;
  //histo_met_mjj->FitSlicesY(doubleExp, 0,-1, 0, "QNR", &aSlices);
}
