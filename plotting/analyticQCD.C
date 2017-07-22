void analyticQCD(bool drawFits=true) { 
  const int nBinsMjj=50, nMetPars=4;
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  system("mkdir -p MitVBFAnalysis/plots");

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls( 1e8);
  TFile *file = TFile::Open("/data/t3home000/snarayan/store/panda/vbf_004_qcd/forD.root", "READ");
  TH2F *histo_met_mjj = ((TH2F*) file->Get("h")->Clone("histo_met_mjj")); histo_met_mjj->SetDirectory(0);
  histo_met_mjj->Rebin2D(1,1); histo_met_mjj->GetYaxis()->SetRangeUser(48,400);
  TF1 *doubleExp = new TF1("doubleExp", "[0]*TMath::Exp(-[1]*x) + [2]*TMath::Exp(-[3]*x)",48,400);
  doubleExp->SetParNames("A_{1}", "k_{1}", "A_{2}", "k_{2}");
  TCanvas *c1=new TCanvas ("c1","c1");
  c1->SetLogy();
  TH1D *histo_projections[nBinsMjj], *histo_parameters[nMetPars];
  TCanvas *canvases[nBinsMjj];
  for(int iPar=0; iPar<nMetPars; iPar++) {
    histo_parameters[iPar] = histo_met_mjj->ProjectionX(Form("h_%s",doubleExp->GetParName(iPar)));
    histo_parameters[iPar]->Clear();
  }
  for(int nBinMjj=1; nBinMjj<nBinsMjj; nBinMjj++) {
    histo_projections[nBinMjj-1] = histo_met_mjj->ProjectionY(Form("metProjection_mjj%dto%d",(int)histo_met_mjj->GetXaxis()->GetBinLowEdge(nBinMjj),(int)histo_met_mjj->GetXaxis()->GetBinLowEdge(nBinMjj+1)), nBinMjj);
  
    // Assign a 10% "theory" uncertainty to each bin (very bad hack)
    for(unsigned int i=1; i<=histo_projections[nBinMjj-1]->GetNbinsX(); i++) {
      histo_projections[nBinMjj-1]->SetBinError(i, sqrt(pow(histo_projections[nBinMjj-1]->GetBinError(i),2) + pow(0.10*histo_projections[nBinMjj-1]->GetBinContent(i),2)));
    }
    double integral=histo_projections[nBinMjj-1]->Integral();
    doubleExp->SetParameter(0, integral);
    doubleExp->SetParameter(1, 0.05);
    doubleExp->SetParameter(2, 0.001*integral);
    doubleExp->SetParameter(3, 0.02);
    doubleExp->SetParLimits(0, integral, integral*integral*integral);
    doubleExp->SetParLimits(1, 0.001, 0.2);
    doubleExp->SetParLimits(2, 0.01, 0.95*integral);
    doubleExp->SetParLimits(3, 0.0001, 0.2);
    if(drawFits) {
      canvases[nBinsMjj-1] = new TCanvas(Form("c_metProjection_mjj%dto%d",(int)histo_met_mjj->GetXaxis()->GetBinLowEdge(nBinMjj),(int)histo_met_mjj->GetXaxis()->GetBinLowEdge(nBinMjj+1)), 
      Form("c_metProjection_mjj%dto%d",(int)histo_met_mjj->GetXaxis()->GetBinLowEdge(nBinMjj),(int)histo_met_mjj->GetXaxis()->GetBinLowEdge(nBinMjj+1)));
      histo_projections[nBinMjj-1]->Draw();
    }
    histo_projections[nBinMjj-1]->Fit(doubleExp,"ME","",48,400);
    if(drawFits) {
      canvases[nBinsMjj-1]->Print( Form("MitVBFAnalysis/plots/c_metProjection_mjj%dto%d.png",(int)histo_met_mjj->GetXaxis()->GetBinLowEdge(nBinMjj),(int)histo_met_mjj->GetXaxis()->GetBinLowEdge(nBinMjj+1)));
      canvases[nBinsMjj-1]->SetLogy();
      canvases[nBinsMjj-1]->Update();
      canvases[nBinsMjj-1]->Print( Form("MitVBFAnalysis/plots/c_log_metProjection_mjj%dto%d.png",(int)histo_met_mjj->GetXaxis()->GetBinLowEdge(nBinMjj),(int)histo_met_mjj->GetXaxis()->GetBinLowEdge(nBinMjj+1)));
      delete canvases[nBinsMjj-1];
    }
    for(int iPar=0; iPar<nMetPars; iPar++) {
      histo_parameters[iPar]->SetBinContent(nBinMjj, doubleExp->GetParameter(iPar));
      histo_parameters[iPar]->SetBinError(nBinMjj, doubleExp->GetParError(iPar));
    }
  }
  
  TCanvas *c_A1 = new TCanvas("c_A1", "c_A1");
  c_A1->SetLogy();
  histo_parameters[0]->SetMinimum(1);
  histo_parameters[0]->SetMaximum(1e10);
  histo_parameters[0]->Draw();
  // need to plot other parameters
}
