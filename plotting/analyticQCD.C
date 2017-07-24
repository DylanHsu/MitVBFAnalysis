void analyticQCD(bool drawFits=true, double E1=50., double E2=1000.) { 
  bool integralConstraint=true;
  const int nBinsMjj=9, nMetPars=6;
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  system("mkdir -p MitVBFAnalysis/plots");

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls( 1e8);
  TFile *file = TFile::Open("MitVBFAnalysis/data/qcd_met_mjj.root", "READ");
  TH2F *histo_met_mjj = ((TH2F*) file->Get("h_mjj_met")->Clone("histo_met_mjj")); histo_met_mjj->SetDirectory(0);
  histo_met_mjj->Rebin2D(1,1); histo_met_mjj->GetYaxis()->SetRangeUser(E1,E2);
  TF1 *doubleExp;
  if(integralConstraint) { 
    doubleExp = new TF1("doubleExp", "([0]-[2]/[3]*(TMath::Exp(-[3]*[6])-TMath::Exp(-[3]*[7]))-[4]/[5]*(TMath::Exp(-[5]*[6])-TMath::Exp(-[5]*[7])))*([1]/(TMath::Exp(-[1]*[6])-TMath::Exp(-[1]*[7])))*TMath::Exp(-[1]*x) + [2]*TMath::Exp(-[3]*x) + [4]*TMath::Exp(-[5]*x)",E1,E2);
    doubleExp->SetParNames("Integral", "k_{1}", "A_{2}", "k_{2}", "A_{3}", "k_{3}", "E_{1}", "E_{2}");
  } else {
    doubleExp = new TF1("doubleExp", "[0]*TMath::Exp(-[1]*x) + [2]*TMath::Exp(-[3]*x+[4]*x*x)",E1,E2);
    doubleExp->SetParNames("A_{1}", "k_{1}", "A_{2}", "k_{2}", "q");
  }
  TCanvas *c1=new TCanvas ("c1","c1");
  c1->SetLogy();
  TH1D *histo_projections[nBinsMjj], *histo_parameters[nMetPars];
  TCanvas *canvases[nBinsMjj];
  for(int iPar=0; iPar<nMetPars; iPar++) {
    histo_parameters[iPar] = histo_met_mjj->ProjectionX(Form("h_%s",doubleExp->GetParName(iPar)));
    histo_parameters[iPar]->Clear();
  }
  for(int nBinMjj=1; nBinMjj<=nBinsMjj; nBinMjj++) {
    histo_projections[nBinMjj-1] = histo_met_mjj->ProjectionY(Form("metProjection_mjj%dto%d",(int)histo_met_mjj->GetXaxis()->GetBinLowEdge(nBinMjj),(int)histo_met_mjj->GetXaxis()->GetBinLowEdge(nBinMjj+1)), nBinMjj, nBinMjj, "e");
    double integral=0; 
    for(unsigned int i=histo_projections[nBinMjj-1]->FindBin(E1); i<=histo_projections[nBinMjj-1]->FindBin(E2-0.001); i++) {
      histo_projections[nBinMjj-1]->SetBinContent(i, histo_projections[nBinMjj-1]->GetBinContent(i)/histo_projections[nBinMjj-1]->GetBinWidth(i)); 
      // FIX: assign a 10% "theory" uncertainty to each bin (very bad hack)
      double fakeUnc=0.0;
      histo_projections[nBinMjj-1]->SetBinError(i, sqrt(pow(histo_projections[nBinMjj-1]->GetBinError(i)/histo_projections[nBinMjj-1]->GetBinWidth(i),2) + pow(fakeUnc*histo_projections[nBinMjj-1]->GetBinContent(i),2)));
      
      integral+=histo_projections[nBinMjj-1]->GetBinContent(i)*histo_projections[nBinMjj-1]->GetBinWidth(i);
    }
    if(integralConstraint) {
      double A2_estimation = TMath::Exp(0.02*300.)*histo_projections[nBinMjj-1]->GetBinContent(histo_projections[nBinMjj-1]->FindBin(300.));
      double A3_estimation = TMath::Exp(0.01*700.)*histo_projections[nBinMjj-1]->GetBinContent(histo_projections[nBinMjj-1]->FindBin(700.));
      doubleExp->FixParameter(0, integral);
      doubleExp->SetParameter(1, 0.05);
      doubleExp->SetParLimits(1, 0.00001, 0.5);
      doubleExp->SetParameter(2, A2_estimation );
      doubleExp->SetParLimits(2, 0.01*A2_estimation, 100.*A2_estimation);
      doubleExp->SetParameter(3, 0.02);
      doubleExp->SetParLimits(3, 0.00001, 0.5);
      doubleExp->SetParameter(4, A3_estimation );
      doubleExp->SetParLimits(4, 0, A2_estimation);
      doubleExp->SetParameter(5, 0.005);
      doubleExp->SetParLimits(5, 0.00001, 0.5);
      doubleExp->FixParameter(6, E1);
      doubleExp->FixParameter(7, E2);
    } else {
      doubleExp->SetParameter(0, integral);
      doubleExp->SetParLimits(0, .01*integral, 100*integral);
      doubleExp->SetParameter(1, 0.05);
      doubleExp->SetParLimits(1, 0.001, 0.2);
      doubleExp->SetParameter(2, 0.01*integral);
      doubleExp->SetParLimits(2, 0.0001*integral, integral);
      doubleExp->SetParameter(3, 0.02);
      doubleExp->SetParLimits(3, 0.0001, 0.2);
      doubleExp->SetParameter(4, 1e-5);
    }
    if(drawFits) {
      canvases[nBinsMjj-1] = new TCanvas(Form("c_metProjection_mjj%dto%d",(int)histo_met_mjj->GetXaxis()->GetBinLowEdge(nBinMjj),(int)histo_met_mjj->GetXaxis()->GetBinLowEdge(nBinMjj+1)), 
      Form("c_metProjection_mjj%dto%d",(int)histo_met_mjj->GetXaxis()->GetBinLowEdge(nBinMjj),(int)histo_met_mjj->GetXaxis()->GetBinLowEdge(nBinMjj+1)));
      histo_projections[nBinMjj-1]->Draw();
    }
    histo_projections[nBinMjj-1]->Fit(doubleExp,"IME","",E1,E2);
    printf("histo integral %f, function integral %f\n", integral, doubleExp->Integral(E1,E2));
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
  if(!integralConstraint) {
    TCanvas *c_A1 = new TCanvas("c_A1", "c_A1");
    c_A1->SetLogy();
    histo_parameters[0]->SetTitle("Parameter A_{1} vs. dijet mass");
    histo_parameters[0]->SetMinimum(1);
    histo_parameters[0]->SetMaximum(1e10);
    histo_parameters[0]->Draw();
    c_A1->Update();
    c_A1->Print("MitVBFAnalysis/plots/c_A1.png");
  }
  TCanvas *c_A2 = new TCanvas("c_A2", "c_A2");
  c_A2->SetLogy();
  histo_parameters[2]->SetTitle("Parameter A_{2} vs. dijet mass");
  histo_parameters[2]->SetMinimum(1e-2);
  histo_parameters[2]->SetMaximum(5e2);
  histo_parameters[2]->Draw();
  histo_parameters[2]->GetXaxis()->SetTitle("GeV");
  c_A2->Update();
  c_A2->Print("MitVBFAnalysis/plots/c_A2.png");
  TCanvas *c_k1 = new TCanvas("c_k1", "c_k1");
  histo_parameters[1]->SetTitle("Parameter k_{1} vs. dijet mass");
  histo_parameters[1]->SetMinimum(0.);
  histo_parameters[1]->SetMaximum(0.12);
  histo_parameters[1]->Draw();
  histo_parameters[1]->GetXaxis()->SetTitle("GeV");
  c_k1->Update();
  c_k1->Print("MitVBFAnalysis/plots/c_k1.png");
  TCanvas *c_k2 = new TCanvas("c_k2", "c_k2");
  histo_parameters[3]->SetTitle("Parameter k_{2} vs. dijet mass");
  histo_parameters[3]->SetMinimum(0.);
  histo_parameters[3]->SetMaximum(0.03);
  histo_parameters[3]->Draw();
  histo_parameters[3]->GetXaxis()->SetTitle("GeV");
  c_k2->Update();
  c_k2->Print("MitVBFAnalysis/plots/c_k2.png");
}
