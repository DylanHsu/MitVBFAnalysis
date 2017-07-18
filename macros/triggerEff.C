#include "MitVBFAnalysis/macros/triggerEff.h"

// Last modified by DGH July 18 2017
// This code produces a flat tree for trigger studies from Panda files (provided as inputPath var)
// Compiles through ACLiC and runs on condor via the shell scripts
// Analyze the output using plotting scripts in  MitVBFAnalysis/plotting

typedef std::map<UInt_t,std::vector<std::pair <UInt_t, UInt_t> > > MapType;

void triggerEff(
TString mode="electrons",
TString inputPath="",
TString batchName="",
bool type1=false,
bool debug=false
) {
  string jsonFile = "PandaAnalysis/data/certs/Cert_294927-297723_13TeV_PromptReco_Collisions17_JSON.txt";
  assert(mode=="electrons" || mode=="muons");
  
  // initialize output tree
  TFile *outputFile = TFile::Open(("vbf_batchTree_triggerEff"+(batchName == "" ? "" : "_"+batchName)+".root"), "RECREATE", "", ROOT::CompressionSettings(ROOT::kZLIB,9));
  bool passMetTriggers, passJetMetCleaning, jet1ExtraClean, jet2ExtraClean, jet3ExtraClean, jet4ExtraClean;
  float leptonPt, leptonEta, leptonPhi, pfMetNoMu, pfMetNoMuPhi, jet1Pt, jet1Eta, jet1Phi, jet1M, jet1PtUp, jet1PtDown, jet2Pt, jet2Eta, jet2Phi, jet2M, jet2PtUp, jet2PtDown, jet3Pt, jet3Eta, jet3Phi, jet3M, jet3PtUp, jet3PtDown, jet4Pt, jet4Eta, jet4Phi, jet4M, jet4PtUp, jet4PtDown;
  unsigned nTau;
  UInt_t runNumber, lumiSection;
  ULong64_t eventNumber;
  TTree *effTree = new TTree("effTree", "effTree");
  effTree->Branch("runNumber"          , &runNumber         );
  effTree->Branch("lumiSection"        , &lumiSection       );
  effTree->Branch("eventNumber"        , &eventNumber       );
  effTree->Branch("passMetTriggers"    , &passMetTriggers   );
  effTree->Branch("passJetMetCleaning" , &passJetMetCleaning);
  effTree->Branch("leptonPt"           , &leptonPt          );
  effTree->Branch("leptonEta"          , &leptonEta         );
  effTree->Branch("leptonPhi"          , &leptonPhi         );
  effTree->Branch("pfMetNoMu"          , &pfMetNoMu         );
  effTree->Branch("pfMetNoMuPhi"       , &pfMetNoMuPhi      );
  effTree->Branch("jet1Pt"             , &jet1Pt            );
  effTree->Branch("jet1Eta"            , &jet1Eta           );
  effTree->Branch("jet1Phi"            , &jet1Phi           );
  effTree->Branch("jet1M"              , &jet1Phi           );
  //effTree->Branch("jet1PtUp"           , &jet1PtUp          );
  //effTree->Branch("jet1PtDown"         , &jet1PtDown        );
  effTree->Branch("jet1ExtraClean"     , &jet1ExtraClean    );
  effTree->Branch("jet2Pt"             , &jet2Pt            );
  effTree->Branch("jet2Eta"            , &jet2Eta           );
  effTree->Branch("jet2Phi"            , &jet2Phi           );
  effTree->Branch("jet2M"              , &jet2Phi           );
  //effTree->Branch("jet2PtUp"           , &jet2PtUp          );
  //effTree->Branch("jet2PtDown"         , &jet2PtDown        );
  effTree->Branch("jet2ExtraClean"     , &jet2ExtraClean    );
  effTree->Branch("jet3Pt"             , &jet2Pt            );
  effTree->Branch("jet3Eta"            , &jet2Eta           );
  effTree->Branch("jet3Phi"            , &jet2Phi           );
  effTree->Branch("jet3M"              , &jet2Phi           );
  //effTree->Branch("jet3PtUp"           , &jet2PtUp          );
  //effTree->Branch("jet3PtDown"         , &jet2PtDown        );
  effTree->Branch("jet3ExtraClean"     , &jet2ExtraClean    );
  effTree->Branch("jet4Pt"             , &jet2Pt            );
  effTree->Branch("jet4Eta"            , &jet2Eta           );
  effTree->Branch("jet4Phi"            , &jet2Phi           );
  effTree->Branch("jet4M"              , &jet2Phi           );
  //effTree->Branch("jet4PtUp"           , &jet2PtUp          );
  //effTree->Branch("jet4PtDown"         , &jet2PtDown        );
  effTree->Branch("jet4ExtraClean"     , &jet2ExtraClean    );
  effTree->Branch("nTau"               , &nTau              );
  
  TChain input("events");
  if(inputPath!="") {
    //input.Add(inputPath);
    TFileCollection fc("fc","",inputPath);
    input.AddFileInfoList((TCollection*)fc.GetList());
  } else {
    if(mode=="electrons") {
      //input.Add("/mnt/hadoop/scratch/yiiyama/pandaf/006/SingleElectron+Run2017A-PromptReco-v2+MINIAOD/*.root");
      //input.Add("/mnt/hadoop/scratch/yiiyama/pandaf/006/SingleElectron+Run2017A-PromptReco-v3+MINIAOD/*.root");
      //input.Add("/mnt/hadoop/scratch/yiiyama/pandaf/006/SingleElectron+Run2017B-PromptReco-v1+MINIAOD/*.root");
      TFileCollection fc("fc","","PandaAnalysis/VBF/triggers/files_SingleElectron2017_RunAB_july18.txt");
      input.AddFileInfoList((TCollection*)fc.GetList());
    } else if(mode=="muons") {
      TFileCollection fc("fc","","PandaAnalysis/VBF/triggers/files_SingleMuon2017_RunAB_july18.txt");
      input.AddFileInfoList((TCollection*)fc.GetList());
    }
  }
  
  //read json file into boost property tree
  MapType fMap;
  boost::property_tree::ptree jsonTree;
  boost::property_tree::read_json(jsonFile.c_str(),jsonTree);
  
  //loop through boost property tree and fill the MapType structure with the list of good lumi
  //ranges for each run
  for (boost::property_tree::ptree::const_iterator it = jsonTree.begin(); it!=jsonTree.end(); ++it) {
    UInt_t runNumber = boost::lexical_cast<UInt_t>(it->first);
    MapType::mapped_type &lumiPairList = fMap[runNumber];
    boost::property_tree::ptree lumiPairListTree = it->second;
    for (boost::property_tree::ptree::const_iterator jt = lumiPairListTree.begin(); jt!=lumiPairListTree.end(); ++jt) {
      boost::property_tree::ptree lumiPairTree = jt->second;
      if (lumiPairTree.size()==2) {
        UInt_t firstLumi = boost::lexical_cast<UInt_t>(lumiPairTree.begin()->second.data());
        UInt_t lastLumi = boost::lexical_cast<UInt_t>((++lumiPairTree.begin())->second.data());
        lumiPairList.push_back(std::pair<UInt_t,UInt_t>(firstLumi,lastLumi));
      }
    }
  }

  //dump run and lumi ranges from MapType structure to verify correct json parsing
  if (debug) {
    printf("Iterating over parsed JSON:\n");
    for (MapType::const_iterator it = fMap.begin(); it != fMap.end(); ++it) {
      printf("  Run %u:\n",it->first);
      for (MapType::mapped_type::const_iterator jt = it->second.begin(); jt < it->second.end(); ++jt) {
        printf("    Lumis %u - %u\n",jt->first,jt->second);
      }
    }

  }
 
  panda::Event event;
  event.setStatus(input, {"!*"});
  event.setAddress(input, {"caloMet", "chsAK4Jets", "electrons", "muons", "pfMet", "photons", "taus", "triggers", "triggerObjects", "runNumber", "lumiNumber", "eventNumber", "metFilters"});
  long iEntry = 0;
  vector<unsigned> lepTriggerTokens, metTriggerTokens;
  if (mode=="electrons") {
    lepTriggerTokens.push_back( event.registerTrigger("HLT_Ele35_WPTight_Gsf") );
    if(type1) {
      metTriggerTokens.push_back( event.registerTrigger("HLT_PFMETTypeOne110_PFMHT110_IDTight") );
      metTriggerTokens.push_back( event.registerTrigger("HLT_PFMETTypeOne120_PFMHT120_IDTight") );
      metTriggerTokens.push_back( event.registerTrigger("HLT_PFMETTypeOne130_PFMHT130_IDTight") );
      metTriggerTokens.push_back( event.registerTrigger("HLT_PFMETTypeOne140_PFMHT140_IDTight") );
    } else {
      metTriggerTokens.push_back( event.registerTrigger("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight") );
      metTriggerTokens.push_back( event.registerTrigger("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight") );
    }
  } else if(mode=="muons") {
    lepTriggerTokens.push_back( event.registerTrigger("HLT_IsoMu24") );
    metTriggerTokens.push_back( event.registerTrigger("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight") );
    metTriggerTokens.push_back( event.registerTrigger("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight") );
  }
  
  while (event.getEntry(input, iEntry++) > 0) {
    // Reset branches
    passMetTriggers=false; passJetMetCleaning=false; jet1ExtraClean=false; jet2ExtraClean=false; jet3ExtraClean=false; jet4ExtraClean=false;
    leptonPt=-1; leptonEta=-1; leptonPhi=-1; pfMetNoMu=-1; pfMetNoMuPhi=-1; jet1Pt=-1; jet1Eta=-1; jet1Phi=-1; jet1PtUp=-1; jet1PtDown=-1; jet2Pt=-1; jet2Eta=-1; jet2Phi=-1; jet2PtUp=-1; jet2PtDown=-1;
    jet3Pt=-1; jet3Eta=-1; jet3Phi=-1; jet3PtUp=-1; jet3PtDown=-1; jet4Pt=-1; jet4Eta=-1; jet4Phi=-1; jet4PtUp=-1; jet4PtDown=-1;
    nTau=0;

    // Check data certification
    bool certifiedEvent=false;
    std::pair<unsigned int, unsigned int> runLumi(event.runNumber, event.lumiNumber);      
    MapType::const_iterator it = fMap.find(runLumi.first);
    if (it!=fMap.end()) {
      //check lumis
      const MapType::mapped_type &lumiPairList = it->second;
      for (MapType::mapped_type::const_iterator jt = lumiPairList.begin(); jt<lumiPairList.end(); ++jt) {
        if (runLumi.second >= jt->first && runLumi.second <= jt->second) {
          //found lumi in accepted range
          certifiedEvent=true;
        }
      }
    }
    if(!certifiedEvent) { if(debug) printf("failed golden json\n"); continue; }

    // Met quality requirements
    bool passMetFilters=(
      !event.metFilters.globalHalo16 &&
      !event.metFilters.hbhe &&
      !event.metFilters.hbheIso &&
      !event.metFilters.ecalDeadCell &&
      !event.metFilters.badsc &&
      !event.metFilters.goodVertices &&
      !event.metFilters.badPFMuons &&
      !event.metFilters.badChargedHadrons
    );
    if(!passMetFilters) { if(debug) printf("failed met filters\n"); continue; }
    bool passPfCaloBalance=(
      TMath::Abs(event.pfMet.pt - event.caloMet.pt) / event.pfMet.pt < 0.5
    );
    if(!passPfCaloBalance) { if (debug) printf("failed pf calo balance\n"); continue; }

    // Require lepton trigger
    bool passLepTriggers=false;
    for(unsigned i=0; i<lepTriggerTokens.size() && !passLepTriggers; i++) if(event.triggerFired(lepTriggerTokens[i])) passLepTriggers=true;
    if(!passLepTriggers) {  if(debug) printf("failed lepton triggers\n"); continue; }
    
    // Begin finding leptons
    std::vector<panda::Lepton*> looseLeps, tightLeps; // Fakeable object and tight leptons
    std::vector<panda::Particle*> matchLeps, matchPhos;
    std::vector<bool> tightLepIsTriggerMatched;

    // Find veto and tight electrons
    if(mode=="electrons") for (auto& ele : event.electrons) {
     float pt = ele.pt(); float eta = ele.eta(); float aeta = fabs(eta);
      if (pt<=10 || aeta>2.5 || (aeta>1.4442 && aeta<1.566)) // electron acceptance cuts
        continue;
      if (!ElectronIP(ele.eta(),ele.dxy,ele.dz)) continue;
      if (!ele.veto) continue;
      looseLeps.push_back(&ele); // passes MFO definition
      if(ele.tight) {
        tightLeps.push_back(&ele);
        tightLepIsTriggerMatched.push_back(ele.triggerMatch[panda::Electron::fEl35Tight]);
        matchLeps.push_back(&ele);
      }
    }
    
    // Find loose and tight muons
    TVector2 vMET, vMETNoMu;
    vMET.SetMagPhi(event.pfMet.pt, event.pfMet.phi);
    vMETNoMu = vMET; // initialize Met no Mu vector as just the PF MET
    if(mode=="muons") for (auto& mu : event.muons) {
      float pt = mu.pt(); float eta = mu.eta(); float aeta = fabs(eta);
      if (pt<=10 || aeta>2.4) continue; //muon acceptance cuts
      if (!mu.loose)  continue;
      if (!MuonIsolation(pt,eta,mu.combIso(),panda::kLoose)) continue;
      looseLeps.push_back(&mu); //passes MFO definition
      if(mu.tight && MuonIsolation(mu.pt(),mu.eta(),mu.combIso(),panda::kTight)) {
        tightLeps.push_back(&mu);
        tightLepIsTriggerMatched.push_back(mu.triggerMatch[panda::Muon::fIsoMu24]);
        matchLeps.push_back(&mu);
      }

      // add the muon pT vector to the Met no Mu
      TVector2 vMu; vMu.SetMagPhi(pt,mu.phi());
      vMETNoMu += vMu;
    }

    // Find hadronic taus
    for (auto& tau : event.taus) {
      if (!tau.decayMode || !tau.decayModeNew) continue;
      if (!tau.looseIsoMVAOld) continue;
      if (tau.pt()<18 || fabs(tau.eta())>2.3)  continue;
      if (IsMatched(&matchLeps,0.16,tau.eta(),tau.phi()))  continue;
      nTau++;
    }
    if(nTau>1) { if(debug) printf("failed tau veto"); continue; }
    // End finding leptons
    
    // Only consider events with exactly one lepton which is tight, apply extra lepton veto
    if(tightLeps.size()!=1 || looseLeps.size()!=tightLeps.size()) continue;

    // Require a trigger matched lepton
    bool matchedTriggerObject= (tightLepIsTriggerMatched[0]);
    if(!matchedTriggerObject) { if(debug) printf("failed trigger object matching\n"); continue; }
    
    // find photons (needed for jet cleaning)
    for (auto& pho : event.photons) {
      if (!pho.loose || !pho.csafeVeto)
        continue;
      float pt = pho.pt();
      if (pt<1) continue;
      float eta = pho.eta(), phi = pho.phi();
      if (pt<15 || fabs(eta)>2.5)
        continue;
      /*
      if (IsMatched(&matchEles,0.16,eta,phi))
        continue;
      */
      if ( pho.medium &&
        pt>175 /*&& fabs(eta)<1.4442*/ ) { // apply eta cut offline
        matchPhos.push_back(&pho);
      }
    }

    // Look for jets 
    panda::JetCollection* jets(0);
    jets = &event.chsAK4Jets;
    vector<panda::Jet*> centralAndForwardJets;
    TLorentzVector vJet;
    panda::Jet *jet1=0, *jet2=0, *jet3=0, *jet4=0, *jetUp1=0, *jetUp2=0, *jetDown1=0, *jetDown2=0; 
    unsigned char iJet=0; passJetMetCleaning=true;
    for (auto& jet : *jets) {
      // jet cleaning
      if (IsMatched(&matchLeps,0.16,jet.eta(),jet.phi()) ||
          IsMatched(&matchPhos,0.16,jet.eta(),jet.phi()))
        continue;
      vJet.SetPtEtaPhiM(jet.pt(),jet.eta(),jet.phi(),jet.m());
      if(iJet<4) passJetMetCleaning &= (vMET.DeltaPhi(vJet.Vect().XYvector()) > 0.5);
      if (jet.pt()>30) { // consider only central jets here
        if (fabs(jet.eta())<4.7) {
          centralAndForwardJets.push_back(&jet);
          if (centralAndForwardJets.size()==1) {
            jet1 = &jet;
            jet1Pt = jet.pt();
            jet1Eta = jet.eta();
            jet1Phi = jet.phi();
            jet1M = jet.m();
            jet1ExtraClean = jet.monojet;
          } else if (centralAndForwardJets.size()==2) {
            jet2 = &jet;
            jet2Pt = jet.pt();
            jet2Eta = jet.eta();
            jet2Phi = jet.phi();
            jet2M = jet.m();
            jet2ExtraClean = jet.monojet;
          } else if (centralAndForwardJets.size()==3) {
            jet3 = &jet;
            jet3Pt = jet.pt();
            jet3Eta = jet.eta();
            jet3Phi = jet.phi();
            jet3M = jet.m();
            jet3ExtraClean = jet.monojet;
          } else if (centralAndForwardJets.size()==4) {
            jet4 = &jet;
            jet4Pt = jet.pt();
            jet4Eta = jet.eta();
            jet4Phi = jet.phi();
            jet4M = jet.m();
            jet4ExtraClean = jet.monojet;
          }
        }
      }
      if(jet1Pt<80.) continue;

      /*
      // do jes variation OUTSIDE of pt>30 check
      if (jet.ptCorrUp>30) {
        if (jet.ptCorrUp > jet1PtUp) {
          if (jetUp1) {
            jetUp2 = jetUp1;
            jet2PtUp  = jet1PtUp;
          }
          jetUp1 = &jet;
          jet1PtUp = jet.ptCorrUp;
        } else if (jet.ptCorrUp > jet2PtUp) {
          jetUp2 = &jet;
          jet2PtUp = jet.ptCorrUp;
        }
      }
      if (jet.ptCorrDown>30) {
        if (jet.ptCorrDown > jet1PtDown) {
          if (jetDown1) {
            jetDown2 = jetDown1;
            jet2PtDown  = jet1PtDown;
          }
          jetDown1 = &jet;
          jet1PtDown = jet.ptCorrDown;
        } else if (jet.ptCorrDown > jet2PtDown) {
          jetDown2 = &jet;
          jet2PtDown = jet.ptCorrDown;
        }
      }*/
      iJet++;
    } //end loop over jets
    
    // Record if the event pass the MET triggers
    passMetTriggers=false;
    for(unsigned i=0; i<metTriggerTokens.size() && !passMetTriggers; i++) if(event.triggerFired(metTriggerTokens[i])) passMetTriggers=true;
    
    leptonPt=tightLeps[0]->pt();
    leptonEta=tightLeps[0]->eta();
    leptonPhi=tightLeps[0]->phi();
    pfMetNoMu=vMETNoMu.Mod();
    pfMetNoMuPhi=vMETNoMu.Phi();
    
    runNumber = event.runNumber;
    lumiSection = event.lumiNumber;
    eventNumber = event.eventNumber;
    effTree->Fill();  
    //printf("run %d ls %d evt %lld %s\n", event.runNumber, event.lumiNumber, event.eventNumber, ((TFile*)input.GetCurrentFile())->GetName()); 
    printf("run %d ls %d evt %u %s\n", event.runNumber, event.lumiNumber, event.eventNumber, ((TFile*)input.GetCurrentFile())->GetName()); 
    
  }
  outputFile->cd();
  effTree->Write();
  outputFile->Close();
}
