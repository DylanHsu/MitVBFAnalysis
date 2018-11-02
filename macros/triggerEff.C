#include "MitVBFAnalysis/macros/triggerEff.h"

// Last modified by DGH July 18 2017
// This code produces a flat tree for trigger studies from Panda files (provided as inputPath var)
// Compiles through ACLiC and runs on condor via the shell scripts
// Analyze the output using plotting scripts in  MitVBFAnalysis/plotting

typedef std::map<UInt_t,std::vector<std::pair <UInt_t, UInt_t> > > MapType;
string jsonFile = "MitVBFAnalysis/data/Cert_294927-300575_13TeV_PromptReco_Collisions17_JSON.txt";

void triggerEff(
TString mode="metElectronRef",
TString inputPath="",
TString batchName="",
bool type1=false,
bool debug=false,
TString customFilter=""
) {
  assert(mode=="metElectronRef" || mode=="metMuonRef" || mode=="electronTnP" || mode=="photonPFHTRef" || mode=="photonL1Ref");

  UInt_t runNumber, lumiSection;
  ULong64_t eventNumber;
  //Read json file into boost property tree
  MapType fMap;
  boost::property_tree::ptree jsonTree;
  boost::property_tree::read_json(jsonFile.c_str(),jsonTree);
  
  //Loop through boost property tree and fill the MapType structure with the list of good lumi
  //ranges for each run
  for (boost::property_tree::ptree::const_iterator it = jsonTree.begin(); it!=jsonTree.end(); ++it) {
    runNumber = boost::lexical_cast<UInt_t>(it->first);
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
  //If running in debug mode, dump run and lumi ranges from MapType structure to verify correct json parsing
  if (debug) {  printf("Iterating over parsed JSON:\n"); for (MapType::const_iterator it = fMap.begin(); it != fMap.end(); ++it) {
    printf("  Run %u:\n",it->first);
    for (MapType::mapped_type::const_iterator jt = it->second.begin(); jt < it->second.end(); ++jt) printf("    Lumis %u - %u\n",jt->first,jt->second);
  }}
  
  TFile *outputFile = TFile::Open(("vbf_batchTree_triggerEff"+(batchName == "" ? "" : "_"+batchName)+".root"), "RECREATE", "", ROOT::CompressionSettings(ROOT::kZLIB,9));
  // Declare the Tree and set up the branches.
  TTree *effTree = new TTree("effTree", "effTree");
  effTree->Branch("runNumber"          , &runNumber         );
  effTree->Branch("lumiSection"        , &lumiSection       );
  effTree->Branch("eventNumber"        , &eventNumber       );
  // Declare the simple variables for the tree branches
  bool passMetTriggers, passJetMetCleaning, jet1ExtraClean, jet2ExtraClean, jet3ExtraClean, jet4ExtraClean;
  bool passElectronTrigger, passPhotonTriggers;
  float lepton1Pt, lepton1Eta, lepton1Phi, lepton2Pt, lepton2Eta, lepton2Phi;
  float photon1Pt, photon1Eta, photon1Phi;
  float pfMet, pfMetPhi, pfMetNoMu, pfMetNoMuPhi;
  float jet1Pt, jet1Eta, jet1Phi, jet1M, jet1PtUp, jet1PtDown, jet2Pt, jet2Eta, jet2Phi, jet2M, jet2PtUp, jet2PtDown, jet3Pt, jet3Eta, jet3Phi, jet3M, jet3PtUp, jet3PtDown, jet4Pt, jet4Eta, jet4Phi, jet4M, jet4PtUp, jet4PtDown;
  unsigned nTau;
  
  // Create maps (filter name) -> (trigger object coordinate)
  // The filter names are currently hard coded, need to store this in a data file
  // The addresses to the values will be passed to effTree branches
  // For now, not storing Eta coordinates, and the filters we want Phi information for must also be in the Pt map
  std::map< TString, float > triggerObjectPts, triggerObjectEtas, triggerObjectPhis;
  
  // Initialize mode-specific maps and branches
  if(mode=="metElectronRef"||mode=="metMuonRef") {
    // trigger object map for met triggers
    triggerObjectPts["hltMET80"]=0;
    triggerObjectPts["hltMET90"]=0;
    triggerObjectPts["hltMETClean70"]=0;
    triggerObjectPts["hltMETClean80"]=0;
    triggerObjectPts["hltMHT80"]=0;
    triggerObjectPts["hltMHT90"]=0;
    triggerObjectPts["hltPFMHTNoMuTightID110"]=0;
    triggerObjectPts["hltPFMHTNoMuTightID120"]=0;
    triggerObjectPts["hltPFMETNoMu110"]=0;
    triggerObjectPts["hltPFMETNoMu120"]=0;
    triggerObjectPhis["hltMETClean70"]=0;
    triggerObjectPhis["hltMETClean80"]=0;
    triggerObjectPhis["hltMHT80"]=0;
    triggerObjectPhis["hltMHT90"]=0;
    // specific tree branches
    effTree->Branch("passMetTriggers"    , &passMetTriggers   );
    effTree->Branch("passJetMetCleaning" , &passJetMetCleaning);
    effTree->Branch("lepton1Pt"           , &lepton1Pt          );
    effTree->Branch("lepton1Eta"          , &lepton1Eta         );
    effTree->Branch("lepton1Phi"          , &lepton1Phi         );
    effTree->Branch("pfMet"              , &pfMet             );
    effTree->Branch("pfMetPhi"           , &pfMetPhi          );
    effTree->Branch("pfMetNoMu"          , &pfMetNoMu         );
    effTree->Branch("pfMetNoMuPhi"       , &pfMetNoMuPhi      );
    effTree->Branch("jet1Pt"             , &jet1Pt            );
    effTree->Branch("jet1Eta"            , &jet1Eta           );
    effTree->Branch("jet1Phi"            , &jet1Phi           );
    effTree->Branch("jet1M"              , &jet1Phi           );
    effTree->Branch("jet1ExtraClean"     , &jet1ExtraClean    );
    effTree->Branch("jet2Pt"             , &jet2Pt            );
    effTree->Branch("jet2Eta"            , &jet2Eta           );
    effTree->Branch("jet2Phi"            , &jet2Phi           );
    effTree->Branch("jet2M"              , &jet2Phi           );
    effTree->Branch("jet2ExtraClean"     , &jet2ExtraClean    );
    effTree->Branch("jet3Pt"             , &jet2Pt            );
    effTree->Branch("jet3Eta"            , &jet2Eta           );
    effTree->Branch("jet3Phi"            , &jet2Phi           );
    effTree->Branch("jet3M"              , &jet2Phi           );
    effTree->Branch("jet3ExtraClean"     , &jet2ExtraClean    );
    effTree->Branch("jet4Pt"             , &jet2Pt            );
    effTree->Branch("jet4Eta"            , &jet2Eta           );
    effTree->Branch("jet4Phi"            , &jet2Phi           );
    effTree->Branch("jet4M"              , &jet2Phi           );
    effTree->Branch("jet4ExtraClean"     , &jet2ExtraClean    );
    effTree->Branch("nTau"               , &nTau              );
    effTree->Branch("hltMET80",                &triggerObjectPts["hltMET80"]               );
    effTree->Branch("hltMET90",                &triggerObjectPts["hltMET90"]               );
    effTree->Branch("hltMETClean70",           &triggerObjectPts["hltMETClean70"]          );
    effTree->Branch("hltMETClean80",           &triggerObjectPts["hltMETClean80"]          );
    effTree->Branch("hltMHT80",                &triggerObjectPts["hltMHT80"]               );
    effTree->Branch("hltMHT90",                &triggerObjectPts["hltMHT90"]               );
    effTree->Branch("hltPFMHTNoMuTightID110",  &triggerObjectPts["hltPFMHTNoMuTightID110"] );
    effTree->Branch("hltPFMHTNoMuTightID120",  &triggerObjectPts["hltPFMHTNoMuTightID120"] );
    effTree->Branch("hltPFMETNoMu110",         &triggerObjectPts["hltPFMETNoMu110"]        );
    effTree->Branch("hltPFMETNoMu120",         &triggerObjectPts["hltPFMETNoMu120"]        );
    effTree->Branch("hltMETClean70Phi",        &triggerObjectPhis["hltMETClean70"]      );
    effTree->Branch("hltMETClean80Phi",        &triggerObjectPhis["hltMETClean80"]      );
    effTree->Branch("hltMHT80Phi",             &triggerObjectPhis["hltMHT80"]           );
    effTree->Branch("hltMHT90Phi",             &triggerObjectPhis["hltMHT90"]           );
  } else if(mode=="electronTnP") {
    effTree->Branch("passElectronTrigger", &passElectronTrigger);
    effTree->Branch("lepton1Pt"          , &lepton1Pt         );
    effTree->Branch("lepton1Eta"         , &lepton1Eta        );
    effTree->Branch("lepton1Phi"         , &lepton1Phi        );
    effTree->Branch("lepton2Pt"          , &lepton2Pt         );
    effTree->Branch("lepton2Eta"         , &lepton2Eta        );
    effTree->Branch("lepton2Phi"         , &lepton2Phi        );
    effTree->Branch("pfMet"              , &pfMet             );
    effTree->Branch("pfMetPhi"           , &pfMetPhi          );
    effTree->Branch("pfMetNoMu"          , &pfMetNoMu         );
    effTree->Branch("pfMetNoMuPhi"       , &pfMetNoMuPhi      );
    effTree->Branch("jet1Pt"             , &jet1Pt            );
    effTree->Branch("jet1Eta"            , &jet1Eta           );
    effTree->Branch("jet1Phi"            , &jet1Phi           );
    effTree->Branch("jet1M"              , &jet1Phi           );
    effTree->Branch("jet1ExtraClean"     , &jet1ExtraClean    );
    effTree->Branch("jet2Pt"             , &jet2Pt            );
    effTree->Branch("jet2Eta"            , &jet2Eta           );
    effTree->Branch("jet2Phi"            , &jet2Phi           );
    effTree->Branch("jet2M"              , &jet2Phi           );
    effTree->Branch("jet2ExtraClean"     , &jet2ExtraClean    );
    effTree->Branch("nTau"               , &nTau              );
  } else if(mode=="photonPFHTRef") {
    effTree->Branch("passPhotonTriggers" , &passPhotonTriggers);
    effTree->Branch("photon1Pt"          , &photon1Pt         );
    effTree->Branch("photon1Eta"         , &photon1Eta        );
    effTree->Branch("photon1Phi"         , &photon1Phi        );
    effTree->Branch("pfMet"              , &pfMet             );
    effTree->Branch("pfMetPhi"           , &pfMetPhi          );
    effTree->Branch("pfMetNoMu"          , &pfMetNoMu         );
    effTree->Branch("pfMetNoMuPhi"       , &pfMetNoMuPhi      );
    effTree->Branch("jet1Pt"             , &jet1Pt            );
    effTree->Branch("jet1Eta"            , &jet1Eta           );
    effTree->Branch("jet1Phi"            , &jet1Phi           );
    effTree->Branch("jet1M"              , &jet1Phi           );
    effTree->Branch("jet1ExtraClean"     , &jet1ExtraClean    );
    effTree->Branch("jet2Pt"             , &jet2Pt            );
    effTree->Branch("jet2Eta"            , &jet2Eta           );
    effTree->Branch("jet2Phi"            , &jet2Phi           );
    effTree->Branch("jet2M"              , &jet2Phi           );
    effTree->Branch("jet2ExtraClean"     , &jet2ExtraClean    );
    effTree->Branch("nTau"               , &nTau              );
  }
  if(debug) printf("Initialized the tree successfully\n");
  TChain input("events");
  if(inputPath!="") {
    // Use a TChain to read one or several files as listed in the file located at inputPath
    TFileCollection fc("fc","",inputPath);
    input.AddFileInfoList((TCollection*)fc.GetList());
  } else {
    // If no inputPath is provided, loop over all the files for given dataset
    if(mode=="metElectronRef") {
      TFileCollection fc("fc","","MitVBFAnalysis/data/files_SingleElectron2017_RunAB_july18.txt");
      input.AddFileInfoList((TCollection*)fc.GetList());
    } else if(mode=="metMuonRef") {
      TFileCollection fc("fc","","MitVBFAnalysis/data/files_SingleMuon2017_RunAB_july18.txt");
      input.AddFileInfoList((TCollection*)fc.GetList());
    } else return;
  }
  
 
  // Initialize the Panda event object
  panda::Event event;
  event.setStatus(input, {"!*"});
  event.setAddress(input, {"caloMet", "chsAK4Jets", "electrons", "muons", "pfMet", "photons", "taus", "triggers", "triggerObjects", "runNumber", "lumiNumber", "eventNumber", "metFilters"});
  
  // Register the trigger tokens for the reference and test triggers
  vector<unsigned> refTriggerTokens, testTriggerTokens;
  if (mode=="metElectronRef") {
    refTriggerTokens.push_back( event.registerTrigger("HLT_Ele35_WPTight_Gsf") );
    if(type1) {
      testTriggerTokens.push_back( event.registerTrigger("HLT_PFMETTypeOne110_PFMHT110_IDTight") );
      testTriggerTokens.push_back( event.registerTrigger("HLT_PFMETTypeOne120_PFMHT120_IDTight") );
      testTriggerTokens.push_back( event.registerTrigger("HLT_PFMETTypeOne130_PFMHT130_IDTight") );
      testTriggerTokens.push_back( event.registerTrigger("HLT_PFMETTypeOne140_PFMHT140_IDTight") );
    } else {
      testTriggerTokens.push_back( event.registerTrigger("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight") );
      testTriggerTokens.push_back( event.registerTrigger("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight") );
      //testTriggerTokens.push_back( event.registerTrigger("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight") );
      //testTriggerTokens.push_back( event.registerTrigger("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight") );
    }
  } else if(mode=="metMuonRef") {
    refTriggerTokens.push_back( event.registerTrigger("HLT_IsoMu24") );
    testTriggerTokens.push_back( event.registerTrigger("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight") );
    testTriggerTokens.push_back( event.registerTrigger("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight") );
    //testTriggerTokens.push_back( event.registerTrigger("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight") );
    //testTriggerTokens.push_back( event.registerTrigger("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight") );
  } else if(mode=="electronTnP") {
    refTriggerTokens.push_back( event.registerTrigger("HLT_Ele35_WPTight_Gsf") );
    testTriggerTokens.push_back( event.registerTrigger("HLT_Ele115_CaloIdVT_GsfTrkIdT") );
  } else if(mode=="photonPFHTRef") {
    refTriggerTokens.push_back( event.registerTrigger("HLT_PFHT590"));
    refTriggerTokens.push_back( event.registerTrigger("HLT_PFHT680"));
    refTriggerTokens.push_back( event.registerTrigger("HLT_PFHT780"));
    testTriggerTokens.push_back( event.registerTrigger("HLT_Photon200"));
  }
  
  // Loop over all of the events in the TChain
  long iEntry = 0;
  while (event.getEntry(input, iEntry++) > 0) {
    // Reset branches
    passMetTriggers=false; passJetMetCleaning=false; jet1ExtraClean=false; jet2ExtraClean=false; jet3ExtraClean=false; jet4ExtraClean=false;
    lepton1Pt=-1; lepton1Eta=-1; lepton1Phi=-1; pfMetNoMu=-1; pfMetNoMuPhi=-1; jet1Pt=-1; jet1Eta=-1; jet1Phi=-1; jet1PtUp=-1; jet1PtDown=-1; jet2Pt=-1; jet2Eta=-1; jet2Phi=-1; jet2PtUp=-1; jet2PtDown=-1;
    jet3Pt=-1; jet3Eta=-1; jet3Phi=-1; jet3PtUp=-1; jet3PtDown=-1; jet4Pt=-1; jet4Eta=-1; jet4Phi=-1; jet4PtUp=-1; jet4PtDown=-1;
    nTau=0;

    if(debug) printf("Run %u, LS %u, evt %lld:\n", event.runNumber, event.lumiNumber, event.eventNumber);
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
    if((mode=="metElectronRef"||mode=="metMuonRef") && !passMetFilters) { if(debug) printf("failed met filters\n"); continue; }
    bool passPfCaloBalance=(
      TMath::Abs(event.pfMet.pt - event.caloMet.pt) / event.pfMet.pt < 0.5
    );
    if((mode=="metElectronRef"||mode=="metMuonRef") && !passPfCaloBalance) { if (debug) printf("failed pf calo balance\n"); continue; }

    // Require reference trigger
    bool passRefTriggers=false;
    for(unsigned i=0; i<refTriggerTokens.size() && !passRefTriggers; i++) if(event.triggerFired(refTriggerTokens[i])) passRefTriggers=true;
    if(!passRefTriggers) {  if(debug) printf("failed reference triggers\n"); continue; }

    // Begin finding leptons
    std::vector<panda::Lepton*> looseLeps, tightLeps; // Fakeable object and tight leptons
    std::vector<panda::Particle*> matchLeps, matchPhos;
    std::vector<panda::Photon*> tightPhos;
    std::vector<bool> tightLepIsTriggerMatched, tightPhoIsTriggerMatched;

    // Find veto and tight electrons
    if(mode=="metElectronRef"||mode=="electronTnP") for (auto& ele : event.electrons) {
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
    if(mode=="metMuonRef") for (auto& mu : event.muons) {
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
    if((mode=="metElectronRef"||mode=="metMuonRef") && nTau>0) { if(debug) printf("failed tau veto"); continue; }
    // End finding leptons
    
    // Only consider events with exactly one lepton which is tight, apply extra lepton veto
    if((mode=="metElectronRef" || mode=="metMuonRef") && (tightLeps.size()!=1 || looseLeps.size()!=tightLeps.size())) continue;
    if(mode=="electronTnP"  && (tightLeps.size()!=2 || looseLeps.size()!=tightLeps.size())) continue;

    // Require a trigger matched lepton
    bool matchedTriggerObject=false;
    if(mode=="metElectronRef"||mode=="metMuonRef") matchedTriggerObject = (tightLepIsTriggerMatched[0]);
    if(mode=="electronTnP") matchedTriggerObject = (tightLepIsTriggerMatched[0] || tightLepIsTriggerMatched[1]);
    if(mode=="photonPFHTRef") matchedTriggerObject=true;
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
      if(mode=="photonPFHTRef" && pho.medium && fabs(eta)<1.4442) {
        tightPhos.push_back(&pho);
        tightPhoIsTriggerMatched.push_back(pho.triggerMatch[panda::Photon::fPh200]);
      }
    }

    if(mode=="photonPFHTRef" && tightPhos.size()==0) { if(debug) printf("no suitable photons\n"); continue; }
    
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
    
    pfMetNoMu=vMETNoMu.Mod();
    pfMetNoMuPhi=vMETNoMu.Phi();
    pfMet=event.pfMet.pt;
    pfMetPhi=event.pfMet.phi;
    runNumber = event.runNumber;
    lumiSection = event.lumiNumber;
    eventNumber = event.eventNumber;
    
    if(mode=="metElectronRef" || mode=="metMuonRef") {
      // July 17 "Update V2 to GRun/V86 HLT/V293 "  CMSSW_9_3_X
      // https://github.com/cms-sw/cmssw/commit/eed72a09f9fc4a775f959d7a556abc6fa017a6d7
      // fragment.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v12 = cms.Path( fragment.HLTBeginSequence + fragment.hltL1sAllETMHadSeeds + fragment.hltPrePFMETNoMu110PFMHTNoMu110IDTight + fragment.HLTRecoMETSequence + fragment.hltMET80 + fragment.HLTHBHENoiseCleanerSequence + fragment.hltMetClean + fragment.hltMETClean70 + fragment.HLTAK4CaloJetsSequence + fragment.hltMht + fragment.hltMHT80 + fragment.HLTAK4PFJetsSequence + fragment.hltPFMHTNoMuTightID + fragment.hltPFMHTNoMuTightID110 + fragment.hltParticleFlowNoMu + fragment.hltPFMETNoMuProducer + fragment.hltPFMETNoMu110 + fragment.HLTEndSequence )
      // fragment.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v12 = cms.Path( fragment.HLTBeginSequence + fragment.hltL1sAllETMHadSeeds + fragment.hltPrePFMETNoMu120PFMHTNoMu120IDTight + fragment.HLTRecoMETSequence + fragment.hltMET90 + fragment.HLTHBHENoiseCleanerSequence + fragment.hltMetClean + fragment.hltMETClean80 + fragment.HLTAK4CaloJetsSequence + fragment.hltMht + fragment.hltMHT90 + fragment.HLTAK4PFJetsSequence + fragment.hltPFMHTNoMuTightID + fragment.hltPFMHTNoMuTightID120 + fragment.hltParticleFlowNoMu + fragment.hltPFMETNoMuProducer + fragment.hltPFMETNoMu120 + fragment.HLTEndSequence )
      // fragment.HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v11 = cms.Path( fragment.HLTBeginSequence + fragment.hltL1sAllETMHadSeeds + fragment.hltPrePFMETNoMu130PFMHTNoMu130IDTight + fragment.HLTRecoMETSequence + fragment.hltMET100 + fragment.HLTHBHENoiseCleanerSequence + fragment.hltMetClean + fragment.hltMETClean90 + fragment.HLTAK4CaloJetsSequence + fragment.hltMht + fragment.hltMHT100 + fragment.HLTAK4PFJetsSequence + fragment.hltPFMHTNoMuTightID + fragment.hltPFMHTNoMuTightID130 + fragment.hltParticleFlowNoMu + fragment.hltPFMETNoMuProducer + fragment.hltPFMETNoMu130 + fragment.HLTEndSequence )
      // fragment.HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v11 = cms.Path( fragment.HLTBeginSequence + fragment.hltL1sAllETMHadSeeds + fragment.hltPrePFMETNoMu140PFMHTNoMu140IDTight + fragment.HLTRecoMETSequence + fragment.hltMET110 + fragment.HLTHBHENoiseCleanerSequence + fragment.hltMetClean + fragment.hltMETClean100 + fragment.HLTAK4CaloJetsSequence + fragment.hltMht + fragment.hltMHT110 + fragment.HLTAK4PFJetsSequence + fragment.hltPFMHTNoMuTightID + fragment.hltPFMHTNoMuTightID140 + fragment.hltParticleFlowNoMu + fragment.hltPFMETNoMuProducer + fragment.hltPFMETNoMu140 + fragment.HLTEndSequence )

      // Loop over the map of trigger object Pts and store them
      // Store Phi as well for the ones we care about
      // Currently getting only the First trigger object for MET-like quantities,
      // but this logic should be reworked if you use it for a filter with more than 1 trigger object
      if(debug) printf("Starting the filter objects loop\n");
      for( std::map<TString,float>::iterator iter = triggerObjectPts.begin(); iter != triggerObjectPts.end(); ++iter) {
        TString theFilterName = iter->first; // Get the filter name from the first half of the (Key)->(Value) pair
        if(debug) printf("Getting trigger objects for filter %s\n", theFilterName.Data());
        panda::HLTObjectStore::HLTObjectVector theTriggerObjects = event.triggerObjects.filterObjects(theFilterName.Data()); // filterObjects argument is const char*
        if(debug) printf("setting triggerObjectPts[%s] to -1\n", theFilterName.Data());
        triggerObjectPts[theFilterName] = -1.; // Set the dummy value
        bool recordPhi = (triggerObjectPhis.find(theFilterName) != triggerObjectPhis.end()); // Check if we are recording the phi for this filter
        if(recordPhi && debug) printf("setting triggerObjectPhis[%s] to -1\n", theFilterName.Data());
        if(recordPhi) triggerObjectPhis[theFilterName]=-777.; // Set the dummy value
        if(theTriggerObjects.size() == 0) continue; // If no trigger objects for this filter, go onto the nextone
        
        if(debug) printf("setting triggerObjectPts[%s] to %f\n", theFilterName.Data(), theTriggerObjects[0]->pt());
        triggerObjectPts[theFilterName]=theTriggerObjects[0]->pt();
        if(recordPhi && debug) printf("setting triggerObjectPhis[%s] to %f\n", theFilterName.Data(), theTriggerObjects[0]->phi());
        if(recordPhi) triggerObjectPhis[theFilterName]=theTriggerObjects[0]->phi();
        if(theTriggerObjects.size() > 1) printf("WARNING: got multiple (%zu) trigger objects for filter %s\n", theTriggerObjects.size(), theFilterName.Data());
      }

      // Record if the event pass the MET triggers
      passMetTriggers=false;
      for(unsigned i=0; i<testTriggerTokens.size() && !passMetTriggers; i++) if(event.triggerFired(testTriggerTokens[i])) passMetTriggers=true;
      lepton1Pt=tightLeps[0]->pt();
      lepton1Eta=tightLeps[0]->eta();
      lepton1Phi=tightLeps[0]->phi();
      lepton2Pt=-1;
      lepton2Eta=-10;
      lepton2Phi=-10;
      effTree->Fill();  
      if(debug) printf("recorded event: run %d ls %d evt %u in file %s\n", event.runNumber, event.lumiNumber, event.eventNumber, ((TFile*)input.GetCurrentFile())->GetName()); 
    } else if(mode=="electronTnP") {
      // Try to trigger match at least one electron
      /*
      bool isTriggerMatched=false;
      passElectronTrigger=false; // test leg passes trigger
      panda::HLTObjectStore::HLTObjectVector theTriggerObjects = event.triggerObjects.filterObjects(theFilterName.Data());
      if(theTriggerObjects.size()==0) continue;
      TLorentzVector triggerObject, recoElectron;
      // NEED TO CHECK THIS LOGIC
      for(unsigned i=0; i<theTriggerObjects.size(); i++) {
        triggerObject.SetPtEtaPhiM(theTriggerObjects[i]->pt(), theTriggerObjects[i]->eta(), theTriggerObjects[i]->phi(), theTriggerObjects[i]->m());
        for(unsigned j=0; j<tightLeps.size(); j++) {
          recoElectron.SetPtEtaPhiM(tightLeps[j]->pt(), tightLeps[j]->eta(), tightLeps[j]->phi(), 0.000511);
          double dr=recoElectron.DeltaR(triggerObject);
          if(isTriggerMatched && dr<0.1) passElectronTrigger=true;
          else if(dr<0.1) isTriggerMatched=true;
        }
      }*/
      //TString theFilterName="hltEle35noerWPTightGsfTrackIsoFilter"; // last filter for Ele35, hardcoded
      TString theFilterName="hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter";
      bool useFilterName=true;
      //if(customFilter!="") theFilterName=customFilter; useFilterName=true;
      panda::HLTObjectStore::HLTObjectVector theTriggerObjects;
      bool triggerFired;
      if(useFilterName) {
        triggerFired=false;
        // check if the trigger fired at all
        for(unsigned i=0; i<testTriggerTokens.size() && !triggerFired; i++) {
          if(event.triggerFired(testTriggerTokens[i])) triggerFired=true;
        }
        if(!triggerFired) { if(debug) printf("test triggers did not fire\n");
        } else {
          theTriggerObjects=event.triggerObjects.filterObjects(theFilterName.Data());
          if(theTriggerObjects.size()==0) { if(debug) printf("no trigger objects to match to for filter %s\n", theFilterName.Data()); }
          else { if(debug) printf("found %d trigger objects for filter %s\n", theTriggerObjects.size(), theFilterName.Data() ); }
        }
      }
      for(unsigned i=0; i<tightLeps.size(); i++) { // tag loop
        if(!tightLepIsTriggerMatched[i]) continue;
        for(unsigned j=0; j<tightLeps.size(); j++) { // probe loop
          passElectronTrigger=false;
          if(j==i) continue;
          if(!useFilterName) passElectronTrigger=tightLepIsTriggerMatched[j];
          
          else if(theTriggerObjects.size()>0 && triggerFired) {
            TLorentzVector triggerObject, recoElectron;
            recoElectron.SetPtEtaPhiM(tightLeps[j]->pt(), tightLeps[j]->eta(), tightLeps[j]->phi(), 0.000511); // create 4-vector for probe lepton
            // NEED TO CHECK THIS LOGIC
            for(unsigned k=0; k<theTriggerObjects.size(); k++) {
              // create 4-vector for trigger object
              triggerObject.SetPtEtaPhiM(theTriggerObjects[k]->pt(), theTriggerObjects[k]->eta(), theTriggerObjects[k]->phi(), theTriggerObjects[k]->m());
              double dr=recoElectron.DeltaR(triggerObject);
              if(dr<0.1) passElectronTrigger=true;
            }
          }
          lepton1Pt  = tightLeps[i]->pt();
          lepton1Eta = tightLeps[i]->eta();
          lepton1Phi = tightLeps[i]->phi();
          lepton2Pt  = tightLeps[j]->pt();
          lepton2Eta = tightLeps[j]->eta();
          lepton2Phi = tightLeps[j]->phi();
          effTree->Fill();
          if(debug) {
            printf("recorded event: run %d ls %d evt %u in file %s\n", event.runNumber, event.lumiNumber, event.eventNumber, ((TFile*)input.GetCurrentFile())->GetName()); 
            printf("tag pt/eta/phi: (%f,%f,%f) probe pt/eta/phi: (%f,%f,%f)\n", tightLeps[i]->pt(), tightLeps[i]->eta(), tightLeps[i]->phi(), tightLeps[j]->pt(), tightLeps[j]->eta(), tightLeps[j]->phi());
          }
        }
      }
    } else if(mode=="photonPFHTRef") {
      passPhotonTriggers=false;
      for(unsigned i=0; i<testTriggerTokens.size() && !passPhotonTriggers; i++) if(event.triggerFired(testTriggerTokens[i])) passPhotonTriggers=true;
      photon1Pt = tightPhos[0]->pt();
      photon1Eta = tightPhos[0]->eta();
      photon1Phi = tightPhos[0]->phi();
      effTree->Fill();
      if(debug) printf("recorded event: run %d ls %d evt %u in file %s\n", event.runNumber, event.lumiNumber, event.eventNumber, ((TFile*)input.GetCurrentFile())->GetName()); 

    }
    //printf("run %d ls %d evt %lld %s\n", event.runNumber, event.lumiNumber, event.eventNumber, ((TFile*)input.GetCurrentFile())->GetName()); 
    
  }
  outputFile->cd();
  effTree->Write();
  outputFile->Close();
}
