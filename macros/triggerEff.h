#include <map>

// Root libraries
#include <Compression.h>
#include <TChain.h>
#include <TFile.h>
#include <TFileCollection.h>
#include <TString.h>
#include <TTree.h>

// Panda
#include "PandaCore/Tools/interface/Common.h"
#include "PandaCore/Tools/interface/DataTools.h"
#include "PandaCore/Tools/interface/JERReader.h"
#include "PandaTree/Objects/interface/Event.h"

//Boost
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/lexical_cast.hpp>


////////////////////////////////////////////////////////////////////////////////////

inline double TTNLOToNNLO(double pt) {
    double a = 0.1102;
    double b = 0.1566;
    double c = -3.685e-4;
    double d = 1.098;

    return TMath::Min(1.25,
                        a*TMath::Exp(-b*pow(pt,2)+1) + c*pt + d);
}

////////////////////////////////////////////////////////////////////////////////////

class LumiRange {
public:
    LumiRange(int l0_,int l1_):
        l0(l0_),
        l1(l1_)
     { }
    ~LumiRange() {}
    bool Contains(int l) {
        return l0<=l && l<=l1;
    }
private:
    int l0, l1;
};


////////////////////////////////////////////////////////////////////////////////////
template <typename T>
class THCorr {
public:
    // wrapper around TH* to do corrections
    THCorr(T *h_) {
        h = h_;
        dim = h->GetDimension();
        TAxis *thurn = h->GetXaxis(); 
        lo1 = thurn->GetBinCenter(1);
        hi1 = thurn->GetBinCenter(thurn->GetNbins());
        if (dim>1) {
            TAxis *taxis = h->GetYaxis();
            lo2 = taxis->GetBinCenter(1);
            hi2 = taxis->GetBinCenter(taxis->GetNbins());
        }
    }
    ~THCorr() {} // does not own histogram!
    double Eval(double x) {
        if (dim!=1) {
          PError("THCorr1::Eval",
              TString::Format("Trying to access a non-1D histogram (%s)!",h->GetName()));
          return -1;
        }
        return getVal(h,bound(x,lo1,hi1));
    }

    double Eval(double x, double y) {
        if (dim!=2) {
          PError("THCorr1::Eval",
             TString::Format("Trying to access a non-2D histogram (%s)!",h->GetName()));
          return -1;
        }
        return getVal(h,bound(x,lo1,hi1),bound(y,lo2,hi2));
    }

    T *GetHist() { return h; }

private:
    T *h;
    int dim;
    double lo1, lo2, hi1, hi2;
};

typedef THCorr<TH1D> THCorr1;
typedef THCorr<TH2D> THCorr2;

////////////////////////////////////////////////////////////////////////////////////

namespace panda {
  enum IDWorkingPoint {
    kVeto,
    kLoose,
    kMedium,
    kTight,
    nIDWorkingPoints
  };
}

inline bool MuonIsolation(double pt, double eta, double iso, panda::IDWorkingPoint isoType) {
    float maxIso=0;
    maxIso = (isoType == panda::kTight) ? 0.15 : 0.25;
    return (iso < pt*maxIso);
}

inline bool ElectronIP(double eta, double dxy, double dz) {
  double aeta = fabs(eta);
  if (aeta<1.4442) {
    return (dxy < 0.05 && dz < 0.10) ;
  } else {
    return (dxy < 0.10 && dz < 0.20);
  }
}

////////////////////////////////////////////////////////////////////////////////////


inline bool IsMatched(std::vector<panda::Particle*>*objects,
               double deltaR2, double eta, double phi) {
  for (auto *x : *objects) {
    if (x->pt()>0) {
      if ( DeltaR2(x->eta(),x->phi(),eta,phi) < deltaR2 )
        return true;
    }
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////////////
/*
Bool_t RunLumiRangeMap::HasRunLumi(const std::pair<UInt_t, UInt_t> &runLumi) const
{
  // Check if a given run,lumi pair is included in the mapped lumi ranges

  //check if run is included in the map
  MapType::const_iterator it = fMap.find(runLumi.first);
  if (it!=fMap.end()) {
    //check lumis
    const MapType::mapped_type &lumiPairList = it->second;
    for (MapType::mapped_type::const_iterator jt = lumiPairList.begin(); jt<lumiPairList.end(); ++jt) {
      if (runLumi.second >= jt->first && runLumi.second <= jt->second) {
        //found lumi in accepted range
        return kTRUE;
      }
    }
  }

  return kFALSE;

}

//--------------------------------------------------------------------------------------------------
void AddJSONFile(const std::string &filepath, std::vector<std::pair<UInt_t,UInt_t> > lumiPairList) 
{

  //read json file into boost property tree
  boost::property_tree::ptree jsonTree;
  boost::property_tree::read_json(filepath,jsonTree);
  
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
  if (0) {
    printf("Iterating over parsed JSON:\n");
    for (MapType::const_iterator it = fMap.begin(); it != fMap.end(); ++it) {
      printf("  Run %u:\n",it->first);
      for (MapType::mapped_type::const_iterator jt = it->second.begin(); jt < it->second.end(); ++jt) {
        printf("    Lumis %u - %u\n",jt->first,jt->second);
      }
    }

  }

}
*/
