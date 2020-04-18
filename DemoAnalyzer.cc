// -*- C++ -*-
//
// Package:    Demo/DemoAnalyzer
// Class:      DemoAnalyzer
//
/**\class DemoAnalyzer DemoAnalyzer.cc Demo/DemoAnalyzer/plugins/DemoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Dener De Souza Lemos
//         Created:  Wed, 12 Feb 2020 13:05:50 GMT
//
//



#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//v0 candidates
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"

#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"

//
// class declaration
//
#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

//
#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>
#include "TVector3.h"
#include <vector>
#include <map>

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class DemoAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit DemoAnalyzer(const edm::ParameterSet&);
      ~DemoAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      static bool vtxSort( const reco::Vertex &  a, const reco::Vertex & b );
      
	  void initHistos(const edm::Service<TFileService> & fs);  

      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
      edm::EDGetTokenT<reco::TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
      edm::EDGetTokenT<reco::Centrality> CentralityTag_;
      edm::EDGetTokenT<int> CentralityBinTag_;
  	  edm::EDGetTokenT<CaloTowerCollection> _towerSrc;
  	  
  	  
  	  double vtx_x, vtx_y, vtx_z, vtx_rho, vtx_xError, vtx_yError, vtx_zError;  //primary vertex
  	  int hiBin;
  	  double hiHF, hiNpix, hiNtracks;
     //add histograms
     
     TH1D* hiHF_hist;
     TH1D* etHFtowerSum_hist;
     TH2D* hiHF_vs_etHFtowerSum;

	 TH1D*hiNpix_hist;
	 TH1D*hiNtracks_hist;	
     TH2D* hiHF_vs_hiNpix;
     TH2D* hiHF_vs_hiNtracks;

     TH1D* hiHF_plus_hist;
     TH1D* hiHF_minus_hist;
     TH2D* hiHF_plus_vs_minus;

//control plots

     TH1D* trk_pt;
     TH1D* trk_eta;
     TH1D* trk_phi;

     TH1D* tower_pt;
     TH1D* tower_eta;
     TH1D* tower_phi;
     

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig)
 :
  vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"))),
  tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
  CentralityTag_(consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("CentralitySrc"))),
  CentralityBinTag_(consumes<int>(iConfig.getParameter<edm::InputTag>("CentralityBinSrc"))),
  _towerSrc(consumes<CaloTowerCollection>(iConfig.getParameter<edm::InputTag>("towerSrc")))
{
   //now do what ever initialization is needed
}


DemoAnalyzer::~DemoAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;
   using namespace std;

//Vertex selection

  edm::Handle<std::vector<reco::Vertex> > vertexCollection;
  iEvent.getByToken(vertexToken_,vertexCollection);
  std::vector<reco::Vertex> vtx_sorted = *vertexCollection;
//  const reco::Vertex & vtx_0 = (*vertexCollection)[0];
  std::sort( vtx_sorted.begin(), vtx_sorted.end(), DemoAnalyzer::vtxSort );
  if(vtx_sorted.size() == 0) return;
  vtx_x = (double)vtx_sorted.begin()->position().x();// VertexX_->Fill(vtx_x);
  vtx_y = (double)vtx_sorted.begin()->position().y();// VertexY_->Fill(vtx_y);
  vtx_z = (double)vtx_sorted.begin()->position().z(); 
  vtx_rho= (double)vtx_sorted.begin()->position().Rho(); 
  vtx_xError = (double)vtx_sorted.begin()->xError();
  vtx_yError = (double)vtx_sorted.begin()->yError();
  vtx_zError = (double)vtx_sorted.begin()->zError();
  math::XYZPoint vtxx(vtx_x,vtx_y,vtx_z);

  if(vtx_z > 15.)return;

  double etHFtowerSumPlus=0;
  double etHFtowerSumMinus=0;
  double etHFtowerSum=0;
  Handle<CaloTowerCollection> towers;
  iEvent.getByToken(_towerSrc,towers);
  for( size_t i = 0; i<towers->size(); ++ i){
  const CaloTower & tower = (*towers)[ i ];
  double eta = tower.eta();
  bool isHF = tower.ietaAbs() > 29;
  if(isHF && eta > 0){
  etHFtowerSumPlus += tower.pt();
  tower_pt->Fill(tower.pt());
  tower_eta->Fill(eta);
  tower_phi->Fill(tower.phi());
  }
  if(isHF && eta < 0){
  etHFtowerSumMinus += tower.pt();
  tower_pt->Fill(tower.pt());
  tower_eta->Fill(eta);
  tower_phi->Fill(tower.phi());
  }
  }


  etHFtowerSum=etHFtowerSumPlus + etHFtowerSumMinus;   
  etHFtowerSum_hist->Fill(etHFtowerSum);
  
  edm::Handle<int> cbin_;
  iEvent.getByToken(CentralityBinTag_,cbin_);
  hiBin = *cbin_;
  
  edm::Handle<reco::Centrality> centrality;
  iEvent.getByToken(CentralityTag_, centrality);

  hiHF = centrality->EtHFtowerSum();
  hiHF_hist->Fill(hiHF);
  hiHF_vs_etHFtowerSum->Fill(hiHF,etHFtowerSum);

  hiHF_plus_hist->Fill(centrality->EtHFtowerSumPlus());
  hiHF_minus_hist->Fill(centrality->EtHFtowerSumMinus());
  hiHF_plus_vs_minus->Fill(centrality->EtHFtowerSumPlus(),centrality->EtHFtowerSumMinus());


  hiNpix = centrality->NpixelTracks();
  hiNpix_hist->Fill(hiNpix);
  hiHF_vs_hiNpix->Fill(hiHF,hiNpix);
  
  hiNtracks = centrality->Ntracks();
  hiNtracks_hist->Fill(hiNtracks);
  hiHF_vs_hiNtracks->Fill(hiHF,hiNtracks);
  
  
    Handle<TrackCollection> tracks;
    iEvent.getByToken(tracksToken_, tracks);
    for(TrackCollection::const_iterator iter_tk = tracks->begin(); iter_tk != tracks->end();++iter_tk) {
      double aux_tk_dz_vtx = (double)iter_tk->dz(vtxx);
      double aux_tk_dzError_vtx  = (double)sqrt(iter_tk->dzError()*iter_tk->dzError()+vtx_zError*vtx_zError);
      double aux_tk_dxy_vtx = (double)iter_tk->dxy(vtxx);
      double aux_tk_dxyError_vtx  = (double)sqrt(iter_tk->dxyError()*iter_tk->dxyError()+vtx_xError*vtx_yError);
      if(iter_tk->pt()<0.5)continue;
      if(fabs(iter_tk->eta())>2.4)continue;
      if(!iter_tk->quality(reco::TrackBase::highPurity))continue;
      if(fabs(iter_tk->ptError())/iter_tk->pt()>0.1)continue;
      if(fabs(aux_tk_dz_vtx/aux_tk_dzError_vtx)>3)continue;
      if(fabs(aux_tk_dxy_vtx/aux_tk_dxyError_vtx)>3)continue;
      if(iter_tk->numberOfValidHits()<11)continue;
      if((iter_tk->normalizedChi2()/iter_tk->hitPattern().trackerLayersWithMeasurement())>0.18)continue; 
	  trk_pt->Fill(iter_tk->pt());
	  trk_eta->Fill(iter_tk->eta());
	  trk_phi->Fill(iter_tk->phi());
    }  
}


// ------------ method called once each job just before starting event loop  ------------
void
DemoAnalyzer::beginJob()
{
  std::cout<<"  This is called once for each job: Begin Job " << std::endl;
  edm::Service<TFileService> fs;
  initHistos(fs);
}

// ------------ method called once each job just after ending the event loop  ------------
void
DemoAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

bool DemoAnalyzer::vtxSort( const reco::Vertex &  a, const reco::Vertex & b ){
   if( a.tracksSize() != b.tracksSize() )
      return  a.tracksSize() > b.tracksSize() ? true : false ;
   else
      return  a.chi2() < b.chi2() ? true : false ;
}

void DemoAnalyzer::initHistos(const edm::Service<TFileService> & fs){

   TFileDirectory Inf = fs->mkdir( "STARLight_Hist" );

   hiHF_hist = Inf.make<TH1D>("hiHF_hist","",200,0.0,100.0);
   etHFtowerSum_hist = Inf.make<TH1D>("etHFtowerSum_hist","",200,0.0,100.0);
   hiHF_vs_etHFtowerSum = Inf.make<TH2D>("hiHF_vs_etHFtowerSum","",200,0.0,100.0,200,0.0,100.0);

   hiNpix_hist = Inf.make<TH1D>("hiNpix_hist","",100,0.0,100.0);
   hiNtracks_hist = Inf.make<TH1D>("hiNtracks_hist","",100,0.0,100.0);	
   hiHF_vs_hiNpix = Inf.make<TH2D>("hiHF_vs_hiNpix","",200,0.0,100.0,100,0.0,100.0);
   hiHF_vs_hiNtracks = Inf.make<TH2D>("hiHF_vs_hiNtracks","",200,0.0,100.0,100,0.0,100.0);

   hiHF_plus_hist = Inf.make<TH1D>("hiHF_plus_hist","",200,0.0,100.0);
   hiHF_minus_hist = Inf.make<TH1D>("hiHF_minus_hist","",200,0.0,100.0);
   hiHF_plus_vs_minus = Inf.make<TH2D>("hiHF_plus_vs_minus","",200,0.0,100.0,200,0.0,100.0);

//control plots

   trk_pt=Inf.make<TH1D>("trk_pt","",200,0.0,10.0);
   trk_eta=Inf.make<TH1D>("trk_eta","",48,-2.4,2.4);
   trk_phi=Inf.make<TH1D>("trk_phi","",32,-3.2,3.2);

   tower_pt=Inf.make<TH1D>("tower_pt","",200,0.0,10.0);
   tower_eta=Inf.make<TH1D>("tower_eta","",60,-6.0,6.0);
   tower_phi=Inf.make<TH1D>("tower_phi","",32,-3.2,3.2);


  
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);
