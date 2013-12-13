// -*- C++ -*-
//
// Package:    EleAnalyzer
// Class:      EleAnalyzer
//
/**\class EleAnalyzer EleAnalyzer.cc NWU/EleAnalyzer/src/EleAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andrey Pozdnyakov
//         Created:  Thu Sep 19 15:13:20 CDT 2013

#include <memory>
#include <iostream>
#include <fstream>
#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "TTree.h"
#include "TLorentzVector.h"
#include "HistManager.h"

using namespace std;
using namespace edm;

class EleAnalyzer : public edm::EDAnalyzer {
   public:
      explicit EleAnalyzer(const edm::ParameterSet&);
      ~EleAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
  
  virtual bool TrackSelection(reco::Track );
  // ----------member data ---------------------------

  HistManager *hists;
  TFile *outFile;

  //TTree *ktree;

  TLorentzVector *p_lminus, *p_lplus, *p_gamma;

};
EleAnalyzer::EleAnalyzer(const edm::ParameterSet& iConfig)
{}

EleAnalyzer::~EleAnalyzer()
{
  delete hists;
  delete outFile;
}

void EleAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  if(iEvent.isRealData())
    exit(0);

  UInt_t evNum = iEvent.id().event();
  //cout<<"event  "<<iEvent.id().event()<<endl;

  Handle<reco::GenParticleCollection> genParticleColl;
  iEvent.getByLabel("genParticles", genParticleColl);

  vector<TLorentzVector> mu_plus, mu_minus, el_plus, el_minus;
  TLorentzVector gen_l1, gen_l2,gen_gamma;
  TLorentzVector gen_lPt1, gen_lPt2;

  Int_t zStarId = 23; //normal
  //Int_t zStarId = 3000001; //a hack

  UInt_t nel = 0;
  for (reco::GenParticleCollection::const_iterator it = genParticleColl->begin(); it != genParticleColl->end(); ++it) {
    if (abs(it->pdgId()) == 11 &&   it->mother()->pdgId() == zStarId ) //for pythia6
      {
        //cout<<"its an electrobn, it's status = "<<it->status()<<";\t  it's mother is: "<<it->mother()->pdgId()<<endl;
        if (it->pdgId()==11)
          {
            gen_l1 = TLorentzVector(it->px(), it->py(), it->pz(), it->energy());
            el_minus.push_back(TLorentzVector(it->px(), it->py(), it->pz(), it->energy()));
          }
        if (it->pdgId()==-11)
          {
            gen_l2 = TLorentzVector(it->px(), it->py(), it->pz(), it->energy());
            el_plus.push_back(TLorentzVector(it->px(), it->py(), it->pz(), it->energy()));
          }
        hists->fill1DHist(it->pdgId(),"pdg_el","Number of electrons from a Z", 40,-20,20,  1, "");
      }

    if (it->pdgId() == 22 && (it->mother()->pdgId() == 25 || it->mother()->pdgId() == zStarId)) //gamma from the higgs or a Z!
      {
        //cout<<"its a gamma, it's status = "<<it->status()<<";\t  it's mother is: "<<it->mother()->pdgId()<<endl;
        hists->fill1DHist(it->mass(),"gamma_mass","gamma mass", 100,-2,2,  1, "");
        gen_gamma = TLorentzVector(it->px(), it->py(), it->pz(), it->energy());
      }

  }


  if (gen_l1.Pt() > gen_l2.Pt()){
    gen_lPt1 = gen_l1;
    gen_lPt2 = gen_l2;
  }
  else{
    gen_lPt1 = gen_l2;
    gen_lPt2 = gen_l1;
  }

  nel = el_minus.size()+el_plus.size();

  hists->fill1DHist(nel,"gen_nel","Number of gen electrons from a Z", 10,0,10,  1, "");
  TLorentzVector z(0,0,0,0);
  TLorentzVector h(0,0,0,0);
  TLorentzVector l1,l2;


  if (gen_lPt1.Pt() < 20 || gen_lPt2.Pt() < 5) return;
  if (fabs(gen_lPt1.Eta()) > 2.5 || fabs(gen_lPt2.Eta()) >2.5) return;
  hists->fill1DHist(nel,"gen_nel_acc","Number of gen electrons from a Z", 10,0,10,  1, "");

  Float_t gen_Mll = (gen_l1+gen_l2).M();
  Float_t gen_qT  = (gen_l1+gen_l2).Pt();
  Float_t gen_dR  = gen_l1.DeltaR(gen_l2);

  if (gen_Mll > 2) return;

  hists->fill1DHist(nel,"gen_nel_mll","Number of gen electrons from a Z", 10,0,10,  1, "");


  edm::Handle<reco::ConversionCollection> hConversions;
  iEvent.getByLabel("allConversions", hConversions);

  if (nel>=2 && gen_gamma.Pt()!=0) {

    Handle<reco::GsfElectronCollection > electrons;
    iEvent.getByLabel("gsfElectrons", electrons);

    //Handle<reco::GsfElectronCollection > calibratedElectrons;
    //iEvent.getByLabel(edm::InputTag("calibratedElectrons","calibratedGsfElectrons"), calibratedElectrons);

    //edm::Handle<edm::ValueMap<float>> mvaTrigV0_handle;
    //iEvent.getByLabel("mvaTrigV0", mvaTrigV0_handle);
    //const edm::ValueMap<float> ele_mvaTrigV0 = (*mvaTrigV0_handle.product());

    //edm::Handle<edm::ValueMap<double>> regEne_handle;
    //iEvent.getByLabel(edm::InputTag("eleRegressionEnergy","eneRegForGsfEle"), regEne_handle);
    //const edm::ValueMap<double> ele_regEne = (*regEne_handle.product());

    //edm::Handle<edm::ValueMap<double>> regErr_handle;
    //iEvent.getByLabel(edm::InputTag("eleRegressionEnergy","eneErrorRegForGsfEle"), regErr_handle);
    //const edm::ValueMap<double> ele_regErr = (*regErr_handle.product());


    Handle<reco::TrackCollection> generalTracks;
    iEvent.getByLabel("generalTracks", generalTracks);
    
    Handle<reco::GsfTrackCollection> gsfTracks;
    iEvent.getByLabel("electronGsfTracks", gsfTracks);
    
    cout<<"Event = "<<evNum<<endl;

    cout<<"    gen1 pt="<<gen_lPt1.Pt()<<"  eta="<<gen_lPt1.Eta()<<" phi="<<gen_lPt1.Phi()<<endl;
    cout<<"    gen2 pt="<<gen_lPt2.Pt()<<"  eta="<<gen_lPt2.Eta()<<" phi="<<gen_lPt2.Phi()<<endl;
    cout<<"GEN INFO:  dR = "<<gen_lPt1.DeltaR(gen_lPt2)<<"  Mll= "<<gen_Mll<<endl;



    TLorentzVector el1,el2;

    Int_t eee=0,eee7=0,eeedr=0;
    for (vector<reco::GsfElectron>::const_iterator iElectron = electrons->begin(); iElectron != electrons->end(); ++iElectron) {
      eee++;
      hists->fill1DHist(iElectron->pt(),"el_pt","pt of an ele", 50,0,200,  1, "");
      if (iElectron->pt() < 15) continue;

      eee7++;
      if(eee7==1) 
        el1.SetXYZM(iElectron->px(), iElectron->py(), iElectron->pz(), 0);
      if(eee7==2) 
        el2.SetXYZM(iElectron->px(), iElectron->py(), iElectron->pz(), 0);


      TLorentzVector el(iElectron->px(), iElectron->py(), iElectron->pz(), iElectron->energy());
      if (el.DeltaR(gen_l1) > 0.3 && el.DeltaR(gen_l2) > 0.3) continue;
      eeedr++;
      
      hists->fill1DHist(el.DeltaR(gen_lPt1),  Form("reco_l%i_gen_lPt1_deltaR", eeedr), "reco_gen_lPt1_deltaR",   100,0,0.3, 1,"");
      hists->fill1DHist(el.DeltaR(gen_lPt2),  Form("reco_l%i_gen_lPt2_deltaR", eeedr), "reco_gen_lPt2_deltaR",   100,0,0.3, 1,"");
      //hists->fill1DHist(el.DeltaR(gen_gamma), Form("reco_l%i_gen_gamma_deltaR",eee), "reco_gen_gamma_deltaR",100,0,5, 1,"");

      reco::GsfTrackRef gsf = iElectron->gsfTrack();
      reco::TrackRef    ctf = iElectron->closestTrack(); //this returns the closest ctf track

      TLorentzVector gsf1;
      gsf1.SetXYZM(gsf->px(), gsf->py(), gsf->pz(), 0);

      cout<<eee <<" ELECTRON,  Pt = "<<iElectron->pt()<<endl;
      
      if (gsf.isNonnull())
        cout<<eee<<"  gsf pt="<<gsf->pt()<<"  eta="<<gsf->eta()<<" phi="<<gsf->phi()<<endl;
      if (ctf.isNonnull())
        cout<<eee<<"  ctf pt="<<ctf->pt()<<"  eta="<<ctf->eta()<<" phi="<<ctf->phi()<<endl;

      //TVector3 t1(gsf->px(), gsf->py(), gsf->pz());

      hists->fill1DHist(fabs(gsf->pt() - gen_lPt1.Pt())/gen_lPt1.Pt() , Form("gsf_l%i_track_ptdiff_1", eeedr), "pt diff of gsf track and gen lPt1", 100, 0,0.5, 1, "");
      hists->fill1DHist(fabs(gsf->pt() - gen_lPt2.Pt())/gen_lPt2.Pt() , Form("gsf_l%i_track_ptdiff_2", eeedr), "pt diff of gsf track and gen lPt2", 100, 0,0.5, 1, "");

      Int_t nnn03=0, nnn04=0, ntr=0;

      //reco::GsfTrackRefVector amb = iElectron->ambiguousTracks();
      TLorentzVector amb1,amb2;

      for (reco::GsfTrackRefVector::const_iterator gtr = iElectron->ambiguousGsfTracksBegin(); gtr != iElectron->ambiguousGsfTracksEnd(); ++gtr)
        {
          ntr++;
          if (iElectron->gsfTrack()==(*gtr))
            cout<<"Same ambig track"<<endl;
          cout<<ntr<<" ambigious loop pt="<<(*gtr)->pt()<<" eta="<<(*gtr)->eta()<<" "<<" phi="<<(*gtr)->phi()<<endl;
          
          TLorentzVector t;
          t.SetXYZM((*gtr)->px(), (*gtr)->py(), (*gtr)->pz(), 0);
          hists->fill1DHist(t.DeltaR(gen_lPt1),  Form("reco_ambig%i_forEl%i_gen_lPt1_deltaR", ntr,eeedr), "reco_gen_lPt1_deltaR",   100,0,0.3, 1,"");

          if (ntr==1)
            amb1.SetXYZM((*gtr)->px(), (*gtr)->py(), (*gtr)->pz(), 0);
          if (ntr==2)
            amb2.SetXYZM((*gtr)->px(), (*gtr)->py(), (*gtr)->pz(), 0);
        

        }
      hists->fill1DHist(ntr, Form("nAmbig_GSFtracks_forEl%i", eeedr),"gsf tracks", 10,0,10,1,"");

      if (ntr>1){
        //there are ambigious tracks..
        Float_t tr_Mll = (gsf1+amb1).M();
        Float_t tr_dR  = gsf1.DeltaR(amb1);
        Float_t tr_qT  = (gsf1+amb1).Pt();
        
        hists->fillProfile(gen_Mll, tr_Mll, "gen_tr_Mll",";gen Mll;tr Mll", 100, 0,20,   0,20, 1,"");
        hists->fillProfile(gen_qT, tr_qT,   "gen_tr_qT", ";gen qT;tr qT",   100, 0,100, 0,100, 1,"");
        hists->fillProfile(gen_dR, tr_dR,   "genDr_trDr",";gen dR;tr dR",   100, 0,0.5, 0,0.5, 1,"");
      }

      nnn03=0, nnn04=0, ntr=0;
      for (vector<reco::GsfTrack>::const_iterator gtr = gsfTracks->begin(); gtr != gsfTracks->end(); ++gtr) {
        if (gtr->pt()<4) continue;
        ntr++;
        if (gtr->pt() == gsf->pt()){
          cout<<ntr<<"  <<-- this one!"<<endl;
        }
        cout<<ntr<<"  gsf loop pt="<<gtr->pt()<<"  eta="<<gtr->eta()<<" "<<" phi="<<gtr->phi()<<endl;
    
        TLorentzVector t;
        t.SetXYZM(gtr->px(), gtr->py(), gtr->pz(), 0);

        hists->fill1DHist(t.DeltaR(el), Form("gsftrack_ele%i_dR", eeedr), "track ele dR", 100, 0, 0.5, 1, "");
        if (t.DeltaR(el)< 0.4)
          nnn04++;
        if (t.DeltaR(el)< 0.3)
          nnn03++;
        
      }

      hists->fill1DHist(nnn04, Form("nGSFtracks_in04_forEl%i", eeedr),"gsf tracks in 0.4 around electron", 10,0,10,1,"");
      hists->fill1DHist(nnn03, Form("nGSFtracks_in03_forEl%i", eeedr),"gsf tracks in 0.3 around electron", 10,0,10,1,"");


      nnn03=0;nnn04=0;ntr=0;
      for (vector<reco::Track>::const_iterator gtr = generalTracks->begin(); gtr != generalTracks->end(); ++gtr) {
        if (gtr->pt() < 3) continue;
        if (!gtr->quality(reco::Track::highPurity)) continue;
        if (! EleAnalyzer::TrackSelection(*gtr)) continue;


        TLorentzVector t;
        t.SetXYZM(gtr->px(), gtr->py(), gtr->pz(), 0);
        if (t.DeltaR(el)>0.4) continue;

        ntr++;
        if (ctf.isNonnull() && gtr->pt() == ctf->pt()){
          cout<<ntr<<"  <<-- this one!"<<endl;
        }

        cout<<ntr<<"  ctf  loop pt="<<gtr->pt()<<"  eta="<<gtr->eta()<<" "<<" phi="<<gtr->phi()<<endl;
    
        
        hists->fill1DHist(t.DeltaR(el), Form("ctftrack_ele%i_dR",eeedr), "ctf track ele dR", 100, 0, 0.5, 1, "");

        if (t.DeltaR(el)< 0.4)
          nnn04++;
        if (t.DeltaR(el)< 0.3)
          nnn03++;
        
      }

      hists->fill1DHist(nnn04, "nCTFtracks_in04","ctf tracks in 0.4 around electron", 10,0,10,1,"");
      hists->fill1DHist(nnn03, "nCTFtracks_in03","ctf tracks in 0.3 around electron", 10,0,10,1,"");
      
    } //end of electrons loop


    if (eee7>1){
      Float_t el_dR = el1.DeltaR(el2);
      hists->fillProfile(gen_dR, el_dR, "genDr_elDr",";gen dR;el dR", 100, 0,0.5, 0,0.5, 1,"");
    }


    hists->fill1DHist(electrons->size(),"reco_nel","Number of reco gsf electrons", 10,0,10,  1, "");



    
    // *** Conversions ** //
    // *** *** ** *******//
    
    edm::Handle<reco::BeamSpot> bsHandle;
    iEvent.getByLabel("offlineBeamSpot", bsHandle);
    const reco::BeamSpot &beamspot = *bsHandle.product();

    Handle<reco::VertexCollection> primaryVtcs;
    iEvent.getByLabel("offlinePrimaryVertices", primaryVtcs);
    reco::VertexRef PV(primaryVtcs, 0);



    Float_t conv_vtxProb=0, conv_lxy=0, conv_lxypv=0;
    Int_t conv_nHitsMax=99;

    int iel=-1;
    for(reco::GsfElectronCollection::const_iterator gsfEle = electrons->begin(); gsfEle!=electrons->end(); ++gsfEle) {
      iel++;      
      int iconv=-1;
      for (reco::ConversionCollection::const_iterator conv = hConversions->begin(); conv!= hConversions->end(); ++conv) {
        iconv++;

        reco::Vertex vtx = conv->conversionVertex();
        if (vtx.isValid()) {

          if (ConversionTools::matchesConversion(*gsfEle, *conv)) {
            
            conv_vtxProb = TMath::Prob( vtx.chi2(), vtx.ndof() );
            math::XYZVector mom(conv->refittedPairMomentum());
            double dbsx = vtx.x() - beamspot.position().x();   
            double dbsy = vtx.y() - beamspot.position().y();
            conv_lxy = (mom.x()*dbsx + mom.y()*dbsy)/mom.rho();
            
            double dpvx = vtx.x() - PV->x();   
            double dpvy = vtx.y() - PV->y();
            conv_lxypv = (mom.x()*dpvx + mom.y()*dpvy)/mom.rho();

            conv_nHitsMax=0;
            for (std::vector<uint8_t>::const_iterator it = conv->nHitsBeforeVtx().begin(); it!=conv->nHitsBeforeVtx().end(); ++it) {

              cout<<"\t\t iconv = "<<iconv<<"  Hits before vertex: "<<(int)(*it)<<endl;
              if ((*it)>conv_nHitsMax) conv_nHitsMax = (*it);
            }

            hists->fill1DHist(conv_lxy,      "conv_lxy",     "conv_lxy",      200,0,50,1,"");
            hists->fill1DHist(conv_lxypv,    "conv_lxyPV",   "conv_lxyPV",    200,0,50,1,"");
            hists->fill1DHist(conv_nHitsMax, "conv_nHitsMax","conv_nHitsMax", 15, 0,15,1,"");
            hists->fill1DHist(conv_vtxProb,  "conv_vtxProb", "conv_vtxProv",  100,0,1, 1,"");

            break;

          }
          
        }
      }
    }

    
    /*
      double conv_vtxProb[50];
      double conv_lxy[50];
      double conv_lxypv[50];
      int conv_nHitsMax[50];
      //int conv_eleind[50];

    int iconv=-1;

    for (reco::ConversionCollection::const_iterator conv = hConversions->begin(); conv!= hConversions->end(); ++conv) {
      iconv++;
      conv_vtxProb[iconv]=0.;
      conv_lxy[iconv]=0.;
      conv_lxypv[iconv]=0.;
      conv_nHitsMax[iconv]=99;
      //conv_eleind[iconv] = -1;
      
      reco::Vertex vtx = conv->conversionVertex();
      if (vtx.isValid()) {
        int iel=-1;
        for(reco::GsfElectronCollection::const_iterator gsfEle = electrons->begin(); gsfEle!=electrons->end(); ++gsfEle) {
          iel++;
          if (ConversionTools::matchesConversion(*gsfEle, *conv)) {
            //conv_eleind[iconv] = iel;
            conv_vtxProb[iconv] = TMath::Prob( vtx.chi2(), vtx.ndof() );
            math::XYZVector mom(conv->refittedPairMomentum());
            double dbsx = vtx.x() - beamspot.position().x();   
            double dbsy = vtx.y() - beamspot.position().y();
            conv_lxy[iconv] = (mom.x()*dbsx + mom.y()*dbsy)/mom.rho();

            double dpvx = vtx.x() - PV->x();   
            double dpvy = vtx.y() - PV->y();
            conv_lxypv[iconv] = (mom.x()*dpvx + mom.y()*dpvy)/mom.rho();

            conv_nHitsMax[iconv]=0;
            for (std::vector<uint8_t>::const_iterator it = conv->nHitsBeforeVtx().begin(); it!=conv->nHitsBeforeVtx().end(); ++it) {
              if ((*it)>conv_nHitsMax[iconv]) conv_nHitsMax[iconv] = (*it);
            }

            break;
          }
          else{
            for (reco::GsfTrackRefVector::const_iterator gtr = gsfEle->ambiguousGsfTracksBegin(); gtr != gsfEle->ambiguousGsfTracksEnd(); ++gtr)
              {
                if (ConversionTools::matchesConversion(*gtr, *conv)) {
                  //conv_eleind[iconv] = iel;
                  conv_vtxProb[iconv] = TMath::Prob( vtx.chi2(), vtx.ndof() );
                  math::XYZVector mom(conv->refittedPairMomentum());
                  double dbsx = vtx.x() - beamspot.position().x();   
                  double dbsy = vtx.y() - beamspot.position().y();
                  conv_lxy[iconv] = (mom.x()*dbsx + mom.y()*dbsy)/mom.rho();
                  conv_nHitsMax[iconv]=0;
                  for (std::vector<uint8_t>::const_iterator it = conv->nHitsBeforeVtx().begin(); it!=conv->nHitsBeforeVtx().end(); ++it) {
                    if ((*it)>conv_nHitsMax[iconv]) conv_nHitsMax[iconv] = (*it);
                  }


                  hists->fill1DHist(conv_lxy[iconv],      "conv_amb_lxy",     "conv_lxy",      200,0,50,1,"");
                  hists->fill1DHist(conv_nHitsMax[iconv], "conv_amb_nHitsMax","conv_nHitsMax", 15, 0,15,1,"");
                  hists->fill1DHist(conv_vtxProb[iconv],  "conv_amb_vtxProb", "conv_vtxProv",  100,0,1, 1,"");

                  break;
                }
                
              }
          }
        }

      }
      

      
    }
    */



    Handle<vector<reco::Photon> > photons;
    iEvent.getByLabel("photons", photons);
    Int_t ppp = 0;
    for (vector<reco::Photon>::const_iterator iPhoton = photons->begin(); iPhoton != photons->end() ; ++iPhoton) {
      
      if (iPhoton->pt()>15) continue;
      TLorentzVector ph(iPhoton->px(), iPhoton->py(), iPhoton->pz(), iPhoton->p());
      
      //hists->fill1DHist(ph.DeltaR(gen_lPt1), Form("reco_ph%i_gen_lPt1_deltaR", ppp),"reco_gen_l1_deltaR",   100,0,5, 1,"");
      //hists->fill1DHist(ph.DeltaR(gen_lPt2), Form("reco_ph%i_gen_lPt2_deltaR", ppp),"reco_gen_l2_deltaR",   100,0,5, 1,"");
      //hists->fill1DHist(ph.DeltaR(gen_gamma),Form("reco_ph%i_gen_gamma_deltaR",ppp),"reco_gen_gamma_deltaR",100,0,5, 1,"");

      if (ph.DeltaR(gen_lPt1) < 0.3){
        hists->fill1DHist(iPhoton->nTrkSolidConeDR04(), "ph_fake_nTrkSolidConeDR04","fake nTrkSolidConeDR04", 10,0,10,1,"");
        hists->fill1DHist(iPhoton->nTrkSolidConeDR03(), "ph_fake_nTrkSolidConeDR03","fake nTrkSolidConeDR03", 10,0,10,1,"");
      }
      else 
        if (ph.DeltaR(gen_lPt1) > 1.5)
          {
            hists->fill1DHist(iPhoton->nTrkSolidConeDR04(), "ph_real_nTrkSolidConeDR04","real nTrkSolidConeDR04", 10,0,10,1,"");
            hists->fill1DHist(iPhoton->nTrkSolidConeDR03(), "ph_real_nTrkSolidConeDR03","real nTrkSolidConeDR03", 10,0,10,1,"");
          }

    }

  }






}


void  EleAnalyzer::beginJob()
{
  outFile = new TFile("output.root","recreate");
  outFile->cd();

  /*
    ktree = new TTree("K","Kevins tree");

  p_lminus = new TLorentzVector();
  p_lplus  = new TLorentzVector();
  p_gamma  = new TLorentzVector();
  ktree->Branch("l_minus",&p_lminus,6400,0);
  ktree->Branch("l_plus",&p_lplus,6400,0);
  ktree->Branch("gamma",&p_gamma,6400,0);

  */
  hists   = new HistManager(outFile);


  outFile->cd();

}

// ------------ method called once each job just after ending the event loop  ------------
void EleAnalyzer::endJob()
{
  outFile->cd();
  hists->writeHists(outFile);
  //ktree->Write();
  outFile->Write();
  outFile->Close();
}

void EleAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{}

void EleAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{}

void  EleAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

void EleAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{}

void EleAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}



bool EleAnalyzer::TrackSelection(reco::Track tr){
  Float_t normalizedChi2 = tr.normalizedChi2();
  Float_t numberOfValidTrackerHits = tr.hitPattern().numberOfValidTrackerHits();

  hists->fill1DHist(normalizedChi2,           "tr_normalizedChi2",          "tr_normalizedChi2",           200,0,5, 1,"");
  hists->fill1DHist(numberOfValidTrackerHits, "tr_numberOfValidTrackerHits","tr_numberOfValidTrackerHits", 30,0,30, 1,"");

  if (normalizedChi2  < 3
      && numberOfValidTrackerHits > 5)
    
    return true;

  else 
    return false;
}


//define this as a plug-in
DEFINE_FWK_MODULE(EleAnalyzer);
