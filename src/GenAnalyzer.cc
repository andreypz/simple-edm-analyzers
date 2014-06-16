// -*- C++ -*-
//
// Package:    GenAnalyzer
// Class:      GenAnalyzer
//
/**\class GenAnalyzer GenAnalyzer.cc NWU/GenAnalyzer/src/GenAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andrey Pozdnyakov
//         Created:  Wed Jun 12 11:55:07 CDT 2013
// $Id: GenAnalyzer.cc,v 1.2 2013/06/12 17:22:20 andrey Exp $
//
//


// system include files
#include <memory>
#include <iostream>
#include <fstream>
#include <string>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TTree.h"
#include "TLorentzVector.h"
#include "HistManager.h"
#include "ZGAngles.h"

using namespace std;
using namespace edm;
//typedef math::XYZTLorentzVector LorentzVector;/
//typedef math::XYZVector Vector;

class GenAnalyzer : public edm::EDAnalyzer {
   public:
      explicit GenAnalyzer(const edm::ParameterSet&);
      ~GenAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------

  HistManager *hists;
  ZGAngles *angles;
  TFile *outFile;
  //TTree *fTree;
  TTree *ktree;

  TLorentzVector *p_lminus, *p_lplus, *p_gamma;
};


GenAnalyzer::GenAnalyzer(const edm::ParameterSet& iConfig)
{}


GenAnalyzer::~GenAnalyzer()
{
  delete angles;
  delete hists;
  delete outFile;
}


void GenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if(!iEvent.isRealData()) {
    cout<<"event  "<<iEvent.id().event()<<endl;

    Handle<reco::GenParticleCollection> genParticleColl;
    iEvent.getByLabel("genParticles", genParticleColl);

    vector<TLorentzVector> mu_plus, mu_minus, el_plus, el_minus;
    TLorentzVector gamma;
    TLorentzVector jpsi;
    Int_t zStarId = 23; //normal
    //Int_t zStarId = 443; //jpsi
    //Int_t zStarId = 3000001; //a hack

    UInt_t nmu=0, nel = 0;
    for (reco::GenParticleCollection::const_iterator it = genParticleColl->begin(); it != genParticleColl->end(); ++it) {
      if (abs(it->pdgId()) == 11 &&   it->mother()->pdgId() == zStarId ) //for pythia6
        {
          //cout<<"its an electrobn, it's status = "<<it->status()<<";\t  it's mother is: "<<it->mother()->pdgId()<<endl;
          if (it->pdgId()==11)
            el_plus.push_back(TLorentzVector(it->px(), it->py(), it->pz(), it->energy()));
          if (it->pdgId()==-11)
            el_minus.push_back(TLorentzVector(it->px(), it->py(), it->pz(), it->energy()));
          ++nel;
        }

      if (it->pdgId() == 443)
        {
          cout<<"* \t Hooray!, we found a jPsi; its status = "<<it->status()<<";\t  its mother is: "<<it->mother()->pdgId()<<endl;
          jpsi = TLorentzVector(it->px(), it->py(), it->pz(), it->energy());

        }


      if (abs(it->pdgId()) == 13 && abs(it->mother()->pdgId()) < 24 && it->status()==3)
        //if (abs(it->pdgId()) == 13 && it->mother()->pdgId() == 443)
        {
          cout<<nmu<<"  * \t \t Yes!, we found a muon; its status = "<<it->status()<<";\t  its mother is: "<<it->mother()->pdgId()<<endl;
          if (it->pdgId() == 13)
            mu_plus.push_back(TLorentzVector(it->px(), it->py(), it->pz(), it->energy()));
          if (it->pdgId() == -13)
            mu_minus.push_back(TLorentzVector(it->px(), it->py(), it->pz(), it->energy()));
          ++nmu;
        }

      if (it->pdgId() == zStarId && it->mother()->pdgId()==25) //
        //if (it->pdgId() == zStarId && it->mother()->pdgId()==25) //
        {
          cout<<"It is a Z, its status = "<<it->status()<<";\t  its mother is: "<<it->mother()->pdgId()<<endl;
          hists->fill1DHist(it->mass(),"z_mass","Z mass", 50,0,150,  1, "");
        }

      if (it->pdgId() == 22 && abs(it->mother()->pdgId()) < 24 && it->status()==3)
        //if (it->pdgId() == 22 && it->mother()->pdgId() == 23)
        {
          cout<<"It is gamma! its status = "<<it->status()<<";\t  its mother is: "<<it->mother()->pdgId()<<endl;
          hists->fill1DHist(it->mass(),"gamma_mass_z","gamma mass", 100,-1,1,  1, "");
          gamma = TLorentzVector(it->px(), it->py(), it->pz(), it->energy());
        }


      /*
      if (it->pdgId() == 22 && (it->mother()->pdgId() == 25 || it->mother()->pdgId() == zStarId)) //gamma from the higgs or a Z!
        {
          cout<<"its a gamma, it's status = "<<it->status()<<";\t  its mother is: "<<it->mother()->pdgId()<<endl;
          hists->fill1DHist(it->mass(),"gamma_mass_h","gamma mass", 100,0,10,  1, "");
          gamma = TLorentzVector(it->px(), it->py(), it->pz(), it->energy());
        }
      */

    }

    hists->fill1DHist(nmu,"nmu","Number of muons from a Z", 10,0,10,  1, "");
    hists->fill1DHist(nel,"nel","Number of electrons from a Z", 10,0,10,  1, "");

    hists->fill1DHist((jpsi+gamma).M(),"m_jpsigamma","M_{J/Psi #gamma}", 50, 50,100,  1, "");


    TLorentzVector z(0,0,0,0);
    TLorentzVector h(0,0,0,0);
    TLorentzVector l1,l2;
    if (nmu>=2)
      {
        l1 = mu_plus[0];
        l2 = mu_minus[0];
      }
    if (nel>=2)
      {
        l1 = el_plus[0];
        l2 = el_minus[0];
      }

    TString ll("muon");

    hists->fill1DHist((l1+l2).M(),"ll_mass",    "m_{ll} (GeV)", 50,0,150, 1, "");
    hists->fill1DHist((l1+l2).M(),"ll_mass_low","m_{ll} (GeV)", 50,0,10,  1, "");

    if (nmu>=2)
      hists->fill1DHist((l1+l2+gamma).M(),"llg_mass","m_{ll#gamma} (GeV)", 50,0,150,  1, "");

    double co1,co2,phi,co3;
    angles->GetAngles(l1, l2, gamma,co1,co2,phi,co3);
    Int_t nbins = 30;
    hists->fill1DHist(co3, Form("%s_costheta", ll.Data()),"Polar Angle between Z and Zll direction in Zll center of mass", nbins,-1,1, 1, "");
    hists->fill1DHist(co1, Form("%s_costheta1", ll.Data()),"Polar Angle between positive lepton and Z direction", nbins,-1,1, 1, "");
    hists->fill1DHist(co2, Form("%s_costheta2", ll.Data()),"Polar Angle between negative lepton and Z direction", nbins,-1,1, 1, "");
    hists->fill1DHist(phi, Form("%s_phi", ll.Data()),"Azimuthal angle between one lepton and Z direction", nbins, -TMath::Pi(),TMath::Pi(), 1, "");

    hists->fill1DHist(l1.Pt(),      Form("%s_l1_pt", ll.Data()),"Pt of l+", 50, 0,150,  1, "");
    hists->fill1DHist(l2.Pt(),      Form("%s_l2_pt", ll.Data()),"Pt of l-", 50, 0,150,  1, "");
    hists->fill1DHist(l1.Eta(),     Form("%s_l1_eta", ll.Data()),"Eta of l+", 50, -3,3,  1, "");
    hists->fill1DHist(l2.Eta(),     Form("%s_l2_eta", ll.Data()),"Eta of l-", 50, -3,3,  1, "");
    hists->fill1DHist(l1.Phi(),     Form("%s_l1_phi", ll.Data()),"Phi of l+", 50, -TMath::Pi(),TMath::Pi(),  1, "");
    hists->fill1DHist(l2.Phi(),     Form("%s_l2_phi", ll.Data()),"Phi of l-", 50, -TMath::Pi(),TMath::Pi(),  1, "");
    
    hists->fill1DHist(gamma.Pt(),      Form("%s_gamma_pt", ll.Data()),"Pt of gamma", 50, 0,150,  1, "");
    hists->fill1DHist(gamma.Eta(),     Form("%s_gamma_eta", ll.Data()),"Eta of gamma", 50, -3,3,  1, "");
    hists->fill1DHist(gamma.Phi(),     Form("%s_gamma_phi", ll.Data()),"Phi of gamma", 50, -TMath::Pi(),TMath::Pi(),  1, "");
    
    
    p_lplus->SetPxPyPzE(l1.Px(), l1.Py(), l1.Pz(), l1.E());
    p_lminus->SetPxPyPzE(l2.Px(), l2.Py(), l2.Pz(), l2.E());
    p_gamma->SetPxPyPzE(gamma.Px(), gamma.Py(), gamma.Pz(), gamma.E());
        
  }

  else
    throw cms::Exception("GenAnalyzer")<<"This is supposed to be a Gen level analyzer!!"<<endl;
}

void GenAnalyzer::beginJob()
{
  outFile = new TFile("output.root","recreate");
  outFile->cd();

  ktree = new TTree("K","Kevins tree");

  p_lminus = new TLorentzVector();
  p_lplus  = new TLorentzVector();
  p_gamma  = new TLorentzVector();
  ktree->Branch("l_minus",&p_lminus,6400,0);
  ktree->Branch("l_plus",&p_lplus,6400,0);
  ktree->Branch("gamma",&p_gamma,6400,0);


  hists   = new HistManager(outFile);
  angles  = new ZGAngles();


  outFile->cd();

}

void  GenAnalyzer::endJob()
{
  outFile->cd();
  hists->writeHists(outFile);
  //ktree->Write();
  outFile->Write();
  outFile->Close();
 }

void  GenAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{ }

void  GenAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{ }

void  GenAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{ }

void  GenAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{ }

void GenAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenAnalyzer);
