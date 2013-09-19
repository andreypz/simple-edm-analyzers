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
  TFile *outFile;
  //TTree *fTree;
  TTree *ktree;

  TLorentzVector *p_lminus, *p_lplus, *p_gamma;
};


GenAnalyzer::GenAnalyzer(const edm::ParameterSet& iConfig)
{}


GenAnalyzer::~GenAnalyzer()
{
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

    Int_t zStarId = 23; //normal
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
      if (abs(it->pdgId()) == 13 && it->mother()->pdgId() == zStarId)
        {
          if (it->pdgId() == 13)
            mu_plus.push_back(TLorentzVector(it->px(), it->py(), it->pz(), it->energy()));
          if (it->pdgId() == -13)
            mu_minus.push_back(TLorentzVector(it->px(), it->py(), it->pz(), it->energy()));
          ++nmu;
        }

      if (it->pdgId() == zStarId && it->mother()->pdgId()==25) //
        {          
          //cout<<"its a Z, it's status = "<<it->status()<<";\t  it's mother is: "<<it->mother()->pdgId()<<endl;
          hists->fill1DHist(it->mass(),"z_mass","Z mass", 50,0,150,  1, "");
        }

      if (it->pdgId() == 22 && (it->mother()->pdgId() == 25 || it->mother()->pdgId() == zStarId)) //gamma from the higgs or a Z!
        {          
          //cout<<"its a gamma, it's status = "<<it->status()<<";\t  it's mother is: "<<it->mother()->pdgId()<<endl;
          hists->fill1DHist(it->mass(),"gamma_mass","gamma mass", 100,0,10,  1, "");
          gamma = TLorentzVector(it->px(), it->py(), it->pz(), it->energy()); 
        }
      

    }
    
    hists->fill1DHist(nmu,"nmu","Number of muons from a Z", 10,0,10,  1, "");
    hists->fill1DHist(nel,"nel","Number of electrons from a Z", 10,0,10,  1, "");
    
    
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

    hists->fill1DHist((l1+l2).M(),"ll_mass","ll mass", 50,0,150,  1, "");

    if (nmu>=2)
      //if (nmu>=2 || nel>=2)
      {
        z = l1 + l2;
        h = z + gamma;

        //Boosting int CM frame
        TLorentzVector beamAxis(0,0,1,1);

        TVector3 b1  = -1*h.BoostVector();  //this is the boost vector

        TLorentzVector h_CM(h);
        h_CM.Boost(b1);
        TLorentzVector z_CM(z);
        z_CM.Boost(b1);
        TLorentzVector beamAxis_CM(beamAxis);
        beamAxis_CM.Boost(b1);

        //Rotating the CM frame:
        TVector3 axis_z_CM = z_CM.Vect().Unit();
        TVector3 axis_y_CM = beamAxis_CM.Vect().Cross(z_CM.Vect()).Unit();
        TVector3 axis_x_CM = axis_y_CM.Cross(axis_z_CM).Unit();
        TRotation rotation;
        rotation = rotation.RotateAxes(axis_x_CM, axis_y_CM, axis_z_CM).Inverse();

        TLorentzVector l1_CM(l1), l2_CM(l2);
        l1_CM.Boost(b1);
        l1_CM.Transform(rotation);
        l2_CM.Boost(b1);
        l2_CM.Transform(rotation); 

        Float_t costheta = z_CM.Vect().Unit().Dot(h.Vect().Unit());

        //cout<<"\t\t z_CM before rotation  x = "<<z_CM.X()<<"  y="<<z_CM.Y()<<"  z="<<z_CM.Z()<<"  t="<<z_CM.T()<<endl; 
        z_CM.Transform(rotation);
        //cout<<"\t\t z_CM after rotation  x = "<<z_CM.X()<<"  y="<<z_CM.Y()<<"  z="<<z_CM.Z()<<"  t="<<z_CM.T()<<endl; 

        TLorentzVector hRot(h);
        
        hRot.Transform(rotation);
        Float_t c3 = hRot.CosTheta();
        Float_t c2 = cos(z_CM.Angle(hRot.Vect()));    

        cout<<"higgs not bosstedd  x = "<<h.X()<<"  y="<<h.Y()<<"  z="<<h.Z()<<"  t="<<h.T()<<endl; 
        cout<<" higgs bosstedd  x = "<<h_CM.X()<<"  y="<<h_CM.Y()<<"  z="<<h_CM.Z()<<"  t="<<h_CM.T()<<endl; 
        cout<<"\t\t z bosstedd  x = "<<z_CM.X()<<"  y="<<z_CM.Y()<<"  z="<<z_CM.Z()<<"  t="<<z_CM.T()<<endl; 


          
        cout<<"Debug \n"<<"\t\t c1 = "<<costheta<<"\n \t\t c2 = "<<c2<<"\n \t\t c3 = "<<c3<<endl;

        if (fabs(costheta)>1 || isnan(costheta))
          throw cms::Exception("This is exceptional! ")<<"Cosine is greater than 1 or nan,  wha?"<<endl;
        

        //Boosting to Z frame from a CM frame!
        TVector3 b2 = -1*z_CM.BoostVector();

        TLorentzVector l1_inZFrame(l1_CM);
        TLorentzVector l2_inZFrame(l2_CM);
        l1_inZFrame.Boost(b2);
          
        //Angle between the lepton(+) and a Z direction in the Z rest frame
        Float_t costheta1 = l1_inZFrame.CosTheta();
        Float_t phi       = l1_inZFrame.Phi();
        //cout<<"l1 in Z frame: px = "<<l1_inZFrame.Px()<<"  py = "<<l1_inZFrame.Py()<<endl;
        //cout<<"l2 in Z frame: px = "<<l2_inZFrame.Px()<<"  py = "<<l2_inZFrame.Py()<<endl;

        
        c2 = (l1_CM.E() - l2_CM.E())/ (l1_CM.Vect() + l2_CM.Vect()).Mag(); //from a formula
  
          
        //if (abs(costheta1)>1 || isnan(costheta1))
        // throw cms::Exception("This is exceptional! ")<<"Cosine is greater than 1,  wha?"<<endl;



        //Boosting directly to Z
        TVector3 b3  = -1*z.BoostVector();
        TLorentzVector l1_ZFrame(l1);
        TLorentzVector l2_ZFrame(l2);
        l1_ZFrame.Boost(b3);
        l2_ZFrame.Boost(b3);

        //cout<<"\t l1 in Z frame: px = "<<l1_ZFrame.Px()<<"  py = "<<l1_ZFrame.Py()<<endl;
        //cout<<"\t l2 in Z frame: px = "<<l2_ZFrame.Px()<<"  py = "<<l2_ZFrame.Py()<<endl;

        //c3 = l1_ZFrame.Vect().Unit().Dot(b3.Unit());
        c3 = cos(TMath::Pi() - l1_ZFrame.Angle(b3));
        cout<<"Debug \n"<<"\t\t c1 = "<<costheta1<<"\n \t\t c2 = "<<c2<<"\n \t\t c3 = "<<c3<<"\n \t\t c4 = "<<endl;
        //costheta1 = c3;

        
        TString ll("ll"), lepton("lepton");
        if (nmu>=2)
          ll="mumu";
        if (nel>=2)
          ll="ee";
        hists->fill1DHist(h.M(),    Form("%s_higgs_mass", ll.Data()),"H mass", 50,0,200,  1, "");
        hists->fill1DHist(z.M(),    Form("%s_mass", ll.Data()),"Z mass", 50,0,150,  1, "");

        Int_t nbins = 30;
        hists->fill1DHist(costheta, Form("%s_costheta",  ll.Data()),"Polar Angle between Z and Zll direction in Zll center of mass", nbins,-1,1,  1, "");
        hists->fill1DHist(costheta1,Form("%s_costheta1", ll.Data()),"Polar Angle between one lepton and Z direction", nbins,-1,1,  1, "");
        hists->fill1DHist(phi,      Form("%s_phi", ll.Data()),"Azimuthal angle between one lepton and Z direction", nbins, -TMath::Pi(),TMath::Pi(),  1, "");

        hists->fill2DHist(costheta1, c3,  Form("%s_c1_vs_c3", ll.Data()),"cos theta vs cos theta :-)", 50,-1,1, 50,-1,1,  1, "");


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
        ktree->Fill();
      }
    
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
