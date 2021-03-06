#include "HistManager.h"

HistManager::HistManager(TFile* myFile){
  theFile = myFile;
  //theFile -> Print();
}


HistManager::~HistManager(){  
}


void HistManager::writeHists(TFile* myFile){
  
  myFile->cd();
  std::map<std::string,TH1F*>::const_iterator mapit1;
  std::map<std::string,TH2F*>::const_iterator mapit2;
  std::map<std::string,TProfile*>::const_iterator mapit3;
  for (mapit1 = the1DMap.begin(); mapit1 != the1DMap.end(); ++mapit1){
    (*mapit1).second->Write();
  }
  for (mapit2 = the2DMap.begin(); mapit2 != the2DMap.end(); ++mapit2){
    (*mapit2).second->Write();
  }
  for (mapit3 = theProfMap.begin(); mapit3 != theProfMap.end(); ++mapit3){
    (*mapit3).second->Write();
  }
  myFile->cd();
  //myFile->Close(); //Don't close it! It should be taken care of outside of this plugin!
  
}

//void HistManager::fill1DHist(float x, TString s, std::string title,
//                           int bins, float xmin, float xmax, float weight, std::string folder){
//
//HistManager::fill1DHist(float x, s.Data(), std::string title,
//           int bins, float xmin, float xmax, float weight, std::string folder);
//}

void HistManager::fill1DHist(float x, std::string name, std::string title,
                               int bins, float xmin, float xmax, float weight, std::string folder){
  
  std::map<std::string,TH1F*>::iterator it;
  it = the1DMap.find(name);
  //std::cout<<"fill 1d  "<<name<<std::endl;
  if (it == the1DMap.end()){
    //theFile->ls();
    //std::cout<<"  ** Before the cd into a directory"<<std::endl;
    theFile->cd(folder.c_str());
    //std::cout<<"  ****** After the cd into a directory"<<std::endl;
    //theFile->ls();
    the1DMap[name] = new TH1F(name.c_str(),title.c_str(),bins,xmin,xmax);
    theFile->cd();
  }
  
  the1DMap[name]->Fill(x,weight);
  
}

void HistManager::fill1DHistUnevenBins(float x, std::string name, std::string title,
                                       int bins, float *binEdges, float weight, std::string folder){
  std::map<std::string,TH1F*>::iterator it;
  it = the1DMap.find(name);
  if (it == the1DMap.end()){
    
    theFile->cd(folder.c_str());
    the1DMap[name] = new TH1F(name.c_str(),title.c_str(),bins,binEdges);
    theFile->cd();
  }
  
  
  the1DMap[name]->Fill(x,weight);
  
}




void HistManager::fill2DHist(float x, float y, std::string name, std::string title,
                             int binsx, float xmin, float xmax,
                             int binsy, float ymin, float ymax, float weight, std::string folder){
  
  std::map<std::string,TH2F*>::iterator it;
    it = the2DMap.find(name);
    if (it == the2DMap.end()){
      theFile->cd(folder.c_str());
      the2DMap[name] = new TH2F(name.c_str(),title.c_str(),binsx,xmin,xmax,binsy,ymin,ymax);
      theFile->cd();
    }
    
    the2DMap[name]->Fill(x,y,weight);
    
}

void HistManager::fill2DHistUnevenBins(float x, float y, std::string name, std::string title,
                                       int binsx, float *binEdgesx,
                                       int binsy, float *binEdgesy, float weight,  std::string folder){
  
  std::map<std::string,TH2F*>::iterator it;
  it = the2DMap.find(name);
  if (it == the2DMap.end()){
    theFile->cd(folder.c_str());
    the2DMap[name] = new TH2F(name.c_str(),title.c_str(),binsx,binEdgesx,binsy,binEdgesy);
    theFile->cd();
  }
  
  
  the2DMap[name]->Fill(x,y,weight);
}


void HistManager::fillProfile(float x, float y, std::string name, std::string title,
                              int binsx, float xmin, float xmax,
                              float ymin, float ymax, float weight,  std::string folder){
  
  std::map<std::string,TProfile*>::iterator it;
  it = theProfMap.find(name);
  if (it == theProfMap.end()){
    theFile->cd(folder.c_str());
    theProfMap[name] = new TProfile(name.c_str(),title.c_str(),binsx,xmin,xmax,ymin,ymax);
    theFile->cd();
  }
  
  theProfMap[name]->Fill(x,y,weight);
  
}


