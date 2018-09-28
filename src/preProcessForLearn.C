//cpp dependencies
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"
#include "TDatime.h"
#include "TH1D.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TDirectoryFile.h"

//Local dependencies
#include "include/checkMakeDir.h"
#include "include/etaPhiFunc.h"
#include "include/getLinBins.h"
#include "include/plotUtilities.h"

int learnFlowTMVA(const std::string inFileName)
{
  if(!checkFile(inFileName)){
    std::cout << "Error: Given inFileName \'" << inFileName << "\' is invalid. return 1."  << std::endl;
    return 1;
  }

  TRandom3* randGen_p = new TRandom3(0);  
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("output");
  std::string outFileName = inFileName;
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  while(outFileName.find(".root") != std::string::npos){outFileName.replace(outFileName.find(".root"), std::string(".root").size(), "");}
  outFileName = "output/" + outFileName + "_FlowTMVA_" + dateStr + ".root";

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1D* genPt_h = new TH1D("genPt_h", ";Gen. p_{T};Counts", 100, 0.5, 100.5);
  TH1D* genEta_h = new TH1D("genEta_h", ";Gen. #eta;Counts", 50, -5., 5.);

  const Int_t nPixelEtaBins = 9;
  const Double_t pixelEtaBinsNeg[nPixelEtaBins+1] = {-5.1, -4.7, -4.5, -4.3, -4.1, -3.9, -3.7, -3.5, -3.3, -3.1};
  const Double_t pixelEtaBinsPos[nPixelEtaBins+1] = {3.1, 3.3, 3.5, 3.7, 3.9, 4.1, 4.3, 4.5, 4.7, 5.1};

  const Int_t nPixelPhiBins = 18;
  Double_t pixelPhiBins[nPixelPhiBins+1];
  getLinBins(-TMath::Pi(), TMath::Pi(), nPixelPhiBins, pixelPhiBins);

  std::cout << "PI=" << TMath::Pi() << std::endl;
  
  const Int_t nRCEtaStrips = 7;
  const Double_t rcEtaStrips[nRCEtaStrips+1] = {-5.1, -3.0, -2.1, -1.3, 1.3, 2.1, 3.0, 5.1};
  const Int_t nRCPerEtaStrip = 3;
  const Double_t rcR = 0.4;

  Double_t maxChargePartPt_Eta2p5_;
  Double_t evtPlane2_;
  Double_t rho_[nRCEtaStrips][nRCPerEtaStrip];
  Double_t rhoEta_[nRCEtaStrips][nRCPerEtaStrip];
  Double_t rhoPhi_[nRCEtaStrips][nRCPerEtaStrip];
  Double_t ptSum_[2*nPixelEtaBins][nPixelPhiBins];
  Double_t ptSumFlat_[2*nPixelEtaBins*nPixelPhiBins];
  
  TTree* pixelTree_p = new TTree("pixelTree", "");
  pixelTree_p->Branch("maxChargePartPt_Eta2p5", &maxChargePartPt_Eta2p5_, "maxChargePartPt_Eta2p5/D");
  pixelTree_p->Branch("evtPlane2", &evtPlane2_, "evtPlane2/D");
  pixelTree_p->Branch("ptSum", ptSum_, ("ptSum[" + std::to_string(2*nPixelEtaBins) + "][" + std::to_string(nPixelPhiBins) + "]/D" ).c_str());
  pixelTree_p->Branch("ptSumFlat", ptSumFlat_, ("ptSumFlat[" + std::to_string(2*nPixelEtaBins*nPixelPhiBins) + "]/D" ).c_str());
  pixelTree_p->Branch("rho", rho_, ("rho[" + std::to_string(nRCEtaStrips) + "][" + std::to_string(nRCPerEtaStrip) + "]/D").c_str());
  pixelTree_p->Branch("rhoEta", rhoEta_, ("rhoEta[" + std::to_string(nRCEtaStrips) + "][" + std::to_string(nRCPerEtaStrip) + "]/D").c_str());
  pixelTree_p->Branch("rhoPhi", rhoPhi_, ("rhoPhi[" + std::to_string(nRCEtaStrips) + "][" + std::to_string(nRCPerEtaStrip) + "]/D").c_str());
    
  std::vector<float>* genPt_p = NULL;
  std::vector<float>* genEta_p = NULL;
  std::vector<float>* genPhi_p = NULL;
  
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* genTree_p = (TTree*)inFile_p->Get("genTree");
  
  genTree_p->SetBranchStatus("*", 0);
  genTree_p->SetBranchStatus("genPt", 1);
  genTree_p->SetBranchStatus("genEta", 1);
  genTree_p->SetBranchStatus("genPhi", 1);

  genTree_p->SetBranchAddress("genPt", &(genPt_p));
  genTree_p->SetBranchAddress("genEta", &(genEta_p));
  genTree_p->SetBranchAddress("genPhi", &(genPhi_p));

  const Int_t nEntries = genTree_p->GetEntries();
  const Int_t nDiv = TMath::Max(1, nEntries/20);

  std::cout << "Processing " << nEntries << " events..." << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;
    genTree_p->GetEntry(entry);
    
    for(Int_t peI = 0; peI < 2*nPixelEtaBins; ++peI){
      for(Int_t ppI = 0; ppI < nPixelPhiBins; ++ppI){
	ptSum_[peI][ppI] = 0.0;
      }
    }

    for(Int_t reI = 0; reI < nRCEtaStrips; ++reI){
      std::vector<double> tempEtaVect;
      std::vector<double> tempPhiVect;

      while(tempEtaVect.size() < (unsigned int)nRCPerEtaStrip){
	double tempEta = randGen_p->Uniform(rcEtaStrips[reI], rcEtaStrips[reI+1]);
	double tempPhi = randGen_p->Uniform(-TMath::Pi(), TMath::Pi());

	bool isGood = true;
	for(unsigned int teI = 0; teI < tempEtaVect.size(); ++teI){
	  if(getDR(tempEtaVect.at(teI), tempPhiVect.at(teI), tempEta, tempPhi) < rcR){
	    isGood = false;
	    break;
	  }
	}

	if(!isGood) continue;

	tempEtaVect.push_back(tempEta);
	tempPhiVect.push_back(tempPhi);
      }
      
      for(Int_t rpI = 0; rpI < nRCPerEtaStrip; ++rpI){
	rho_[reI][rpI] = 0;
	rhoEta_[reI][rpI] = tempEtaVect.at(rpI);
	rhoPhi_[reI][rpI] = tempPhiVect.at(rpI);
      }

      tempEtaVect.clear();
      tempPhiVect.clear();
    }

    Double_t tempMaxChargePartPt_Eta2p5_ = -1;
    Double_t evtPlane2Cos_ = 0;
    Double_t evtPlane2Sin_ = 0;
    
    const Int_t nGen_ = genPt_p->size();
    for(Int_t gI = 0; gI < nGen_; ++gI){
      genPt_h->Fill(genPt_p->at(gI));
      genEta_h->Fill(genEta_p->at(gI));

      if(TMath::Abs(genEta_p->at(gI)) < 2.5 &&  genPt_p->at(gI) > tempMaxChargePartPt_Eta2p5_){
	tempMaxChargePartPt_Eta2p5_ = genPt_p->at(gI);
      }
      
      if(TMath::Abs(genEta_p->at(gI)) < 1. && genPt_p->at(gI) >= .3 && genPt_p->at(gI) < 3.){
	evtPlane2Cos_ += TMath::Cos(genPhi_p->at(gI));
	evtPlane2Sin_ += TMath::Sin(genPhi_p->at(gI));	
      }
      
      Int_t partEtaPos = -1;
      for(Int_t reI = 0; reI < nRCEtaStrips; ++reI){
	if(rcEtaStrips[reI] <= genEta_p->at(gI) && rcEtaStrips[reI+1] > genEta_p->at(gI)){
	  partEtaPos = reI;
	  break;
	}
      }
      if(partEtaPos == -1){
	if(genEta_p->at(gI) == rcEtaStrips[nRCEtaStrips]) partEtaPos = nRCEtaStrips-1;
	else{
	  std::cout << "Eta pos of \'" << genEta_p->at(gI) << " not found return 1" << std::endl;
	  return 1;
	}
      }

      for(Int_t rcI = 0; rcI < nRCPerEtaStrip; ++rcI){
	if(getDR(rhoEta_[partEtaPos][rcI], rhoPhi_[partEtaPos][rcI], genEta_p->at(gI), genPhi_p->at(gI)) < rcR){
	  rho_[partEtaPos][rcI] += genPt_p->at(gI);
	  break;
	}
      }
      
      if(TMath::Abs(genEta_p->at(gI)) < pixelEtaBinsPos[0]) continue;
      
      Int_t etaPos = -1;
      Int_t phiPos = -1;

      for(Int_t pI = 0; pI < nPixelPhiBins; ++pI){
	if(genPhi_p->at(gI) >= pixelPhiBins[pI] && genPhi_p->at(gI) < pixelPhiBins[pI+1]){
	  phiPos = pI;
	  break;
	}
      }
      
      if(phiPos == -1){
	if(genPhi_p->at(gI) == pixelPhiBins[nPixelPhiBins]) phiPos = nPixelPhiBins-1;
	else{
	  std::cout << "Warning: phi \'" << genPhi_p->at(gI) << "\' is not found. continuing" << std::endl;
	  continue;	 
	}
      }

      if(genEta_p->at(gI) < 0){
	for(Int_t pI = 0; pI < nPixelEtaBins; ++pI){
	  if(genEta_p->at(gI) >= pixelEtaBinsNeg[pI] && genEta_p->at(gI) < pixelEtaBinsNeg[pI+1]){
	    etaPos = pI;
	    break;
	  }
	}
	
	if(etaPos == -1){
	  if(genEta_p->at(gI) == pixelEtaBinsNeg[nPixelEtaBins]) etaPos = nPixelEtaBins-1;
	  else{
	    std::cout << "Warning: eta \'" << genEta_p->at(gI) << "\' is not found. continuing" << std::endl;
	    continue;	 
	  }
	}
      }
      else{
	for(Int_t pI = 0; pI < nPixelEtaBins; ++pI){
	  if(genEta_p->at(gI) >= pixelEtaBinsPos[pI] && genEta_p->at(gI) < pixelEtaBinsPos[pI+1]){
	    etaPos = pI;
	    break;
	  }
	}
	
	if(etaPos == -1){
	  if(genEta_p->at(gI) == pixelEtaBinsPos[nPixelEtaBins]) etaPos = nPixelEtaBins-1;
	  else{
	    std::cout << "Warning: eta \'" << genEta_p->at(gI) << "\' is not found. continuing" << std::endl;
	    continue;	 
	  }
	}
      }
      
      ptSum_[etaPos][phiPos] += genPt_p->at(gI);
    }

    unsigned int pos = 0;
    for(Int_t eI = 0; eI < 2*nPixelEtaBins; ++eI){
      for(Int_t pI = 0; pI < nPixelPhiBins; ++pI){
	ptSumFlat_[pos] = ptSum_[eI][pI];
	++pos;
      }
    }

    maxChargePartPt_Eta2p5_ = tempMaxChargePartPt_Eta2p5_;
    
    evtPlane2_ = TMath::ATan2(evtPlane2Cos_, evtPlane2Sin_);

    for(Int_t reI = 0; reI < nRCEtaStrips; ++reI){
      for(Int_t rpI = 0; rpI < nRCPerEtaStrip; ++rpI){
	rho_[reI][rpI] /= TMath::Pi()*rcR*rcR;
      }

      for(Int_t rpI = 0; rpI < nRCPerEtaStrip; ++rpI){
	for(Int_t rpI2 = rpI+1; rpI2 < nRCPerEtaStrip; ++rpI2){
	  if(rho_[reI][rpI] < rho_[reI][rpI2]){
	    Double_t rhoTemp = rho_[reI][rpI2];
	    Double_t rhoEtaTemp = rhoEta_[reI][rpI2];
	    Double_t rhoPhiTemp = rhoPhi_[reI][rpI2];

	    rho_[reI][rpI2] = rho_[reI][rpI];
	    rhoEta_[reI][rpI2] = rhoEta_[reI][rpI];
	    rhoPhi_[reI][rpI2] = rhoPhi_[reI][rpI];

	    rho_[reI][rpI] = rhoTemp;
	    rhoEta_[reI][rpI] = rhoEtaTemp;
	    rhoPhi_[reI][rpI] = rhoPhiTemp;
	  }
	}
      }
    }
    
    pixelTree_p->Fill();
  }
  
  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();
  
  genPt_h->Write("", TObject::kOverwrite);
  genEta_h->Write("", TObject::kOverwrite);

  pixelTree_p->Write("", TObject::kOverwrite);
  
  delete genPt_h;
  delete genEta_h;

  delete pixelTree_p;

  TDirectoryFile* paramsDir_p = (TDirectoryFile*)outFile_p->mkdir("params");
  paramsDir_p->cd();

  std::string pixelEtaBinsNegStr = "";
  std::string pixelEtaBinsPosStr = "";
  for(Int_t eI = 0; eI < nPixelEtaBins+1; ++eI){
    pixelEtaBinsNegStr = pixelEtaBinsNegStr + prettyString(pixelEtaBinsNeg[eI], 1, false) + ",";
    pixelEtaBinsPosStr = pixelEtaBinsPosStr + prettyString(pixelEtaBinsPos[eI], 1, false) + ",";
  }
  pixelEtaBinsNegStr.replace(pixelEtaBinsNegStr.size()-1, 1, "");
  pixelEtaBinsPosStr.replace(pixelEtaBinsPosStr.size()-1, 1, "");

  std::string pixelPhiBinsStr = "";
  for(Int_t pI = 0; pI < nPixelPhiBins+1; ++pI){
    if(pI == 0 || pI == nPixelPhiBins) pixelPhiBinsStr = pixelPhiBinsStr + prettyString(pixelPhiBins[pI], 5, false) + ",";
    else pixelPhiBinsStr = pixelPhiBinsStr + prettyString(pixelPhiBins[pI], 2, false) + ",";
  }
  pixelPhiBinsStr.replace(pixelPhiBinsStr.size()-1, 1, "");
  
  std::string rcEtaStripsStr = "";
  for(Int_t rcI = 0; rcI < nRCEtaStrips+1; ++rcI){
    rcEtaStripsStr = rcEtaStripsStr + prettyString(rcEtaStrips[rcI], 1, false) + ",";
  }
  rcEtaStripsStr.replace(rcEtaStripsStr.size()-1, 1, "");

  TNamed nPixelEtaBinsName("nPixelEtaBins", std::to_string(nPixelEtaBins).c_str());
  TNamed pixelEtaBinsNegName("pixelEtaBinsNeg", pixelEtaBinsNegStr.c_str());
  TNamed pixelEtaBinsPosName("pixelEtaBinsPos", pixelEtaBinsPosStr.c_str());
  
  TNamed nPixelPhiBinsName("nPixelPhiBins", std::to_string(nPixelPhiBins).c_str());
  TNamed pixelPhiBinsName("pixelPhiBins", pixelPhiBinsStr.c_str());

  TNamed nRCEtaStripsName("nRCEtaStrips", std::to_string(nRCEtaStrips).c_str());
  TNamed rcEtaStripsName("rcEtaStrips", rcEtaStripsStr.c_str());
  TNamed nRCPerEtaStripName("nRCPerEtaStrip", std::to_string(nRCPerEtaStrip).c_str());
  TNamed rcRName("rcR", prettyString(rcR, 1, false).c_str());

  nPixelEtaBinsName.Write("", TObject::kOverwrite);
  pixelEtaBinsNegName.Write("", TObject::kOverwrite);
  pixelEtaBinsPosName.Write("", TObject::kOverwrite);

  nPixelPhiBinsName.Write("", TObject::kOverwrite);
  pixelPhiBinsName.Write("", TObject::kOverwrite);
  
  nRCEtaStripsName.Write("", TObject::kOverwrite);
  rcEtaStripsName.Write("", TObject::kOverwrite);
  nRCPerEtaStripName.Write("", TObject::kOverwrite);
  rcRName.Write("", TObject::kOverwrite);
  
  outFile_p->Close();
  delete outFile_p;

  delete randGen_p;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/learnFlowTMVA.exe <inFileName>" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += learnFlowTMVA(argv[1]);
  return retVal;
}
