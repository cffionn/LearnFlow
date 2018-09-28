//cpp dependencies
#include <iostream>
#include <string>

//ROOT dependencies
#include "TFile.h"
#include "TTree.h"

//Local dependencies
#include "include/checkMakeDir.h"

int learnRho(const std::string inFileName)
{
  if(!checkFile(inFileName)){
    std::cout << "Given inFileName \'" << inFileName << "\' is not valid. return 1" << std::endl;
    return 1;
  }

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");

  inFile_p->Close();
  delete inFile_p;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/learnRho.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += learnRho(argv[1]);
  return retVal;
}
