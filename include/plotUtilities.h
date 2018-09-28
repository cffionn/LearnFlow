#ifndef PLOTUTILITIES_H
#define PLOTUTILITIES_H

//cpp dependencies
#include <string>

//ROOT dependencies
#include "TCanvas.h"
#include "TH1.h"
#include "TBox.h"
#include "TError.h"

std::string prettyString(const double inVal, const int prec, const bool doDot)
{
  std::string retStr = std::to_string(inVal);
  while(retStr.find(".") < retStr.size()-1-prec){retStr.replace(retStr.size()-1, 1,"");}
  if(doDot) retStr.replace(retStr.find("."), 1, "p");
  return retStr;
}


void prettyCanv(TCanvas* canv_p)
{
  canv_p->SetRightMargin(0.01);
  canv_p->SetLeftMargin(1.5*canv_p->GetLeftMargin());
  canv_p->SetBottomMargin(canv_p->GetLeftMargin());
  canv_p->SetTopMargin(canv_p->GetLeftMargin()/2.);

  return;
}


void prettyTH1(TH1* hist_p, const double size, const int style, const int col)
{
  hist_p->SetMarkerSize(size);
  hist_p->SetMarkerStyle(style);
  hist_p->SetMarkerColor(col);
  hist_p->GetXaxis()->CenterTitle();
  hist_p->GetYaxis()->CenterTitle();

  return;
}


void drawWhiteBox(Double_t x1, Double_t x2, Double_t y1, Double_t y2)
{
  TBox* tempBox_p = new TBox();
  tempBox_p->SetFillColor(0);
  tempBox_p->DrawBox(x1, y1, x2, y2);
  delete tempBox_p;
}


void quietSaveAs(TCanvas* canv_p, const std::string saveName)
{
  Int_t oldLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kWarning;
  canv_p->SaveAs(saveName.c_str());
  gErrorIgnoreLevel = oldLevel;
  
  return;
}


#endif
