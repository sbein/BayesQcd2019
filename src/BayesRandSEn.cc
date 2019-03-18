#include <cmath>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include "TMinuit.h"
#include "TLorentzVector.h"
#include <time.h>
#include <algorithm>
#include <TH1.h>
#include <unistd.h>
#include "src/UsefulJet.h"

using namespace std;


void findJetToPin(std::vector<UsefulJet> jetVec, int nparams, int & ipin, double & cstart )
{
  double ht = getHT(jetVec, JET_PT_THRESH);
  double desiredMht = TMath::Min(120.0, ht/3);// /3 worked better than /2 fyi, for HT300
  //double desiredMht = 120;
  TLorentzVector mhtvec = getMHT(jetVec,30);
  double mhtPt = mhtvec.Pt(); double mhtPhi = mhtvec.Phi();
  //cout << "starting on MHT=" << mhtPt << endl;
  if (mhtPt<desiredMht) {
    ipin = -1; cstart = 1.0;
    return;  
  }
  //cout << "high met event, " << mhtvec.Pt() << endl;
  for (unsigned int i = 0; i<nparams; i++)
    {
      double denom = 2*(-pow(desiredMht,2)+pow(jetVec[i].Pt(),2)+pow(mhtPt,2)+2*(jetVec[i].Pt())*
  			mhtPt*TMath::Cos(mhtPhi-jetVec[i].Phi()));
      double num = 2*pow(jetVec[i].Pt(),2)+2*jetVec[i].Pt()*mhtPt*TMath::Cos(mhtPhi-jetVec[i].Phi());
      double discriminant = 2.0*pow(jetVec[i].Pt(),2)*(2*pow(desiredMht,2)+pow(mhtPt,2)*
  						       (TMath::Cos(2*(mhtPhi-jetVec[i].Phi()))-1));
      double num1 = num+TMath::Sqrt(discriminant);
      double num2 = num-TMath::Sqrt(discriminant);
      double c1 = num1/denom;
      if (fabs(c1-1)<0.8){
  	ipin = i; cstart = c1;
  	return;
      }
      double c2 = num2/denom;
      if (fabs(c2-1)<0.8){
  	ipin = i; cstart = c2;/////////change in hamburg////////
  	return;
      }
    }
  cout << "couldn't find jet to pin, returning 1" << endl;
  ipin = -2; cstart = 1.0;
  return;
}


TemplateSet _Templates_;//global struct
double *_par_; 
UsefulJet _leadjet_; double _iLeadJet_ = 0; double _LeadJetPt_ = 0; 
int _nbjets_ = 0; int _ibjet_;
TGraph * _ActiveMhtPrior_, * _ActiveMhtPhiPrior_; 
TLorentzVector _ActiveMht_; double _ActiveHt_; int _iht_;
double _ActiveLikelihood_; 
double _ptNew_, _ptBinC_, _PtBinOther_, _interpolatedFactor_;
double _a_, _b_; 
int _ietaNew_, _iptNew_, _otherbin_;
double _rmht_, _rdphi_;


void getLeadingBJet_CC(std::vector<UsefulJet> RecoJets, int & ibjet, int & nbjets, UsefulJet & leadingbjet){
  nbjets = 0;
  ibjet = -1;
  leadingbjet = UsefulJet();
  double highestPt = 0;
  for (unsigned int ireco = 0; ireco < RecoJets.size(); ireco++){
    if (!(RecoJets[ireco].csv>BTAG_CSV)) continue;
    double pt = RecoJets[ireco].Pt();
    if (!(pt>_Templates_.lhdMhtThresh)) continue;//////CHANGE IN HAMBURG////////
    //but first, lets see why this is crashing!
    nbjets+=1;
    if (!(pt>highestPt)) continue;
    highestPt = pt;
    ibjet = ireco;
    leadingbjet = RecoJets[ireco].Clone();
  }
  return;
}


void MakeTemplatesGlobal(TemplateSet newtemplates)
{
  _Templates_ = newtemplates;
}

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  _ActiveLikelihood_ = 1.0;
  _iLeadJet_ = 0; _LeadJetPt_ = 0; _nbjets_ = 0; 
  _ActiveMht_.SetPtEtaPhiE(0,0,0,0); _ActiveHt_ = 0;
  for (int i = 0; i < _Templates_.nparams; i++)
    {
      _Templates_.dynamicJets[i]*=(_par_[i]*1.0/par[i]);//my attempt to scale the jets only once per step.
      _par_[i]=par[i];//my attempt to scale the jets only once per step.
      if(_Templates_.dynamicJets[i].Pt()>_LeadJetPt_)
	{ _iLeadJet_ = i; _LeadJetPt_ = _Templates_.dynamicJets[i].Pt();}
      if(fabs(_Templates_.dynamicJets[i].Eta())<2.4)
	{
	  if (_Templates_.dynamicJets[i].Pt()>JET_PT_THRESH) _ActiveHt_+= _Templates_.dynamicJets[i].Pt(); 
	  if(_Templates_.dynamicJets[i].csv>BTAG_CSV &&
	     _Templates_.dynamicJets[i].Pt()>_Templates_.lhdMhtThresh ) _nbjets_+=1;
	}    
      if(_Templates_.dynamicJets[i].Pt()>_Templates_.lhdMhtThresh)
	_ActiveMht_-=_Templates_.dynamicJets[i].tlv;
      //_ptNew_ = TMath::Min(1799.0, _Templates_.dynamicJets[i].Pt());
      _ptNew_ = TMath::Min(1799.0, _Templates_.dynamicJets[i].E());//Hacking to use energy not momentum
      _ietaNew_ = TMath::Min(_Templates_.hEtaTemplate->GetXaxis()->FindBin(fabs(_Templates_.dynamicJets[i].Eta())),13);//this min here until mcresolutions is fixed XXX
      _iptNew_ = _Templates_.hPtTemplate->GetXaxis()->FindBin(_ptNew_);
      _ptBinC_ = _Templates_.hPtTemplate->GetXaxis()->GetBinCenter(_iptNew_);
      _otherbin_ = -1;
      if (_ptNew_-_ptBinC_>0) _otherbin_ = 1;
      _PtBinOther_ = _Templates_.hPtTemplate->GetXaxis()->GetBinCenter(_iptNew_+_otherbin_);
      _a_ = (_PtBinOther_-_ptNew_)/(_PtBinOther_ - _ptBinC_);
      _b_ = (_ptNew_-_ptBinC_)/(_PtBinOther_ - _ptBinC_);
      try
	{
	  _interpolatedFactor_ = 0.5*(_a_*_Templates_.ResponseFunctions[_Templates_.dynamicJets[i].csv>BTAG_CSV][_ietaNew_][_iptNew_]->Eval(par[i],0,"S") +
				      _b_*_Templates_.ResponseFunctions[_Templates_.dynamicJets[i].csv>BTAG_CSV][_ietaNew_][_iptNew_+_otherbin_]->Eval(par[i],0,"S"));
	}
      catch (std::exception& e)
	{
	  _interpolatedFactor_ = _b_*_Templates_.ResponseFunctions[_Templates_.dynamicJets[i].csv>BTAG_CSV][_ietaNew_][_iptNew_+_otherbin_]->Eval(par[i],0,"S");
	}
      _ActiveLikelihood_*=_interpolatedFactor_;
    }
  _ActiveHt_ = TMath::Min(4999.0,_ActiveHt_);
  _iht_ = _Templates_.hHtTemplate->GetXaxis()->FindBin(_ActiveHt_);

  if (_nbjets_ == 0)
    {
      _leadjet_ = _Templates_.dynamicJets[_iLeadJet_];
      _ActiveMhtPrior_ = _Templates_.gGenMhtPtTemplatesB0[_iht_];
      _ActiveMhtPhiPrior_ = _Templates_.gGenMhtDPhiTemplatesB0[_iht_];
    }
  else 
    {
      getLeadingBJet_CC(_Templates_.dynamicJets, _ibjet_, _nbjets_, _leadjet_);
      if(_nbjets_ == 1)
	{
	  _ActiveMhtPrior_ = _Templates_.gGenMhtPtTemplatesB1[_iht_];
	  _ActiveMhtPhiPrior_ = _Templates_.gGenMhtDPhiTemplatesB1[_iht_];
	}
      else if(_nbjets_ == 2)
	{
	  _ActiveMhtPrior_ = _Templates_.gGenMhtPtTemplatesB2[_iht_];
	  _ActiveMhtPhiPrior_ = _Templates_.gGenMhtDPhiTemplatesB2[_iht_];
	}
      else
	{
	  _ActiveMhtPrior_ = _Templates_.gGenMhtPtTemplatesB3[_iht_];
	  _ActiveMhtPhiPrior_ = _Templates_.gGenMhtDPhiTemplatesB3[_iht_];
	}
    }
  _rmht_ = _ActiveMht_.Pt();
  _rdphi_ = _ActiveMht_.DeltaPhi(_leadjet_.tlv);
  _ActiveLikelihood_*=_ActiveMhtPrior_->Eval(_rmht_,0,"S");
  _ActiveLikelihood_*=_ActiveMhtPhiPrior_->Eval(_rdphi_,0,"S");
  f = -fabs(_ActiveLikelihood_);
  return;

}


bool RebalanceJets_BayesFitter(std::vector<UsefulJet> originalJets){

  _Templates_.dynamicJets.clear();
  _Templates_.nparams = 0;

  for (unsigned int j = 0; j < originalJets.size(); j++)
    {
      if (originalJets[j].Pt()>_Templates_.lhdMhtThresh && j<12) _Templates_.nparams+=1; 
      _Templates_.dynamicJets.push_back(originalJets[j].Clone());
    }

  _par_ = new double[_Templates_.nparams]; for(int i=0; i<_Templates_.nparams; i++) _par_[i] = 1.0;

  int ipin; double cstart;
  findJetToPin(_Templates_.dynamicJets, _Templates_.nparams, ipin, cstart);
  //cout << "found jet to pin, " << ipin << ", " << cstart << endl;

  TMinuit * gMinuit = new TMinuit(_Templates_.nparams);
  //gMinuit->SetPrintLevel(0);

  gMinuit->SetFCN(fcn);
  gMinuit->SetPrintLevel(-1);
  Double_t arglist[_Templates_.nparams];
  Int_t ierflg = 0;//Long ierflg = 0;
  if (_Templates_.nparams>0) arglist[0] = 0;

  for(int i = 0; i < _Templates_.nparams; i++)
    {
      char buffer [50]; int name;
      name = sprintf (buffer, "c%d", i);
      //double lowerbound = 0.3; double upperbound = 3.0;
      double lowerbound = 0.01; double upperbound = 5.0;
      //make identical to diss python:
      lowerbound = 0.3; upperbound = 3.5;
      if(i==ipin) gMinuit->mnparm(i, buffer, cstart, 0.05, lowerbound, upperbound, ierflg);
      else  gMinuit->mnparm(i,buffer,1.0,0.05,lowerbound,upperbound,ierflg);
    }
  gMinuit->SetMaxIterations(10000);
  if (_Templates_.nparams>0) arglist[0] = 10000;
  if (_Templates_.nparams>1) arglist[1] = 1;
  gMinuit->mnexcm( "MINIMIZE", arglist, 2, ierflg );
  if (ierflg!=0){
    cout << "failed fit" << endl;
    TLorentzVector omhtVec = getMHT(originalJets, JET_PT_THRESH);
    TLorentzVector rmhtVec = getMHT(_Templates_.dynamicJets, JET_PT_THRESH);
    return false;
  }

  double currentValue(1), currentError(1);
  for (int i=0; i<_Templates_.nparams; i++){
    gMinuit->GetParameter (i, currentValue, currentError);
    _Templates_.dynamicJets[i]*=(originalJets[i].Pt()/_Templates_.dynamicJets[i].Pt()/currentValue);
  }

  // TLorentzVector omhtVec = getMHT(originalJets, JET_PT_THRESH);
  // double omht = omhtVec.Pt();
  // double odphi = omhtVec.Phi();
  //cout << "old mht = " << omht << ", " << odphi << endl;

  TLorentzVector rmhtVec = getMHT(_Templates_.dynamicJets, JET_PT_THRESH);
  double rmht = rmhtVec.Pt();
  double rdphi = rmhtVec.Phi();
  //cout << "new mht = " << rmht << ", " << rdphi << endl;
  
  std::sort(_Templates_.dynamicJets.begin(), _Templates_.dynamicJets.end());

  return true;
}


std::vector<UsefulJet> smearJets_CC(std::vector<UsefulJet> jetVec, int n2smear){
  std::vector<UsefulJet> smearedJets;
  for (unsigned int j=0; j<jetVec.size(); j++){
    smearedJets.push_back(jetVec[j].Clone());
    double pt = smearedJets.back().Pt(); 
    double eta = smearedJets.back().Eta(); 
    if (pt<8 || int(j)>=n2smear) continue;
    int ieta = TMath::Min(_Templates_.hEtaTemplate->GetXaxis()->FindBin(fabs(eta)), 13); //until the new JECs arrive XXX
    int ipt = _Templates_.hPtTemplate->GetXaxis()->FindBin(pt);
    if (!(_Templates_.ResponseHistos.at(smearedJets.back().csv>BTAG_CSV).at(ieta).at(ipt)->Integral()==0))
      {
	double rando = _Templates_.ResponseHistos.at(smearedJets.back().csv>BTAG_CSV).at(ieta).at(ipt)->GetRandom();
	smearedJets.back()*=rando;
      }
  }
  std::sort(smearedJets.begin(), smearedJets.end());
  return smearedJets;
}




void GleanTemplatesFromFile(TFile* ftemplate)
{

  //TH1F* hPtTemplate = (TH1F*)ftemplate->Get("hEnTemplate");//Partial name change now
  TH1F* hPtTemplate = (TH1F*)ftemplate->Get("hPtTemplate");//Partial name change now
  TAxis* templatePtAxis = (TAxis*)hPtTemplate->GetXaxis();
  TH1F* hEtaTemplate = (TH1F*)ftemplate->Get("hEtaTemplate");
  TAxis* templateEtaAxis = (TAxis*)hEtaTemplate->GetXaxis();
  TH1F* hHtTemplate = (TH1F*)ftemplate->Get("hHtTemplate");
  TAxis* templateHtAxis = (TAxis*)hHtTemplate->GetXaxis();

  TH1F* hNull = new TH1F("hNull","hNull",1,-4,-3);
  std::vector<TH1F* >  hNullVec;
  hNullVec.push_back(hNull);
  TGraph* gNull = new TGraph();
  std::vector<TGraph*>  gNullVec;
  gNullVec.push_back(gNull);

  std::vector<std::vector<TH1F*> > hResTemplates_CC;
  hResTemplates_CC.push_back(hNullVec);
  std::vector<std::vector<TH1F*> > hResTemplatesB_CC;
  hResTemplatesB_CC.push_back(hNullVec);

  std::vector<std::vector<TH1F*> > hRebTemplates_CC;
  hRebTemplates_CC.push_back(hNullVec);
  std::vector<std::vector<TH1F*> > hRebTemplatesB_CC;
  hRebTemplatesB_CC.push_back(hNullVec);

  std::vector<std::vector<TGraph*> > gRebTemplates_CC;
  gRebTemplates_CC.push_back(gNullVec);
  std::vector<std::vector<TGraph*> > gRebTemplatesB_CC;
  gRebTemplatesB_CC.push_back(gNullVec);

  int name;
  for(unsigned int ieta = 1; ieta<=templateEtaAxis->GetNbins()+1; ieta++)
    {
      hResTemplates_CC.push_back(hNullVec);
      hRebTemplates_CC.push_back(hNullVec);
      gRebTemplates_CC.push_back(gNullVec);
      hResTemplatesB_CC.push_back(hNullVec);
      hRebTemplatesB_CC.push_back(hNullVec);
      gRebTemplatesB_CC.push_back(gNullVec);
      for (unsigned int ipt=1; ipt <= templatePtAxis->GetNbins()+1; ipt++)
    	{	
    	  char hname[100];
    	  //name = sprintf (hname, "hRTemplate(gEn%2.1f-%2.1f, gEta%2.1f-%2.1f)", templatePtAxis->GetBinLowEdge(ipt), templatePtAxis->GetBinUpEdge(ipt), templateEtaAxis->GetBinLowEdge(ieta), templateEtaAxis->GetBinUpEdge(ieta));
    	  name = sprintf (hname, "hRTemplate(gPt%2.1f-%2.1f, gEta%2.1f-%2.1f)", templatePtAxis->GetBinLowEdge(ipt), templatePtAxis->GetBinUpEdge(ipt), templateEtaAxis->GetBinLowEdge(ieta), templateEtaAxis->GetBinUpEdge(ieta));
    	  TH1F* h = (TH1F*)ftemplate->Get(hname);
    	  hResTemplates_CC.back().push_back(h);
    	  hRebTemplates_CC.back().push_back(h);//once upon a time these were rPt"s

    	  //name = sprintf (hname, "hRTemplate(gEn%2.1f-%2.1f, gEta%2.1f-%2.1f)B", templatePtAxis->GetBinLowEdge(ipt), templatePtAxis->GetBinUpEdge(ipt), templateEtaAxis->GetBinLowEdge(ieta), templateEtaAxis->GetBinUpEdge(ieta));
    	  name = sprintf (hname, "hRTemplate(gPt%2.1f-%2.1f, gEta%2.1f-%2.1f)B", templatePtAxis->GetBinLowEdge(ipt), templatePtAxis->GetBinUpEdge(ipt), templateEtaAxis->GetBinLowEdge(ieta), templateEtaAxis->GetBinUpEdge(ieta));
    	  TH1F* hB = (TH1F*)ftemplate->Get(hname);
    	  hResTemplatesB_CC.back().push_back(hB);//trying regular templates
    	  hRebTemplatesB_CC.back().push_back(hB);//was hB

    	  char gname[100];
    	  //name = sprintf (gname, "splines/hRTemplate(gEn%2.1f-%2.1f, gEta%2.1f-%2.1f)_graph", templatePtAxis->GetBinLowEdge(ipt), templatePtAxis->GetBinUpEdge(ipt), templateEtaAxis->GetBinLowEdge(ieta), templateEtaAxis->GetBinUpEdge(ieta));
    	  name = sprintf (gname, "splines/hRTemplate(gPt%2.1f-%2.1f, gEta%2.1f-%2.1f)_graph", templatePtAxis->GetBinLowEdge(ipt), templatePtAxis->GetBinUpEdge(ipt), templateEtaAxis->GetBinLowEdge(ieta), templateEtaAxis->GetBinUpEdge(ieta));
    	  TGraph* g = (TGraph*)ftemplate->Get(gname);
    	  gRebTemplates_CC.back().push_back(g);

    	  //name = sprintf (gname, "splines/hRTemplate(gEn%2.1f-%2.1f, gEta%2.1f-%2.1f)B_graph", templatePtAxis->GetBinLowEdge(ipt), templatePtAxis->GetBinUpEdge(ipt), templateEtaAxis->GetBinLowEdge(ieta), templateEtaAxis->GetBinUpEdge(ieta));
    	  name = sprintf (gname, "splines/hRTemplate(gPt%2.1f-%2.1f, gEta%2.1f-%2.1f)B_graph", templatePtAxis->GetBinLowEdge(ipt), templatePtAxis->GetBinUpEdge(ipt), templateEtaAxis->GetBinLowEdge(ieta), templateEtaAxis->GetBinUpEdge(ieta));
    	  TGraph* gB = (TGraph*)ftemplate->Get(gname);
    	  gRebTemplatesB_CC.back().push_back(gB);//Instead of gB
    	}
    }


  std::vector<TGraph*> gGenMhtPtTemplatesB0_CC;
  gGenMhtPtTemplatesB0_CC.push_back(gNull);
  std::vector<TGraph*> gGenMhtPtTemplatesB1_CC;
  gGenMhtPtTemplatesB1_CC.push_back(gNull);
  std::vector<TGraph*> gGenMhtPtTemplatesB2_CC;
  gGenMhtPtTemplatesB2_CC.push_back(gNull);
  std::vector<TGraph*> gGenMhtPtTemplatesB3_CC;
  gGenMhtPtTemplatesB3_CC.push_back(gNull);
  std::vector<TGraph*> gGenMhtDPhiTemplatesB0_CC;
  gGenMhtDPhiTemplatesB0_CC.push_back(gNull);
  std::vector<TGraph*> gGenMhtDPhiTemplatesB1_CC;
  gGenMhtDPhiTemplatesB1_CC.push_back(gNull);
  std::vector<TGraph*> gGenMhtDPhiTemplatesB2_CC;
  gGenMhtDPhiTemplatesB2_CC.push_back(gNull);
  std::vector<TGraph*> gGenMhtDPhiTemplatesB3_CC;
  gGenMhtDPhiTemplatesB3_CC.push_back(gNull);

  string keyvar = "Mht";//This could be changed later if the MET is to be used

  for(unsigned int iht = 1; iht < templateHtAxis->GetNbins()+2; iht++)
    {
      char gname[100];
      name = sprintf (gname, "splines/hGen%sPtB0(ght%2.1f-%2.1f)_graph", keyvar.c_str(), templateHtAxis->GetBinLowEdge(iht), templateHtAxis->GetBinUpEdge(iht));
      TGraph* fb0 = (TGraph*)ftemplate->Get(gname);
      gGenMhtPtTemplatesB0_CC.push_back(fb0);
      char gnamePhi[100];
      name = sprintf (gnamePhi, "splines/hGen%sPhiB0(ght%2.1f-%2.1f)_graph", keyvar.c_str(), templateHtAxis->GetBinLowEdge(iht), templateHtAxis->GetBinUpEdge(iht));
      TGraph* fb0phi = (TGraph*)ftemplate->Get(gnamePhi);
      gGenMhtDPhiTemplatesB0_CC.push_back(fb0phi);
      name = sprintf (gname, "splines/hGen%sPtB1(ght%2.1f-%2.1f)_graph", keyvar.c_str(), templateHtAxis->GetBinLowEdge(iht), templateHtAxis->GetBinUpEdge(iht));
      TGraph* fb1 = (TGraph*)ftemplate->Get(gname);
      gGenMhtPtTemplatesB1_CC.push_back(fb1);
      name = sprintf (gnamePhi, "splines/hGen%sPhiB1(ght%2.1f-%2.1f)_graph", keyvar.c_str(), templateHtAxis->GetBinLowEdge(iht), templateHtAxis->GetBinUpEdge(iht));
      TGraph* fb1phi = (TGraph*)ftemplate->Get(gnamePhi);
      gGenMhtDPhiTemplatesB1_CC.push_back(fb1phi);
      name = sprintf (gname, "splines/hGen%sPtB2(ght%2.1f-%2.1f)_graph", keyvar.c_str(), templateHtAxis->GetBinLowEdge(iht), templateHtAxis->GetBinUpEdge(iht));
      TGraph* fb2 = (TGraph*)ftemplate->Get(gname);
      gGenMhtPtTemplatesB2_CC.push_back(fb2);
      name = sprintf (gnamePhi, "splines/hGen%sPhiB2(ght%2.1f-%2.1f)_graph", keyvar.c_str(), templateHtAxis->GetBinLowEdge(iht), templateHtAxis->GetBinUpEdge(iht));
      TGraph* fb2phi = (TGraph*)ftemplate->Get(gnamePhi);
      gGenMhtDPhiTemplatesB2_CC.push_back(fb2phi);
      name = sprintf (gname, "splines/hGen%sPtB3(ght%2.1f-%2.1f)_graph", keyvar.c_str(), templateHtAxis->GetBinLowEdge(iht), templateHtAxis->GetBinUpEdge(iht));
      TGraph* fb3 = (TGraph*)ftemplate->Get(gname);
      gGenMhtPtTemplatesB3_CC.push_back(fb3);
      name = sprintf (gnamePhi, "splines/hGen%sPhiB3(ght%2.1f-%2.1f)_graph", keyvar.c_str(), templateHtAxis->GetBinLowEdge(iht), templateHtAxis->GetBinUpEdge(iht));
      TGraph* fb3phi = (TGraph*)ftemplate->Get(gnamePhi);
      gGenMhtDPhiTemplatesB3_CC.push_back(fb3phi);
    }

  _Templates_.hEtaTemplate = hEtaTemplate;
  _Templates_.hPtTemplate = hPtTemplate;//_Templates_.hResTemplates
  _Templates_.hHtTemplate = hHtTemplate;
  _Templates_.ResponseFunctions.clear();
  _Templates_.ResponseFunctions.push_back(gRebTemplates_CC);
  _Templates_.ResponseFunctions.push_back(gRebTemplatesB_CC);

  _Templates_.ResponseHistos.clear();
  _Templates_.ResponseHistos.push_back(hResTemplates_CC);
  _Templates_.ResponseHistos.push_back(hResTemplatesB_CC);
  _Templates_.gGenMhtPtTemplatesB0 = gGenMhtPtTemplatesB0_CC;
  _Templates_.gGenMhtPtTemplatesB1 = gGenMhtPtTemplatesB1_CC;
  _Templates_.gGenMhtPtTemplatesB2 = gGenMhtPtTemplatesB2_CC;
  _Templates_.gGenMhtPtTemplatesB3 = gGenMhtPtTemplatesB3_CC;
  _Templates_.gGenMhtDPhiTemplatesB0 = gGenMhtDPhiTemplatesB0_CC;
  _Templates_.gGenMhtDPhiTemplatesB1 = gGenMhtDPhiTemplatesB1_CC;
  _Templates_.gGenMhtDPhiTemplatesB2 = gGenMhtDPhiTemplatesB2_CC;
  _Templates_.gGenMhtDPhiTemplatesB3 = gGenMhtDPhiTemplatesB3_CC;
  _Templates_.lhdMhtThresh = lhdMhtThresh;
  _Templates_.nparams = 0;
}


std::vector<double> createMatchedCsvVector(std::vector<TLorentzVector> GenJets, std::vector<UsefulJet> RecoJets)
{
  std::vector<double> matchedCsvs;
  for (unsigned int ig = 0; ig<GenJets.size(); ig++)
    {
      double csv = 0;
      double dR = 0.6;
      for (unsigned int ir = 0; ir<RecoJets.size(); ir++)
	{
	  double dR_ = GenJets[ig].DeltaR(RecoJets[ir].tlv);
	  if (dR_<dR)
	    {
	      dR=dR_;
	      csv = RecoJets[ir].csv;
	    }
	  if (dR_ < 0.25) break;
	}
      matchedCsvs.push_back(csv);
    }
  return matchedCsvs;
}

/*
*/
