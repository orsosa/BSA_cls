#define BSA_survey_cls_cxx
#include "BSA_survey_cls.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void BSA_survey_cls::Loop()
{
//   In a ROOT session, you can do:
//      root> .L BSA_survey_cls_sim.C
//      root> BSA_survey_cls_sim t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0) return;

  std::cout<<"processing...\n";
  std::cout.fill('.');
   
  Long64_t nentries = fChain->GetEntries();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry = 0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    TH1D *h;
    if (!(npart == 1)) continue;

    if ( !(DIS() && eFID_ec() && eFID_dc()) ) continue;
    for (int k = 0; k<npart; k++){
      if ( !(FWD(k) && CF(k) && piFID_ec(k) && piFID_dc(k)) ) continue;
      ofile->GetObject("hM",h);
      h->Fill(M[k]);
      
      if (helic == -1) 
	ofile->GetObject("hp_phiR",h);
      else if(helic == 1)
	ofile->GetObject("hn_phiR",h);
      h->Fill(phiR[k]);

      if (helic == -1) 
	ofile->GetObject("hp_phiH",h);
      else if(helic == 1)
	ofile->GetObject("hn_phiH",h);
      h->Fill(pdata_phiHs[k][0]);

      for (auto& x: bedg){
	TString bn = x.first;
	for (int n=0;n<x.second.size()-1;n++){
	  TString hname = pltv + "_" + bn + Form("_b%d",n);
	  TString ttlsuf =  Form("%.2f<%s<%.2f",x.second[n], bn.Data(), x.second[n+1]);
	  std::cout<<x.second[n]<<"<"<<pdata_phiHs[k][0]<<"<"<<x.second[n+1]<<std::endl;
	  if (x.second[n] < brv[bn][k]&& brv[bn][k]< x.second[n+1]){
	    if (helic == -1) 
	      ofile->GetObject("hp_"+ hname,h);
	    else if(helic == 1)
	      ofile->GetObject("hn_"+ hname,h);
	    h->Fill(pdata_phiHs[k][0]);
	    std::cout<<x.second[n]<<"<"<<pdata_phiHs[k][0]<<"<"<<x.second[n+1]<<std::endl;
	    break;
	  }
	}
      }

    }

    // if (Cut(ientry) < 0) continue;
    std::cout<<std::setw(15)<<float(jentry+1)/nentries*100<<" %"<<"\r";
    std::cout.flush();
      
  }

  std::cout<<getALU("hp_phiR","hn_phiR",pltv,ttlv)<<std::endl;
  std::cout<<getALU("hp_phiH","hn_phiH","phiH","#phi_{H}")<<std::endl;
  for (auto& x: bedg){
    TString bn = x.first;
    for (int k=0;k<x.second.size()-1;k++){
      TString hname = pltv + "_" + bn + Form("_b%d",k);
      TString ttlsuf =  Form("%.2f<%s<%.2f",x.second[k], bn.Data(), x.second[k+1]);
      std::cout<<getALU("hp_"+ hname,"hn_"+ hname,"phiH","#phi_{H}")<<std::endl;
    }
  }
  
  /*
  TH1D *hp,*hn, *hs, *hm, *hALU;
  TFitResultPtr res;
  TF1 *ff = new TF1("ff"," [A]*sin(x*TMath::DegToRad())",0,360);
  
  /// phiR asymmetry ///
  ofile->GetObject("hp_phiR",hp);
  ofile->GetObject("hn_phiR",hn);
  hp->Sumw2();
  hn->Sumw2();
  hs = (TH1D *)hp->Clone("hs_phiR");
  hm = (TH1D *)hp->Clone("hm_phiR");
  hALU = (TH1D * )hp->Clone("hALU_phiR");
  hALU -> SetTitle("A_{LU}^{sin(#phi_{R#perp})} = (N^{+} - N^{-}) / (N^{+} + N^{-})");
  hs->Add(hp,hn,1,1);
  hm->Add(hp,hn,1,-1);
  hALU->Divide(hm,hs);
  ff->SetParameter(0,0.01);
  res = hALU->Fit(ff,"Rs+");
  ///////// end phiR asymmetry /////
  
  /// pip phiH asymmetry ///
  ofile->GetObject("hp_phiH",hp);
  ofile->GetObject("hn_phiH",hn);
  hp->Sumw2();
  hn->Sumw2();
  hs = (TH1D *)hp->Clone("hs_phiH");
  hm = (TH1D *)hp->Clone("hm_phiH");
  hALU = (TH1D * )hp->Clone("hALU_phiH");
  hALU -> SetTitle("A_{LU}^{sin(#phi_{H})} = (N^{+} - N^{-}) / (N^{+} + N^{-})");
  hs->Add(hp,hn,1,1);
  hm->Add(hp,hn,1,-1);
  hALU->Divide(hm,hs);
  ff->SetParameter(0,0.01);
  res = hALU->Fit(ff,"Rs+");
  ///////// end pip phiH asymmetry /////
  */

  
  ofile->Write("",TObject::kOverwrite);
  bm->Show("main");

}

Bool_t BSA_survey_cls::DIS()
{

  return (Q2>1)&&(W>2)&&(y<0.85);
}

Bool_t BSA_survey_cls::eFID_ec()
{
  return (20<e_pcal_lv)&&(20<e_pcal_lw);
}

Bool_t BSA_survey_cls::piFID_ec(int k)
{
  return (20<det_pcal_lv[k][0])&&(20<det_pcal_lw[k][0])&&(20<det_pcal_lv[k][1])&&(20<det_pcal_lw[k][1]);
}

Bool_t BSA_survey_cls::eFID_dc()
{
  return true;
}

Bool_t BSA_survey_cls::piFID_dc(int k)
{
  return true;
}

Bool_t BSA_survey_cls::FWD(int k)
{
  return ( ((int)det_statPart[k][0]%4000)/2000 >= 1 ) && ( ((int)det_statPart[k][1]%4000)/2000 >= 1 );
}

Bool_t BSA_survey_cls::CF(int k)
{
  return (vze<20) && (xFm0[k]>0) && (xFm1[k]>0) && (pdata_e[k][0]/Nu>0.2) && (pdata_e[k][1]/Nu>0.2) && ((pdata_e[k][0] + pdata_e[k][1])/Nu<0.95);

}

Float_t BSA_survey_cls::getALU(TString hpname, TString hnname, TString pv, TString tv){
  TH1D *hp,*hn, *hs, *hm, *hALU;
  TFitResultPtr res;
  TF1 *ff = new TF1("ff"," [A]*sin(x*TMath::DegToRad())",0,360);

  ofile->GetObject(hpname,hp);
  ofile->GetObject(hnname,hn);
  hp->Sumw2();
  hn->Sumw2();
  hs = (TH1D *)hp->Clone("hs_" + pv);
  hm = (TH1D *)hp->Clone("hm_" + pv);
  hALU = (TH1D * )hp->Clone("hALU_" + pv);
  hALU -> SetTitle("A_{LU}^{sin(" + tv + ")} = (N^{+} - N^{-}) / (N^{+} + N^{-})");
  hs->Add(hp,hn,1,1);
  hm->Add(hp,hn,1,-1);
  hALU->Divide(hm,hs);
  ff->SetParameter(0,0.01);
  res = hALU->Fit(ff,"Rs+");
  return ff->GetParameter(0);
}
