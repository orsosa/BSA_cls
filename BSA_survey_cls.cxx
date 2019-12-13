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
  std::cout<<"# trees to be processed: "<<fChain->GetNtrees()<<std::endl;
  
  Long64_t nentries = fChain->GetEntries();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry = 0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    TH1D *h;
    ofile->GetObject("hNpair_all",h);
    h->Fill(npart);
    
    //    if (!(npart == 1)) continue;
    if ( !(DIS() && eFID_ec() && eFID_dc() && ePID()) ) continue;
    fillEvHistos();
    //    for (int k = 0; k<npart; k++){
    for (int k = 0; k<1; k++){
      
      if ( !(FWD(k) && CF(k) && piFID_ec(k) && pipFID_dc(k) && pimFID_dc(k) && pi0PID(k)) ) continue; // kin limits.
      fillPartHistos(k);
      
      if (*helicity == -1) 
	ofile->GetObject("hp_phiR",h);
      else if(*helicity == 1)
	ofile->GetObject("hn_phiR",h);
      h->Fill(phiR[k]);
      
      if (*helicity == -1) 
	ofile->GetObject("hp_phiH",h);
      else if(*helicity == 1)
	ofile->GetObject("hn_phiH",h);
      h->Fill(pdata_phiHs[k][0]);
      for (auto& x: bedg){
	TString bn = x.first;
	brv[bn]->GetNdata();
	for (int n=0;n<x.second.size()-1;n++){
	  TString hname = pltv + "_" + bn + Form("_b%d",n);
	  TString hnamepip = "phiH_" + bn + Form("_b%d",n);
	  TString hsin = "hsin" + pltv + "_" + bn;
	  TString hsinpip = "hsinphiH_" + bn;
	  TString ttlsuf =  Form("%.2f<%s<%.2f",x.second[n], bn.Data(), x.second[n+1]);

	  if (x.second[n] < brv[bn]->EvalInstance(k) && brv[bn]->EvalInstance(k)< x.second[n+1]){
	    fillHist(hnamepip,pdata_phiHs[k][0]);
	    fillHist(hname,phiR[k]);
	    fillHist2D(hsin, brv[bn]->EvalInstance(k), sin(phiR[k]*TMath::DegToRad()));
	    fillHist2D(hsinpip, brv[bn]->EvalInstance(k), sin(pdata_phiHs[k][0]*TMath::DegToRad()));

	    break;
	  }
	}
      }
    }
    // if (Cut(ientry) < 0) continue;
    std::cout<<std::setw(15)<<float(jentry+1)/nentries*100<<" %"<<"\r";
    std::cout.flush();
    
  }
 
  TH1D *h;
  Float_t val,err;
  std::cout<<getALU("hp_phiR","hn_phiR",pltv,ttlv,val,err)<<std::endl;
  std::cout<<getALU("hp_phiH","hn_phiH","phiH","#phi_{H}",val,err)<<std::endl;
 
  for (auto& x: bedg){
    TString bn = x.first;
    for (int k=0;k<x.second.size()-1;k++){
      TString hname = pltv + "_" + bn + Form("_b%d",k);
      TString hnamepip = "phiH_" + bn + Form("_b%d",k);
      TString ttlsuf =  Form("%.2f<%s<%.2f",x.second[k], bn.Data(), x.second[k+1]);
      //// ALU pltv ////
      getALU("hp_"+ hname,"hn_"+ hname, hname, ttlv + ", (" + ttlsuf + ")",val,err);
      ofile->GetObject("hALU_"+pltv + "_" + bn,h);
      h->SetBinContent(k+1,val);
      h->SetBinError(k+1,err);
      configHisto(h,bn,"A_{LU}^{sin(" + ttlv + ")}",kBlack);
      ///// end pip ////
      ///// ALU pip ////
      getALU("hp_"+ hnamepip,"hn_"+ hnamepip, hnamepip, "#phi_{H}, (" + ttlsuf  + ")",val,err);
      ofile->GetObject("hALU_phiH_" + bn,h);
      h->SetBinContent(k+1,val);
      h->SetBinError(k+1,err);
      configHisto(h,bn,"A_{LU}^{sin(#phi_{H})}",kBlack);
      //// end ALU pip ////
    }
    TH1D *halu_p, *halu_n, *halu_s, *halu, *halu_pip_p, *halu_pip_n, *halu_pip_s, *halu_pip;
    getALU2D(bn);
    ofile->GetObject("hpALU_"+ pltv + "_" + bn,halu_p);
    ofile->GetObject("hnALU_"+ pltv + "_" + bn,halu_n);
    ofile->GetObject("hsALU_"+ pltv + "_" + bn,halu_s);
    ofile->GetObject("hpALU_phiH_" + bn,halu_pip_p);
    ofile->GetObject("hnALU_phiH_" + bn,halu_pip_n);
    ofile->GetObject("hsALU_phiH_" + bn,halu_pip_s);
    ofile->GetObject("hALU_"+pltv + "_" + bn,halu);
    ofile->GetObject("hALU_phiH_" + bn,halu_pip);
    halu_p->GetYaxis()->SetRangeUser(minALU,maxALU);
    halu_n->GetYaxis()->SetRangeUser(minALU,maxALU);
    halu_s->GetYaxis()->SetRangeUser(minALU,maxALU);
    halu_pip_p->GetYaxis()->SetRangeUser(minALU,maxALU);
    halu_pip_n->GetYaxis()->SetRangeUser(minALU,maxALU);
    halu_pip_s->GetYaxis()->SetRangeUser(minALU,maxALU);
    halu->GetYaxis()->SetRangeUser(minALU,maxALU);
    halu_pip->GetYaxis()->SetRangeUser(minALU,maxALU);
    TCanvas *c = new TCanvas("can_pip_"+ bn + "_all","can_pip_"+bn,800,600);
    halu_pip->Draw();
    halu_pip_p->Draw("same");
    halu_pip_n->Draw("same");
    halu_pip_s->Draw("same");
    ofile->Add(c);
    c = new TCanvas("can_"+ bn + "_all","can_"+bn,800,600);
    halu->Draw();
    halu_p->Draw("same");
    halu_n->Draw("same");
    halu_s->Draw("same");
    ofile->Add(c);
    /// only sum and fit estimations.
    c = new TCanvas("can_pip_"+ bn,"can_pip_"+bn,800,600);
    halu_pip->Draw();
    halu_pip_s->Draw("same");
    ofile->Add(c);
    c = new TCanvas("can_"+ bn,"can_"+bn,800,600);
    halu->Draw();
    halu_s->Draw("same");
    ofile->Add(c);
  
    
  }
  
  ofile->Write("",TObject::kOverwrite);
  ofile->Close();
  bm->Show("main");

}

Bool_t BSA_survey_cls::ePID()
{
  return  (vze<20) && (Pe>2)&&(-2.5<e_chi2pid&&e_chi2pid<2.5);
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
  Bool_t fid1 = true;
  if (options.Contains("pi0")){
    fid1 = (20<det_pcal_lv[k][1])&&(20<det_pcal_lw[k][1])&&(20<det_pcal_lv[k][2])&&(20<det_pcal_lw[k][2]);;
  }
  else{
    fid1 = (20<det_pcal_lv[k][1])&&(20<det_pcal_lw[k][1]);
  }
  
  return (20<det_pcal_lv[k][0])&&(20<det_pcal_lw[k][0])&&fid1;
}

Bool_t BSA_survey_cls::eFID_dc()
{
  return true;
}

Bool_t BSA_survey_cls::pipFID_dc(int k)
{
  return true;
}

Bool_t BSA_survey_cls::pimFID_dc(int k)
{
  return true;
}

Bool_t BSA_survey_cls::FWD(int k)
{
  bool fwd1 = true;
  if (options.Contains("pi0")){
    fwd1 = ( ((int)det_statPart[k][1]%4000)/2000 >= 1 )&&( ((int)det_statPart[k][2]%4000)/2000 >= 1 );
  }
  else{
    fwd1 = ( ((int)det_statPart[k][1]%4000)/2000 >= 1 );
  }
  
  return ( ((int)det_statPart[k][0]%4000)/2000 >= 1 ) && fwd1;
}

Bool_t BSA_survey_cls::CF(int k)
{
  Float_t xF_1 = -1111;
  Float_t E_1 = -1111;
  if (options.Contains("pi0")){
    xF_1 = getxF(pdata_e[k][1] + pdata_e[k][2], pdata_px[k][1] + pdata_px[k][2], pdata_py[k][1] + pdata_py[k][2], pdata_pz[k][1] + pdata_pz[k][2]);
    E_1 = (pdata_e[k][1] + pdata_e[k][2]);
  }
  else{
    xF_1 = xFm1[k];
    E_1 = pdata_e[k][1];
  }
  
  
  return (xFm0[k]>0) && (xF_1>0) && (pdata_e[k][0]/Nu>0.1) && (E_1/Nu>0.1) && ((pdata_e[k][0] + E_1)/Nu<0.95);

}

Float_t BSA_survey_cls::getALU(TString hpname, TString hnname, TString pv, TString tv, Float_t &val, Float_t &err){
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
  std::cout<<tv<<std::endl;
  res = hALU->Fit(ff,"Rs+0");
  val = ff->GetParameter(0);
  err = ff->GetParError(0);
  return ff->GetParameter(0);
}


Float_t BSA_survey_cls::getALU2D(TString bn){
  TString psname, nsname;
  TH2D *hp,*hn;
  TH1D *halup,*halun,*halus;

  ////  pltv ////
  psname = "hsin"+pltv + "_" + bn + "_p";
  nsname = "hsin"+pltv + "_" + bn + "_n";
  ofile->GetObject(psname,hp);
  ofile->GetObject(nsname,hn);

  hp->Sumw2();
  hn->Sumw2();

  halup = (TH1D*)hp->ProfileX("hpALU_"+ pltv + "_" + bn);
  halun = (TH1D*)hn->ProfileX("hnALU_"+ pltv + "_" + bn);
  halus = (TH1D*)halup->Clone("hsALU_"+ pltv + "_" + bn);
  halup->Scale(2.);
  halun->Scale(2.);
  halus->Add(halup,halun,0.5,0.5);

  halup->SetTitle("2<sin(" + ttlv + ")>^{+}");
  halun->SetTitle("-2<sin(" + ttlv + ")>^{-}");
  halus->SetTitle("(<sin(" + ttlv + ")>^{+} - <sin(" + ttlv+ ")>^{-})");
  configHisto(halup,bn,"A_{LU}(<sin(" + ttlv + ")>^{+})",kRed,kFullTriangleUp);
  configHisto(halun,bn,"A_{LU}(<sin(" + ttlv + ")>^{-})",kGreen+3,kFullTriangleDown);
  configHisto(halus,bn,"A_{LU}(#sum<sin(" + ttlv + ")>^{+/-})",kMagenta+1,kFullStar);
    
  ///// end pltv /////

  ////  phiH ////
  psname = "hsinphiH_" + bn + "_p";
  nsname = "hsinphiH_" + bn + "_n";
  ofile->GetObject(psname,hp);
  ofile->GetObject(nsname,hn);

  hp->Sumw2();
  hn->Sumw2();

  halup = (TH1D*)hp->ProfileX("hpALU_phiH_" + bn);
  halun = (TH1D*)hn->ProfileX("hnALU_phiH_" + bn);
  halus = (TH1D*)halup->Clone("hsALU_phiH_" + bn);
  halup->Scale(2.);
  halun->Scale(2.);
  
  halus->Add(halup,halun,0.5,0.5);

  halup->SetTitle("2<sin(#phi_{H})>^{+}");
  halun->SetTitle("-2<sin(#phi_{H})>^{-}");
  halus->SetTitle("(<sin(#phi_{H})>^{+} - <sin(#phi_{H})>^{-})");

  configHisto(halup,bn,"A_{LU}(<sin(#phi_{H})>^{+})",kRed,kFullTriangleUp);
  configHisto(halun,bn,"A_{LU}(<sin(#phi_{H})>^{-})",kGreen+3,kFullTriangleDown);
  configHisto(halus,bn,"A_{LU}(#sum<sin(#phi_{H})>^{+/-})",kMagenta+1,kFullStar);
  
  
  ///// end phiH /////

  return 0;
}



Int_t BSA_survey_cls::fillHist(TString hname, Float_t value){
  TH1D *h;
  if (*helicity == -1) // positive helicity, it is flipped!
    ofile->GetObject("hp_"+ hname,h);
  else if(*helicity == 1)// negative helicity, it is flipped!
    ofile->GetObject("hn_"+ hname,h);
  else
    return 1;
  h->Fill(value);
  return 0;
}

Int_t BSA_survey_cls::fillHist2D(TString hname, Float_t x, Float_t y){
  TH2D *h;
  if (*helicity == -1) { // positive helicity, it is flipped!
    ofile->GetObject(hname + "_p",h);
    h->Fill(x,y);
  }
  else if(*helicity == 1){ // negative helicity, it is flipped!
    ofile->GetObject(hname + "_n",h);
    h->Fill(x,-y);
  }
  else
    return 1;
  
  return 0;
}

Int_t BSA_survey_cls::configHisto(TH1D *h, TString xtitle, TString ytitle, Color_t c, EMarkerStyle ms){

  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);
  h->SetLineColor(c);
  h->SetMarkerColor(c);
  h->SetMarkerStyle(ms);
  h->SetMarkerSize(1.3);
  
  return 0;
  
}

Int_t BSA_survey_cls::LoadElecFIDPar(){
  std::ifstream fpar("/home/orsosa/Dropbox/INFN_work/Phys_ana/PID/DC_elec_par.txt");

  char junk[100];
  fpar.getline(junk,100);
  fpar.getline(junk,100);
  fpar.getline(junk,100);
  fpar.getline(junk,100);
  fpar.getline(junk,100);
  fpar.getline(junk,100);

  for (int k=0;k<NSECTORS;k++)
  {
    fpar>>pl0_e[k]>>pl1_e[k]>>pl2_e[k]>>pl3_e[k]>>pr0_e[k]>>pr1_e[k]>>pr2_e[k]>>pr3_e[k];
  }
  fpar.close();
  return 0;
}


Int_t BSA_survey_cls::setStyle(){

  myStyle  = new TStyle("orsosaStyle","My Root Styles");
  myStyle->SetHistMinimumZero(0);
  myStyle->SetPalette(1,0);
  //myStyle->SetPalette(kBlueYellow);
  myStyle->SetCanvasBorderMode(0);
  myStyle->SetPadBorderMode(0);
  myStyle->SetPadColor(0);
  myStyle->SetCanvasColor(0);
  myStyle->SetTitleFillColor(0);
  myStyle->SetTitleBorderSize(0);

  myStyle->SetStatColor(0);
  myStyle->SetOptStat("e");

  myStyle->SetLabelSize(0.05,"xyz"); // size of axis value font
  myStyle->SetTitleSize(0.05,"xyz"); // size of axis title font
  myStyle->SetTitleOffset(0.75,"xyz"); // axis title offset 
  myStyle->SetTitleFont(22,"xyz"); // font option
  myStyle->SetTitleFont(22,"a"); // pad font option
  myStyle->SetLabelFont(22,"xyz");


  myStyle->SetLabelSize(0.02,"z"); // size of axis value font
  myStyle->SetLabelOffset(-0.03,"z"); // size of axis value font
  myStyle->SetTickLength(0.002,"z"); // size of axis value font

  
  // Stat and legend fonts
  myStyle->SetStatFont(22); 
  myStyle->SetLegendFont(22); 
  // hiostogram style
  myStyle->SetHistLineWidth(2);
  myStyle->SetCanvasDefH(768);
  myStyle->SetCanvasDefW(1024);

  myStyle->SetPadBottomMargin(0.1);
  myStyle->SetPadTopMargin(0.1);
  myStyle->SetPadLeftMargin(0.1);
  myStyle->SetPadRightMargin(0.075);

  myStyle->SetPadTickX(1);
  myStyle->SetPadTickY(1);

  myStyle->SetFrameBorderMode(0);
  
  myStyle->SetGridStyle(3);
  myStyle->SetGridWidth(2);
  myStyle->SetPadGridX(kTRUE);
  myStyle->SetPadGridY(kTRUE);

  
  gROOT->SetStyle("orsosaStyle"); //uncomment to set this style
 
  
  return 0;
}

Bool_t BSA_survey_cls::pi0PID(Int_t k){
  if (!options.Contains("pi0")) return kTRUE;
  Float_t minth_brem = 10;
  Float_t minE = 0.5;
  Bool_t ret = true;
  ret = ret && acos((pdata_px[k][1]*Pex + pdata_py[k][1]*Pey + pdata_pz[k][1]*Pez)/pdata_e[k][1]/(Nu/y-Nu))*TMath::RadToDeg()>minth_brem;
  ret = ret && acos((pdata_px[k][2]*Pex + pdata_py[k][2]*Pey + pdata_pz[k][2]*Pez)/pdata_e[k][2]/(Nu/y-Nu))*TMath::RadToDeg()>minth_brem; 
  ret = ret && (pdata_e[k][1]>minE)&&(pdata_e[k][2]>minE);
  Float_t mpi0 = sqrt(2*(pdata_e[k][1]*pdata_e[k][2] - pdata_px[k][1]*pdata_px[k][2] - pdata_py[k][1]*pdata_py[k][2] - pdata_pz[k][1]*pdata_pz[k][2]));
  ret = ret &&  0.104<mpi0&&mpi0<0.167;
  return ret;
}


Float_t BSA_survey_cls::getxF(Float_t E, Float_t Px, Float_t Py, Float_t Pz){
  Float_t kMprt = 0.93827;
  Float_t kEbeam = Nu/y;
  Float_t P2 = Px*Px + Py*Py + Pz*Pz;
  Float_t cospq = ((kEbeam-Pez)*Pz - Pex*Px - Pey*Py)/( sqrt((Q2 + Nu*Nu)*P2) );
  Float_t Pt2 = P2*(1-cospq*cospq);
  Float_t Pl2 = P2*cospq*cospq;
  Float_t Pl = sqrt(P2)*cospq;
  ////// LORENTZ BOOST //////////
  Float_t b=TMath::Sqrt(Q2 + Nu*Nu)/(Nu + kMprt);
  Float_t g=(Nu + kMprt)/W;

  Float_t PlCM = g*(Pl - b*E);
      
  Float_t xFm = 2*PlCM/W;
  return xFm;

}
