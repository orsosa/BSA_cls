//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Nov 18 15:49:54 2019 by ROOT version 6.18/02
// from TTree outdata/:pi+_pi-, all combinations
// found on file: ../filtering/outfiles/pippim_sim_ga_1.root
//////////////////////////////////////////////////////////

#ifndef BSA_survey_cls_h
#define BSA_survey_cls_h

#include <TROOT.h>
#include <TChain.h>
#include <TH1F.h>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCut.h"
#include "TROOT.h"
#include "TBenchmark.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTreeFormula.h"

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include "TObject.h"

class BSA_survey_cls {
 public :
  TBenchmark *bm;
  TString OUTDIR;
  TFile *ofile;
  TChain          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  TStyle *myStyle;
  // Fixed size dimensions of array or collections stored in the TTree if any.
  static constexpr Int_t kMaxdet = 50;
  static constexpr Int_t kMaxpdata = 50;
  static constexpr Int_t kMaxmc_pdata = 50;
  static constexpr Float_t minALU = -0.15;
  static constexpr Float_t maxALU = 0.15;

  //// bin variables
  std::vector <Double_t> *be = 0;
  TString pltv = "", ttlv = "";  
  std::map<TString, std::vector <Double_t>> bedg;
  std::map<TString, TTreeFormula *> brv;

  
  // Declaration of leaf types
  Int_t           npart;
  Float_t         M[kMaxpdata];   //[npart]
  Float_t         Phx[kMaxpdata];   //[npart]
  Float_t         Phy[kMaxpdata];   //[npart]
  Float_t         Phz[kMaxpdata];   //[npart]
  Float_t         Z[kMaxpdata];   //[npart]
  Float_t         Cospq[kMaxpdata];   //[npart]
  Float_t         Pt2[kMaxpdata];   //[npart]
  Float_t         Event[kMaxpdata];   //[npart]
  Float_t         M2_01[kMaxpdata];   //[npart]
  Float_t         M2_02[kMaxpdata];   //[npart]
  Float_t         phiH[kMaxpdata];   //[npart]
  Float_t         phiR[kMaxpdata];   //[npart]
  Float_t         helic002_phiR[kMaxpdata];   //[npart]
  Float_t         helic005_phiR[kMaxpdata];   //[npart]
  Float_t         helic010_phiR[kMaxpdata];   //[npart]
  Float_t         helic020_phiR[kMaxpdata];   //[npart]
  Float_t         helic002_phiRst[kMaxpdata];   //[npart]
  Float_t         helic005_phiRst[kMaxpdata];   //[npart]
  Float_t         helic010_phiRst[kMaxpdata];   //[npart]
  Float_t         helic020_phiRst[kMaxpdata];   //[npart]
  Float_t         wUxS_phiR[kMaxpdata];   //[npart]
  Float_t         R[kMaxpdata];   //[npart]
  Float_t         KF[kMaxpdata];   //[npart]
  Float_t         Mx2[kMaxpdata];   //[npart]
  Float_t         xF[kMaxpdata];   //[npart]
  Float_t         xF0[kMaxpdata];   //[npart]
  Float_t         xF1[kMaxpdata];   //[npart]
  Float_t         plcm[kMaxpdata];   //[npart]
  Float_t         plcm0[kMaxpdata];   //[npart]
  Float_t         plcm1[kMaxpdata];   //[npart]
  Float_t         Eh[kMaxpdata];   //[npart]
  Float_t         xFm[kMaxpdata];   //[npart]
  Float_t         xFm0[kMaxpdata];   //[npart]
  Float_t         xFm1[kMaxpdata];   //[npart]
  Float_t         theta0[kMaxpdata];   //[npart]
  Float_t         theta1[kMaxpdata];   //[npart]
  Float_t         cos_theta_P0cm[kMaxpdata];   //[npart]
  Float_t         sin_theta_P0cm[kMaxpdata];   //[npart]
  Float_t         xFo[kMaxpdata];   //[npart]
  Float_t         xFo0[kMaxpdata];   //[npart]
  Float_t         xFo1[kMaxpdata];   //[npart]
  Float_t         phiH_phiR[kMaxpdata];   //[npart]
  Float_t         phiR_cov[kMaxpdata];   //[npart]
  Float_t         p0T2[kMaxpdata];   //[npart]
  Float_t         p1T2[kMaxpdata];   //[npart]
  Float_t         phipq[kMaxpdata];   //[npart]
  Float_t         phT2[kMaxpdata];   //[npart]
  Float_t         etaCM0[kMaxpdata];   //[npart]
  Float_t         etaCM1[kMaxpdata];   //[npart]
  Float_t         etaBF0p[kMaxpdata];   //[npart]
  Float_t         etaBF1p[kMaxpdata];   //[npart]
  Float_t         etaBF0m[kMaxpdata];   //[npart]
  Float_t         etaBF1m[kMaxpdata];   //[npart]
  Float_t         etaBF0[kMaxpdata];   //[npart]
  Float_t         etaBF1[kMaxpdata];   //[npart]
  Float_t         phiR_ha[kMaxpdata];   //[npart]
  Float_t         plcm0_r[kMaxpdata];   //[npart]
  Float_t         plcm1_r[kMaxpdata];   //[npart]
  Float_t         phiR_covH[kMaxpdata];   //[npart]
  Float_t         E0_phcm[kMaxpdata];   //[npart]
  Float_t         E1_phcm[kMaxpdata];   //[npart]
  Int_t           mix_npart;
  Int_t           mc_mix_npart;
  Int_t           det_;
  UInt_t          det_fUniqueID[kMaxdet];   //[det_]
  UInt_t          det_fBits[kMaxdet];   //[det_]
  Float_t         det_beta[kMaxdet][6];   //[det_]
  Float_t         det_m2b[kMaxdet][6];   //[det_]
  Float_t         det_vx[kMaxdet][6];   //[det_]
  Float_t         det_vy[kMaxdet][6];   //[det_]
  Float_t         det_vz[kMaxdet][6];   //[det_]
  Float_t         det_dcx[kMaxdet][6];   //[det_]
  Float_t         det_dcy[kMaxdet][6];   //[det_]
  Float_t         det_dcz[kMaxdet][6];   //[det_]
  Float_t         det_statPart[kMaxdet][6];   //[det_]
  Float_t         det_dc_chi2[kMaxdet][6];   //[det_]
  Float_t         det_dc_ndf[kMaxdet][6];   //[det_]
  Float_t         det_pcal_lu[kMaxdet][6];   //[det_]
  Float_t         det_pcal_lv[kMaxdet][6];   //[det_]
  Float_t         det_pcal_lw[kMaxdet][6];   //[det_]
  Int_t           pdata_;
  UInt_t          pdata_fUniqueID[kMaxpdata];   //[pdata_]
  UInt_t          pdata_fBits[kMaxpdata];   //[pdata_]
  Float_t         pdata_e[kMaxpdata][6];   //[pdata_]
  Float_t         pdata_px[kMaxpdata][6];   //[pdata_]
  Float_t         pdata_py[kMaxpdata][6];   //[pdata_]
  Float_t         pdata_pz[kMaxpdata][6];   //[pdata_]
  Float_t         pdata_helic002_phiH[kMaxpdata][6];   //[pdata_]
  Float_t         pdata_helic005_phiH[kMaxpdata][6];   //[pdata_]
  Float_t         pdata_helic010_phiH[kMaxpdata][6];   //[pdata_]
  Float_t         pdata_helic020_phiH[kMaxpdata][6];   //[pdata_]
  Float_t         pdata_phiHs[kMaxpdata][6];   //[pdata_]
  Int_t           mc_pdata_;
  UInt_t          mc_pdata_fUniqueID[kMaxmc_pdata];   //[mc_pdata_]
  UInt_t          mc_pdata_fBits[kMaxmc_pdata];   //[mc_pdata_]
  Float_t         mc_pdata_e[kMaxmc_pdata][6];   //[mc_pdata_]
  Float_t         mc_pdata_px[kMaxmc_pdata][6];   //[mc_pdata_]
  Float_t         mc_pdata_py[kMaxmc_pdata][6];   //[mc_pdata_]
  Float_t         mc_pdata_pz[kMaxmc_pdata][6];   //[mc_pdata_]
  Float_t         mc_pdata_helic002_phiH[kMaxmc_pdata][6];   //[mc_pdata_]
  Float_t         mc_pdata_helic005_phiH[kMaxmc_pdata][6];   //[mc_pdata_]
  Float_t         mc_pdata_helic010_phiH[kMaxmc_pdata][6];   //[mc_pdata_]
  Float_t         mc_pdata_helic020_phiH[kMaxmc_pdata][6];   //[mc_pdata_]
  Float_t         mc_pdata_phiHs[kMaxmc_pdata][6];   //[mc_pdata_]
  Int_t           mc_npart;
  Float_t         mc_M[kMaxpdata];   //[mc_npart]
  Float_t         mc_Phx[kMaxpdata];   //[mc_npart]
  Float_t         mc_Phy[kMaxpdata];   //[mc_npart]
  Float_t         mc_Phz[kMaxpdata];   //[mc_npart]
  Float_t         mc_Z[kMaxpdata];   //[mc_npart]
  Float_t         mc_Cospq[kMaxpdata];   //[mc_npart]
  Float_t         mc_Pt2[kMaxpdata];   //[mc_npart]
  Float_t         mc_Event[kMaxpdata];   //[mc_npart]
  Float_t         mc_M2_01[kMaxpdata];   //[mc_npart]
  Float_t         mc_M2_02[kMaxpdata];   //[mc_npart]
  Float_t         mc_phiH[kMaxpdata];   //[mc_npart]
  Float_t         mc_phiR[kMaxpdata];   //[mc_npart]
  Float_t         mc_helic002_phiR[kMaxpdata];   //[mc_npart]
  Float_t         mc_helic005_phiR[kMaxpdata];   //[mc_npart]
  Float_t         mc_helic010_phiR[kMaxpdata];   //[mc_npart]
  Float_t         mc_helic020_phiR[kMaxpdata];   //[mc_npart]
  Float_t         mc_helic002_phiRst[kMaxpdata];   //[mc_npart]
  Float_t         mc_helic005_phiRst[kMaxpdata];   //[mc_npart]
  Float_t         mc_helic010_phiRst[kMaxpdata];   //[mc_npart]
  Float_t         mc_helic020_phiRst[kMaxpdata];   //[mc_npart]
  Float_t         mc_wUxS_phiR[kMaxpdata];   //[mc_npart]
  Float_t         mc_Mx2[kMaxpdata];   //[mc_npart]
  Float_t         mc_xF[kMaxpdata];   //[mc_npart]
  Float_t         mc_xF0[kMaxpdata];   //[mc_npart]
  Float_t         mc_xF1[kMaxpdata];   //[mc_npart]
  Float_t         mc_plcm[kMaxpdata];   //[mc_npart]
  Float_t         mc_plcm0[kMaxpdata];   //[mc_npart]
  Float_t         mc_plcm1[kMaxpdata];   //[mc_npart]
  Float_t         mc_Eh[kMaxpdata];   //[mc_npart]
  Float_t         mc_xFm[kMaxpdata];   //[mc_npart]
  Float_t         mc_xFm0[kMaxpdata];   //[mc_npart]
  Float_t         mc_xFm1[kMaxpdata];   //[mc_npart]
  Float_t         mc_theta0[kMaxpdata];   //[mc_npart]
  Float_t         mc_theta1[kMaxpdata];   //[mc_npart]
  Float_t         mc_cos_theta_P0cm[kMaxpdata];   //[mc_npart]
  Float_t         mc_sin_theta_P0cm[kMaxpdata];   //[mc_npart]
  Float_t         mc_xFo[kMaxpdata];   //[mc_npart]
  Float_t         mc_xFo0[kMaxpdata];   //[mc_npart]
  Float_t         mc_xFo1[kMaxpdata];   //[mc_npart]
  Float_t         mc_phiH_phiR[kMaxpdata];   //[mc_npart]
  Float_t         mc_phiR_cov[kMaxpdata];   //[mc_npart]
  Float_t         mc_p0T2[kMaxpdata];   //[mc_npart]
  Float_t         mc_p1T2[kMaxpdata];   //[mc_npart]
  Float_t         mc_phipq[kMaxpdata];   //[mc_npart]
  Float_t         mc_phT2[kMaxpdata];   //[mc_npart]
  Float_t         mc_etaCM0[kMaxpdata];   //[mc_npart]
  Float_t         mc_etaCM1[kMaxpdata];   //[mc_npart]
  Float_t         mc_etaBF0p[kMaxpdata];   //[mc_npart]
  Float_t         mc_etaBF1p[kMaxpdata];   //[mc_npart]
  Float_t         mc_etaBF0m[kMaxpdata];   //[mc_npart]
  Float_t         mc_etaBF1m[kMaxpdata];   //[mc_npart]
  Float_t         mc_etaBF0[kMaxpdata];   //[mc_npart]
  Float_t         mc_etaBF1[kMaxpdata];   //[mc_npart]
  Float_t         mc_phiR_ha[kMaxpdata];   //[mc_npart]
  Float_t         mc_plcm0_r[kMaxpdata];   //[mc_npart]
  Float_t         mc_plcm1_r[kMaxpdata];   //[mc_npart]
  Float_t         mc_phiR_covH[kMaxpdata];   //[mc_npart]
  Float_t         mc_E0_phcm[kMaxpdata];   //[mc_npart]
  Float_t         mc_E1_phcm[kMaxpdata];   //[mc_npart]
  Float_t         mc_R[kMaxpdata];   //[mc_npart]
  Float_t         mc_KF[kMaxpdata];   //[mc_npart]
  Float_t         Q2;
  Float_t         W;
  Float_t         Nu;
  Float_t         Xb;
  Float_t         vxec;
  Float_t         vyec;
  Float_t         vzec;
  Float_t         vxe;
  Float_t         vye;
  Float_t         vze;
  Float_t         Pex;
  Float_t         Pey;
  Float_t         Pez;
  Float_t         event;
  Float_t         Pe;
  Float_t         Ee;
  Float_t         e_Ein;
  Float_t         e_Eout;
  Float_t         e_Epcal;
  Float_t         e_npheltcc;
  Float_t         e_nphehtcc;
  Float_t         helic;
  Float_t         e_chi2pid;
  Float_t         e_pcal_lu;
  Float_t         e_pcal_lv;
  Float_t         e_pcal_lw;
  Float_t         e_ecin_lu;
  Float_t         e_ecin_lv;
  Float_t         e_ecin_lw;
  Float_t         e_ecout_lu;
  Float_t         e_ecout_lv;
  Float_t         e_ecout_lw;
  Float_t         e_pcal_hx;
  Float_t         e_pcal_hy;
  Float_t         e_pcal_hz;
  Float_t         e_ecin_hx;
  Float_t         e_ecin_hy;
  Float_t         e_ecin_hz;
  Float_t         e_ecout_hx;
  Float_t         e_ecout_hy;
  Float_t         e_ecout_hz;
  Float_t         e_trajx_sl0;
  Float_t         e_trajx_sl1;
  Float_t         e_trajx_sl2;
  Float_t         e_trajx_sl3;
  Float_t         e_trajx_sl4;
  Float_t         e_trajx_sl5;
  Float_t         e_trajy_sl0;
  Float_t         e_trajy_sl1;
  Float_t         e_trajy_sl2;
  Float_t         e_trajy_sl3;
  Float_t         e_trajy_sl4;
  Float_t         e_trajy_sl5;
  Float_t         e_trajz_sl0;
  Float_t         e_trajz_sl1;
  Float_t         e_trajz_sl2;
  Float_t         e_trajz_sl3;
  Float_t         e_trajz_sl4;
  Float_t         e_trajz_sl5;
  Float_t         e_pathtof;
  Float_t         e_timetof;
  Float_t         e_sector_tof;
  Float_t         e_Beta;
  Float_t         STTime;
  Float_t         RFTime;
  Float_t         e_dcx_rot_0;
  Float_t         e_dcy_rot_0;
  Float_t         e_dcx_rot_1;
  Float_t         e_dcy_rot_1;
  Float_t         e_dcx_rot_2;
  Float_t         e_dcy_rot_2;
  Float_t         e_sector_ltcc;
  Float_t         e_sector_htcc;
  Float_t         e_sector_ecal;
  Float_t         revent;
  Float_t         e_dc_chi2;
  Float_t         e_ftof1ax;
  Float_t         e_ftof1ay;
  Float_t         e_ftof1az;
  Float_t         e_pcalx;
  Float_t         e_pcaly;
  Float_t         e_pcalz;
  Float_t         e_ecalx;
  Float_t         e_ecaly;
  Float_t         e_ecalz;
  Float_t         e_ltccx;
  Float_t         e_ltccy;
  Float_t         e_ltccz;
  Float_t         e_htccx;
  Float_t         e_htccy;
  Float_t         e_htccz;
  Float_t         e_ftof1bx;
  Float_t         e_ftof1by;
  Float_t         e_ftof1bz;
  Float_t         e_ftof2x;
  Float_t         e_ftof2y;
  Float_t         e_ftof2z;
  Float_t         helonline_hel;
  Float_t         helonline_helRaw;
  Float_t         helflip_hel;
  Float_t         helflip_helRaw;
  Float_t         helflip_event;
  Float_t         e_dc_status;
  Float_t         e_dc_ndf;
  Float_t         e_sector_dc;
  Float_t         e_statPart;
  Float_t         e_DCPx;
  Float_t         e_DCPy;
  Float_t         e_DCPz;
  Float_t         y;
  Float_t         th_e;
  Float_t         phi_e;
  Float_t         helicRaw;
  Float_t         fA;
  Float_t         fB;
  Float_t         fC;
  Float_t         fV;
  Float_t         fW;
  Float_t         mc_Q2;
  Float_t         mc_W;
  Float_t         mc_Nu;
  Float_t         mc_Xb;
  Float_t         mc_vxe;
  Float_t         mc_vye;
  Float_t         mc_vze;
  Float_t         mc_Pex;
  Float_t         mc_Pey;
  Float_t         mc_Pez;
  Float_t         mc_event;
  Float_t         e_mcmass;
  Float_t         mc_Pe;
  Float_t         mc_Ee;
  Float_t         mc_revent;
  Float_t         mc_y;
  Float_t         mc_th_e;
  Float_t         mc_phi_e;
  Float_t         mc_e_Beta;
  Float_t         mc_helic;
  Float_t         mc_fA;
  Float_t         mc_fB;
  Float_t         mc_fC;
  Float_t         mc_fV;
  Float_t         mc_fW;

  // List of branches
  TBranch        *b_npart;   //!
  TBranch        *b_M;   //!
  TBranch        *b_Phx;   //!
  TBranch        *b_Phy;   //!
  TBranch        *b_Phz;   //!
  TBranch        *b_Z;   //!
  TBranch        *b_Cospq;   //!
  TBranch        *b_Pt2;   //!
  TBranch        *b_Event;   //!
  TBranch        *b_M2_01;   //!
  TBranch        *b_M2_02;   //!
  TBranch        *b_phiH;   //!
  TBranch        *b_phiR;   //!
  TBranch        *b_helic002_phiR;   //!
  TBranch        *b_helic005_phiR;   //!
  TBranch        *b_helic010_phiR;   //!
  TBranch        *b_helic020_phiR;   //!
  TBranch        *b_helic002_phiRst;   //!
  TBranch        *b_helic005_phiRst;   //!
  TBranch        *b_helic010_phiRst;   //!
  TBranch        *b_helic020_phiRst;   //!
  TBranch        *b_wUxS_phiR;   //!
  TBranch        *b_R;   //!
  TBranch        *b_KF;   //!
  TBranch        *b_Mx2;   //!
  TBranch        *b_xF;   //!
  TBranch        *b_xF0;   //!
  TBranch        *b_xF1;   //!
  TBranch        *b_plcm;   //!
  TBranch        *b_plcm0;   //!
  TBranch        *b_plcm1;   //!
  TBranch        *b_Eh;   //!
  TBranch        *b_xFm;   //!
  TBranch        *b_xFm0;   //!
  TBranch        *b_xFm1;   //!
  TBranch        *b_theta0;   //!
  TBranch        *b_theta1;   //!
  TBranch        *b_cos_theta_P0cm;   //!
  TBranch        *b_sin_theta_P0cm;   //!
  TBranch        *b_xFo;   //!
  TBranch        *b_xFo0;   //!
  TBranch        *b_xFo1;   //!
  TBranch        *b_phiH_phiR;   //!
  TBranch        *b_phiR_cov;   //!
  TBranch        *b_p0T2;   //!
  TBranch        *b_p1T2;   //!
  TBranch        *b_phipq;   //!
  TBranch        *b_phT2;   //!
  TBranch        *b_etaCM0;   //!
  TBranch        *b_etaCM1;   //!
  TBranch        *b_etaBF0p;   //!
  TBranch        *b_etaBF1p;   //!
  TBranch        *b_etaBF0m;   //!
  TBranch        *b_etaBF1m;   //!
  TBranch        *b_etaBF0;   //!
  TBranch        *b_etaBF1;   //!
  TBranch        *b_phiR_ha;   //!
  TBranch        *b_plcm0_r;   //!
  TBranch        *b_plcm1_r;   //!
  TBranch        *b_phiR_covH;   //!
  TBranch        *b_E0_phcm;   //!
  TBranch        *b_E1_phcm;   //!
  TBranch        *b_mix_npart;   //!
  TBranch        *b_mc_mix_npart;   //!
  TBranch        *b_det_;   //!
  TBranch        *b_det_fUniqueID;   //!
  TBranch        *b_det_fBits;   //!
  TBranch        *b_det_beta;   //!
  TBranch        *b_det_m2b;   //!
  TBranch        *b_det_vx;   //!
  TBranch        *b_det_vy;   //!
  TBranch        *b_det_vz;   //!
  TBranch        *b_det_dcx;   //!
  TBranch        *b_det_dcy;   //!
  TBranch        *b_det_dcz;   //!
  TBranch        *b_det_statPart;   //!
  TBranch        *b_det_dc_chi2;   //!
  TBranch        *b_det_dc_ndf;   //!
  TBranch        *b_det_pcal_lu;   //!
  TBranch        *b_det_pcal_lv;   //!
  TBranch        *b_det_pcal_lw;   //!
  TBranch        *b_pdata_;   //!
  TBranch        *b_pdata_fUniqueID;   //!
  TBranch        *b_pdata_fBits;   //!
  TBranch        *b_pdata_e;   //!
  TBranch        *b_pdata_px;   //!
  TBranch        *b_pdata_py;   //!
  TBranch        *b_pdata_pz;   //!
  TBranch        *b_pdata_helic002_phiH;   //!
  TBranch        *b_pdata_helic005_phiH;   //!
  TBranch        *b_pdata_helic010_phiH;   //!
  TBranch        *b_pdata_helic020_phiH;   //!
  TBranch        *b_pdata_phiHs;   //!
  TBranch        *b_mc_pdata_;   //!
  TBranch        *b_mc_pdata_fUniqueID;   //!
  TBranch        *b_mc_pdata_fBits;   //!
  TBranch        *b_mc_pdata_e;   //!
  TBranch        *b_mc_pdata_px;   //!
  TBranch        *b_mc_pdata_py;   //!
  TBranch        *b_mc_pdata_pz;   //!
  TBranch        *b_mc_pdata_helic002_phiH;   //!
  TBranch        *b_mc_pdata_helic005_phiH;   //!
  TBranch        *b_mc_pdata_helic010_phiH;   //!
  TBranch        *b_mc_pdata_helic020_phiH;   //!
  TBranch        *b_mc_pdata_phiHs;   //!
  TBranch        *b_mc_npart;   //!
  TBranch        *b_mc_M;   //!
  TBranch        *b_mc_Phx;   //!
  TBranch        *b_mc_Phy;   //!
  TBranch        *b_mc_Phz;   //!
  TBranch        *b_mc_Z;   //!
  TBranch        *b_mc_Cospq;   //!
  TBranch        *b_mc_Pt2;   //!
  TBranch        *b_mc_Event;   //!
  TBranch        *b_mc_M2_01;   //!
  TBranch        *b_mc_M2_02;   //!
  TBranch        *b_mc_phiH;   //!
  TBranch        *b_mc_phiR;   //!
  TBranch        *b_mc_helic002_phiR;   //!
  TBranch        *b_mc_helic005_phiR;   //!
  TBranch        *b_mc_helic010_phiR;   //!
  TBranch        *b_mc_helic020_phiR;   //!
  TBranch        *b_mc_helic002_phiRst;   //!
  TBranch        *b_mc_helic005_phiRst;   //!
  TBranch        *b_mc_helic010_phiRst;   //!
  TBranch        *b_mc_helic020_phiRst;   //!
  TBranch        *b_mc_wUxS_phiR;   //!
  TBranch        *b_mc_Mx2;   //!
  TBranch        *b_mc_xF;   //!
  TBranch        *b_mc_xF0;   //!
  TBranch        *b_mc_xF1;   //!
  TBranch        *b_mc_plcm;   //!
  TBranch        *b_mc_plcm0;   //!
  TBranch        *b_mc_plcm1;   //!
  TBranch        *b_mc_Eh;   //!
  TBranch        *b_mc_xFm;   //!
  TBranch        *b_mc_xFm0;   //!
  TBranch        *b_mc_xFm1;   //!
  TBranch        *b_mc_theta0;   //!
  TBranch        *b_mc_theta1;   //!
  TBranch        *b_mc_cos_theta_P0cm;   //!
  TBranch        *b_mc_sin_theta_P0cm;   //!
  TBranch        *b_mc_xFo;   //!
  TBranch        *b_mc_xFo0;   //!
  TBranch        *b_mc_xFo1;   //!
  TBranch        *b_mc_phiH_phiR;   //!
  TBranch        *b_mc_phiR_cov;   //!
  TBranch        *b_mc_p0T2;   //!
  TBranch        *b_mc_p1T2;   //!
  TBranch        *b_mc_phipq;   //!
  TBranch        *b_mc_phT2;   //!
  TBranch        *b_mc_etaCM0;   //!
  TBranch        *b_mc_etaCM1;   //!
  TBranch        *b_mc_etaBF0p;   //!
  TBranch        *b_mc_etaBF1p;   //!
  TBranch        *b_mc_etaBF0m;   //!
  TBranch        *b_mc_etaBF1m;   //!
  TBranch        *b_mc_etaBF0;   //!
  TBranch        *b_mc_etaBF1;   //!
  TBranch        *b_mc_phiR_ha;   //!
  TBranch        *b_mc_plcm0_r;   //!
  TBranch        *b_mc_plcm1_r;   //!
  TBranch        *b_mc_phiR_covH;   //!
  TBranch        *b_mc_E0_phcm;   //!
  TBranch        *b_mc_E1_phcm;   //!
  TBranch        *b_mc_R;   //!
  TBranch        *b_mc_KF;   //!
  TBranch        *b_Q2;   //!
  TBranch        *b_W;   //!
  TBranch        *b_Nu;   //!
  TBranch        *b_Xb;   //!
  TBranch        *b_vxec;   //!
  TBranch        *b_vyec;   //!
  TBranch        *b_vzec;   //!
  TBranch        *b_vxe;   //!
  TBranch        *b_vye;   //!
  TBranch        *b_vze;   //!
  TBranch        *b_Pex;   //!
  TBranch        *b_Pey;   //!
  TBranch        *b_Pez;   //!
  TBranch        *b_event;   //!
  TBranch        *b_Pe;   //!
  TBranch        *b_Ee;   //!
  TBranch        *b_e_Ein;   //!
  TBranch        *b_e_Eout;   //!
  TBranch        *b_e_Epcal;   //!
  TBranch        *b_e_npheltcc;   //!
  TBranch        *b_e_nphehtcc;   //!
  TBranch        *b_helic;   //!
  TBranch        *b_e_chi2pid;   //!
  TBranch        *b_e_pcal_lu;   //!
  TBranch        *b_e_pcal_lv;   //!
  TBranch        *b_e_pcal_lw;   //!
  TBranch        *b_e_ecin_lu;   //!
  TBranch        *b_e_ecin_lv;   //!
  TBranch        *b_e_ecin_lw;   //!
  TBranch        *b_e_ecout_lu;   //!
  TBranch        *b_e_ecout_lv;   //!
  TBranch        *b_e_ecout_lw;   //!
  TBranch        *b_e_pcal_hx;   //!
  TBranch        *b_e_pcal_hy;   //!
  TBranch        *b_e_pcal_hz;   //!
  TBranch        *b_e_ecin_hx;   //!
  TBranch        *b_e_ecin_hy;   //!
  TBranch        *b_e_ecin_hz;   //!
  TBranch        *b_e_ecout_hx;   //!
  TBranch        *b_e_ecout_hy;   //!
  TBranch        *b_e_ecout_hz;   //!
  TBranch        *b_e_trajx_sl0;   //!
  TBranch        *b_e_trajx_sl1;   //!
  TBranch        *b_e_trajx_sl2;   //!
  TBranch        *b_e_trajx_sl3;   //!
  TBranch        *b_e_trajx_sl4;   //!
  TBranch        *b_e_trajx_sl5;   //!
  TBranch        *b_e_trajy_sl0;   //!
  TBranch        *b_e_trajy_sl1;   //!
  TBranch        *b_e_trajy_sl2;   //!
  TBranch        *b_e_trajy_sl3;   //!
  TBranch        *b_e_trajy_sl4;   //!
  TBranch        *b_e_trajy_sl5;   //!
  TBranch        *b_e_trajz_sl0;   //!
  TBranch        *b_e_trajz_sl1;   //!
  TBranch        *b_e_trajz_sl2;   //!
  TBranch        *b_e_trajz_sl3;   //!
  TBranch        *b_e_trajz_sl4;   //!
  TBranch        *b_e_trajz_sl5;   //!
  TBranch        *b_e_pathtof;   //!
  TBranch        *b_e_timetof;   //!
  TBranch        *b_e_sector_tof;   //!
  TBranch        *b_e_Beta;   //!
  TBranch        *b_STTime;   //!
  TBranch        *b_RFTime;   //!
  TBranch        *b_e_dcx_rot_0;   //!
  TBranch        *b_e_dcy_rot_0;   //!
  TBranch        *b_e_dcx_rot_1;   //!
  TBranch        *b_e_dcy_rot_1;   //!
  TBranch        *b_e_dcx_rot_2;   //!
  TBranch        *b_e_dcy_rot_2;   //!
  TBranch        *b_e_sector_ltcc;   //!
  TBranch        *b_e_sector_htcc;   //!
  TBranch        *b_e_sector_ecal;   //!
  TBranch        *b_revent;   //!
  TBranch        *b_e_dc_chi2;   //!
  TBranch        *b_e_ftof1ax;   //!
  TBranch        *b_e_ftof1ay;   //!
  TBranch        *b_e_ftof1az;   //!
  TBranch        *b_e_pcalx;   //!
  TBranch        *b_e_pcaly;   //!
  TBranch        *b_e_pcalz;   //!
  TBranch        *b_e_ecalx;   //!
  TBranch        *b_e_ecaly;   //!
  TBranch        *b_e_ecalz;   //!
  TBranch        *b_e_ltccx;   //!
  TBranch        *b_e_ltccy;   //!
  TBranch        *b_e_ltccz;   //!
  TBranch        *b_e_htccx;   //!
  TBranch        *b_e_htccy;   //!
  TBranch        *b_e_htccz;   //!
  TBranch        *b_e_ftof1bx;   //!
  TBranch        *b_e_ftof1by;   //!
  TBranch        *b_e_ftof1bz;   //!
  TBranch        *b_e_ftof2x;   //!
  TBranch        *b_e_ftof2y;   //!
  TBranch        *b_e_ftof2z;   //!
  TBranch        *b_helonline_hel;   //!
  TBranch        *b_helonline_helRaw;   //!
  TBranch        *b_helflip_hel;   //!
  TBranch        *b_helflip_helRaw;   //!
  TBranch        *b_helflip_event;   //!
  TBranch        *b_e_dc_status;   //!
  TBranch        *b_e_dc_ndf;   //!
  TBranch        *b_e_sector_dc;   //!
  TBranch        *b_e_statPart;   //!
  TBranch        *b_e_DCPx;   //!
  TBranch        *b_e_DCPy;   //!
  TBranch        *b_e_DCPz;   //!
  TBranch        *b_y;   //!
  TBranch        *b_th_e;   //!
  TBranch        *b_phi_e;   //!
  TBranch        *b_helicRaw;   //!
  TBranch        *b_fA;   //!
  TBranch        *b_fB;   //!
  TBranch        *b_fC;   //!
  TBranch        *b_fV;   //!
  TBranch        *b_fW;   //!
  TBranch        *b_mc_Q2;   //!
  TBranch        *b_mc_W;   //!
  TBranch        *b_mc_Nu;   //!
  TBranch        *b_mc_Xb;   //!
  TBranch        *b_mc_vxe;   //!
  TBranch        *b_mc_vye;   //!
  TBranch        *b_mc_vze;   //!
  TBranch        *b_mc_Pex;   //!
  TBranch        *b_mc_Pey;   //!
  TBranch        *b_mc_Pez;   //!
  TBranch        *b_mc_event;   //!
  TBranch        *b_e_mcmass;   //!
  TBranch        *b_mc_Pe;   //!
  TBranch        *b_mc_Ee;   //!
  TBranch        *b_mc_revent;   //!
  TBranch        *b_mc_y;   //!
  TBranch        *b_mc_th_e;   //!
  TBranch        *b_mc_phi_e;   //!
  TBranch        *b_mc_e_Beta;   //!
  TBranch        *b_mc_helic;   //!
  TBranch        *b_mc_fA;   //!
  TBranch        *b_mc_fB;   //!
  TBranch        *b_mc_fC;   //!
  TBranch        *b_mc_fV;   //!
  TBranch        *b_mc_fW;   //!

  BSA_survey_cls(TString infile = "", TString binfo = "binning_info.txt");
  virtual ~BSA_survey_cls();
  virtual Int_t   Cut(Long64_t entry);
  virtual Bool_t  DIS();
  virtual Bool_t  eFID_ec();
  virtual Bool_t  eFID_dc();
  virtual Bool_t  piFID_ec(int k = 0);
  virtual Bool_t  pipFID_dc(int k = 0);
  virtual Bool_t  pimFID_dc(int k = 0);
  virtual Bool_t  FWD(int k = 0);
  virtual Bool_t  CF(int k = 0);
  virtual Int_t   GetEntry(Long64_t entry);
  virtual Float_t getALU(TString hpname, TString hnname, TString pv, TString tv, Float_t &val, Float_t &err);
  virtual Float_t getALU2D(TString bn);
    
  virtual Int_t   fillHist(TString hname, Float_t value = -111111);
  virtual Int_t   fillHist2D(TString hname, Float_t x, Float_t y);
  virtual Int_t   configHisto(TH1D *h, TString xtitle, TString ytitle,Color_t c = kBlack, EMarkerStyle ms = kFullDotLarge);
  virtual Int_t   setStyle();
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void    Init(TChain *tree, TString binfo);
  virtual void    Loop();
  virtual Bool_t  Notify();
  virtual void    Show(Long64_t entry = -1);
};

#endif

#ifdef BSA_survey_cls_cxx
BSA_survey_cls::BSA_survey_cls(TString infile, TString binfo) : fChain(0) 
{
  TChain *tch = new TChain();
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  TFile *f;
  if (infile == "") {
    tch->Add("../filtering/outfiles/pippim_sim_ga_1.root/outdata");
  }

  else {
    tch->Add(infile + "/outdata");
  }
  Init(tch,binfo);
}

BSA_survey_cls::~BSA_survey_cls()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t BSA_survey_cls::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t BSA_survey_cls::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
    for (auto& x: brv){
      x.second->UpdateFormulaLeaves();
    }
  }
  return centry;
}

void BSA_survey_cls::Init(TChain *tree,TString binfo)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  bm = new TBenchmark();
  bm->Start("main");
  // Set branch addresses and branch pointers
  if (!tree) return;
  setStyle();
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  
  fChain->SetBranchAddress("npart", &npart, &b_npart);
  fChain->SetBranchAddress("M", M, &b_M);
  fChain->SetBranchAddress("Phx", Phx, &b_Phx);
  fChain->SetBranchAddress("Phy", Phy, &b_Phy);
  fChain->SetBranchAddress("Phz", Phz, &b_Phz);
  fChain->SetBranchAddress("Z", Z, &b_Z);
  fChain->SetBranchAddress("Cospq", Cospq, &b_Cospq);
  fChain->SetBranchAddress("Pt2", Pt2, &b_Pt2);
  fChain->SetBranchAddress("Event", Event, &b_Event);
  fChain->SetBranchAddress("M2_01", M2_01, &b_M2_01);
  fChain->SetBranchAddress("M2_02", M2_02, &b_M2_02);
  fChain->SetBranchAddress("phiH", phiH, &b_phiH);
  fChain->SetBranchAddress("phiR", phiR, &b_phiR);
  fChain->SetBranchAddress("helic002_phiR", helic002_phiR, &b_helic002_phiR);
  fChain->SetBranchAddress("helic005_phiR", helic005_phiR, &b_helic005_phiR);
  fChain->SetBranchAddress("helic010_phiR", helic010_phiR, &b_helic010_phiR);
  fChain->SetBranchAddress("helic020_phiR", helic020_phiR, &b_helic020_phiR);
  fChain->SetBranchAddress("helic002_phiRst", helic002_phiRst, &b_helic002_phiRst);
  fChain->SetBranchAddress("helic005_phiRst", helic005_phiRst, &b_helic005_phiRst);
  fChain->SetBranchAddress("helic010_phiRst", helic010_phiRst, &b_helic010_phiRst);
  fChain->SetBranchAddress("helic020_phiRst", helic020_phiRst, &b_helic020_phiRst);
  fChain->SetBranchAddress("wUxS_phiR", wUxS_phiR, &b_wUxS_phiR);
  fChain->SetBranchAddress("R", R, &b_R);
  fChain->SetBranchAddress("KF", KF, &b_KF);
  fChain->SetBranchAddress("Mx2", Mx2, &b_Mx2);
  fChain->SetBranchAddress("xF", xF, &b_xF);
  fChain->SetBranchAddress("xF0", xF0, &b_xF0);
  fChain->SetBranchAddress("xF1", xF1, &b_xF1);
  fChain->SetBranchAddress("plcm", plcm, &b_plcm);
  fChain->SetBranchAddress("plcm0", plcm0, &b_plcm0);
  fChain->SetBranchAddress("plcm1", plcm1, &b_plcm1);
  fChain->SetBranchAddress("Eh", Eh, &b_Eh);
  fChain->SetBranchAddress("xFm", xFm, &b_xFm);
  fChain->SetBranchAddress("xFm0", xFm0, &b_xFm0);
  fChain->SetBranchAddress("xFm1", xFm1, &b_xFm1);
  fChain->SetBranchAddress("theta0", theta0, &b_theta0);
  fChain->SetBranchAddress("theta1", theta1, &b_theta1);
  fChain->SetBranchAddress("cos_theta_P0cm", cos_theta_P0cm, &b_cos_theta_P0cm);
  fChain->SetBranchAddress("sin_theta_P0cm", sin_theta_P0cm, &b_sin_theta_P0cm);
  fChain->SetBranchAddress("xFo", xFo, &b_xFo);
  fChain->SetBranchAddress("xFo0", xFo0, &b_xFo0);
  fChain->SetBranchAddress("xFo1", xFo1, &b_xFo1);
  fChain->SetBranchAddress("phiH_phiR", phiH_phiR, &b_phiH_phiR);
  fChain->SetBranchAddress("phiR_cov", phiR_cov, &b_phiR_cov);
  fChain->SetBranchAddress("p0T2", p0T2, &b_p0T2);
  fChain->SetBranchAddress("p1T2", p1T2, &b_p1T2);
  fChain->SetBranchAddress("phipq", phipq, &b_phipq);
  fChain->SetBranchAddress("phT2", phT2, &b_phT2);
  fChain->SetBranchAddress("etaCM0", etaCM0, &b_etaCM0);
  fChain->SetBranchAddress("etaCM1", etaCM1, &b_etaCM1);
  fChain->SetBranchAddress("etaBF0p", etaBF0p, &b_etaBF0p);
  fChain->SetBranchAddress("etaBF1p", etaBF1p, &b_etaBF1p);
  fChain->SetBranchAddress("etaBF0m", etaBF0m, &b_etaBF0m);
  fChain->SetBranchAddress("etaBF1m", etaBF1m, &b_etaBF1m);
  fChain->SetBranchAddress("etaBF0", etaBF0, &b_etaBF0);
  fChain->SetBranchAddress("etaBF1", etaBF1, &b_etaBF1);
  fChain->SetBranchAddress("phiR_ha", phiR_ha, &b_phiR_ha);
  fChain->SetBranchAddress("plcm0_r", plcm0_r, &b_plcm0_r);
  fChain->SetBranchAddress("plcm1_r", plcm1_r, &b_plcm1_r);
  fChain->SetBranchAddress("phiR_covH", phiR_covH, &b_phiR_covH);
  fChain->SetBranchAddress("E0_phcm", E0_phcm, &b_E0_phcm);
  fChain->SetBranchAddress("E1_phcm", E1_phcm, &b_E1_phcm);
  fChain->SetBranchAddress("mix_npart", &mix_npart, &b_mix_npart);
  fChain->SetBranchAddress("mc_mix_npart", &mc_mix_npart, &b_mc_mix_npart);
  fChain->SetBranchAddress("det", &det_, &b_det_);
  fChain->SetBranchAddress("det.fUniqueID", det_fUniqueID, &b_det_fUniqueID);
  fChain->SetBranchAddress("det.fBits", det_fBits, &b_det_fBits);
  fChain->SetBranchAddress("det.beta[6]", det_beta, &b_det_beta);
  fChain->SetBranchAddress("det.m2b[6]", det_m2b, &b_det_m2b);
  fChain->SetBranchAddress("det.vx[6]", det_vx, &b_det_vx);
  fChain->SetBranchAddress("det.vy[6]", det_vy, &b_det_vy);
  fChain->SetBranchAddress("det.vz[6]", det_vz, &b_det_vz);
  fChain->SetBranchAddress("det.dcx[6]", det_dcx, &b_det_dcx);
  fChain->SetBranchAddress("det.dcy[6]", det_dcy, &b_det_dcy);
  fChain->SetBranchAddress("det.dcz[6]", det_dcz, &b_det_dcz);
  fChain->SetBranchAddress("det.statPart[6]", det_statPart, &b_det_statPart);
  fChain->SetBranchAddress("det.dc_chi2[6]", det_dc_chi2, &b_det_dc_chi2);
  fChain->SetBranchAddress("det.dc_ndf[6]", det_dc_ndf, &b_det_dc_ndf);
  fChain->SetBranchAddress("det.pcal_lu[6]", det_pcal_lu, &b_det_pcal_lu);
  fChain->SetBranchAddress("det.pcal_lv[6]", det_pcal_lv, &b_det_pcal_lv);
  fChain->SetBranchAddress("det.pcal_lw[6]", det_pcal_lw, &b_det_pcal_lw);
  fChain->SetBranchAddress("pdata", &pdata_, &b_pdata_);
  fChain->SetBranchAddress("pdata.fUniqueID", pdata_fUniqueID, &b_pdata_fUniqueID);
  fChain->SetBranchAddress("pdata.fBits", pdata_fBits, &b_pdata_fBits);
  fChain->SetBranchAddress("pdata.e[6]", pdata_e, &b_pdata_e);
  fChain->SetBranchAddress("pdata.px[6]", pdata_px, &b_pdata_px);
  fChain->SetBranchAddress("pdata.py[6]", pdata_py, &b_pdata_py);
  fChain->SetBranchAddress("pdata.pz[6]", pdata_pz, &b_pdata_pz);
  fChain->SetBranchAddress("pdata.helic002_phiH[6]", pdata_helic002_phiH, &b_pdata_helic002_phiH);
  fChain->SetBranchAddress("pdata.helic005_phiH[6]", pdata_helic005_phiH, &b_pdata_helic005_phiH);
  fChain->SetBranchAddress("pdata.helic010_phiH[6]", pdata_helic010_phiH, &b_pdata_helic010_phiH);
  fChain->SetBranchAddress("pdata.helic020_phiH[6]", pdata_helic020_phiH, &b_pdata_helic020_phiH);
  fChain->SetBranchAddress("pdata.phiHs[6]", pdata_phiHs, &b_pdata_phiHs);
  fChain->SetBranchAddress("mc_pdata", &mc_pdata_, &b_mc_pdata_);
  fChain->SetBranchAddress("mc_pdata.fUniqueID", mc_pdata_fUniqueID, &b_mc_pdata_fUniqueID);
  fChain->SetBranchAddress("mc_pdata.fBits", mc_pdata_fBits, &b_mc_pdata_fBits);
  fChain->SetBranchAddress("mc_pdata.e[6]", mc_pdata_e, &b_mc_pdata_e);
  fChain->SetBranchAddress("mc_pdata.px[6]", mc_pdata_px, &b_mc_pdata_px);
  fChain->SetBranchAddress("mc_pdata.py[6]", mc_pdata_py, &b_mc_pdata_py);
  fChain->SetBranchAddress("mc_pdata.pz[6]", mc_pdata_pz, &b_mc_pdata_pz);
  fChain->SetBranchAddress("mc_pdata.helic002_phiH[6]", mc_pdata_helic002_phiH, &b_mc_pdata_helic002_phiH);
  fChain->SetBranchAddress("mc_pdata.helic005_phiH[6]", mc_pdata_helic005_phiH, &b_mc_pdata_helic005_phiH);
  fChain->SetBranchAddress("mc_pdata.helic010_phiH[6]", mc_pdata_helic010_phiH, &b_mc_pdata_helic010_phiH);
  fChain->SetBranchAddress("mc_pdata.helic020_phiH[6]", mc_pdata_helic020_phiH, &b_mc_pdata_helic020_phiH);
  fChain->SetBranchAddress("mc_pdata.phiHs[6]", mc_pdata_phiHs, &b_mc_pdata_phiHs);
  fChain->SetBranchAddress("mc_npart", &mc_npart, &b_mc_npart);
  fChain->SetBranchAddress("mc_M", mc_M, &b_mc_M);
  fChain->SetBranchAddress("mc_Phx", mc_Phx, &b_mc_Phx);
  fChain->SetBranchAddress("mc_Phy", mc_Phy, &b_mc_Phy);
  fChain->SetBranchAddress("mc_Phz", mc_Phz, &b_mc_Phz);
  fChain->SetBranchAddress("mc_Z", mc_Z, &b_mc_Z);
  fChain->SetBranchAddress("mc_Cospq", mc_Cospq, &b_mc_Cospq);
  fChain->SetBranchAddress("mc_Pt2", mc_Pt2, &b_mc_Pt2);
  fChain->SetBranchAddress("mc_Event", mc_Event, &b_mc_Event);
  fChain->SetBranchAddress("mc_M2_01", mc_M2_01, &b_mc_M2_01);
  fChain->SetBranchAddress("mc_M2_02", mc_M2_02, &b_mc_M2_02);
  fChain->SetBranchAddress("mc_phiH", mc_phiH, &b_mc_phiH);
  fChain->SetBranchAddress("mc_phiR", mc_phiR, &b_mc_phiR);
  fChain->SetBranchAddress("mc_helic002_phiR", mc_helic002_phiR, &b_mc_helic002_phiR);
  fChain->SetBranchAddress("mc_helic005_phiR", mc_helic005_phiR, &b_mc_helic005_phiR);
  fChain->SetBranchAddress("mc_helic010_phiR", mc_helic010_phiR, &b_mc_helic010_phiR);
  fChain->SetBranchAddress("mc_helic020_phiR", mc_helic020_phiR, &b_mc_helic020_phiR);
  fChain->SetBranchAddress("mc_helic002_phiRst", mc_helic002_phiRst, &b_mc_helic002_phiRst);
  fChain->SetBranchAddress("mc_helic005_phiRst", mc_helic005_phiRst, &b_mc_helic005_phiRst);
  fChain->SetBranchAddress("mc_helic010_phiRst", mc_helic010_phiRst, &b_mc_helic010_phiRst);
  fChain->SetBranchAddress("mc_helic020_phiRst", mc_helic020_phiRst, &b_mc_helic020_phiRst);
  fChain->SetBranchAddress("mc_wUxS_phiR", mc_wUxS_phiR, &b_mc_wUxS_phiR);
  fChain->SetBranchAddress("mc_Mx2", mc_Mx2, &b_mc_Mx2);
  fChain->SetBranchAddress("mc_xF", mc_xF, &b_mc_xF);
  fChain->SetBranchAddress("mc_xF0", mc_xF0, &b_mc_xF0);
  fChain->SetBranchAddress("mc_xF1", mc_xF1, &b_mc_xF1);
  fChain->SetBranchAddress("mc_plcm", mc_plcm, &b_mc_plcm);
  fChain->SetBranchAddress("mc_plcm0", mc_plcm0, &b_mc_plcm0);
  fChain->SetBranchAddress("mc_plcm1", mc_plcm1, &b_mc_plcm1);
  fChain->SetBranchAddress("mc_Eh", mc_Eh, &b_mc_Eh);
  fChain->SetBranchAddress("mc_xFm", mc_xFm, &b_mc_xFm);
  fChain->SetBranchAddress("mc_xFm0", mc_xFm0, &b_mc_xFm0);
  fChain->SetBranchAddress("mc_xFm1", mc_xFm1, &b_mc_xFm1);
  fChain->SetBranchAddress("mc_theta0", mc_theta0, &b_mc_theta0);
  fChain->SetBranchAddress("mc_theta1", mc_theta1, &b_mc_theta1);
  fChain->SetBranchAddress("mc_cos_theta_P0cm", mc_cos_theta_P0cm, &b_mc_cos_theta_P0cm);
  fChain->SetBranchAddress("mc_sin_theta_P0cm", mc_sin_theta_P0cm, &b_mc_sin_theta_P0cm);
  fChain->SetBranchAddress("mc_xFo", mc_xFo, &b_mc_xFo);
  fChain->SetBranchAddress("mc_xFo0", mc_xFo0, &b_mc_xFo0);
  fChain->SetBranchAddress("mc_xFo1", mc_xFo1, &b_mc_xFo1);
  fChain->SetBranchAddress("mc_phiH_phiR", mc_phiH_phiR, &b_mc_phiH_phiR);
  fChain->SetBranchAddress("mc_phiR_cov", mc_phiR_cov, &b_mc_phiR_cov);
  fChain->SetBranchAddress("mc_p0T2", mc_p0T2, &b_mc_p0T2);
  fChain->SetBranchAddress("mc_p1T2", mc_p1T2, &b_mc_p1T2);
  fChain->SetBranchAddress("mc_phipq", mc_phipq, &b_mc_phipq);
  fChain->SetBranchAddress("mc_phT2", mc_phT2, &b_mc_phT2);
  fChain->SetBranchAddress("mc_etaCM0", mc_etaCM0, &b_mc_etaCM0);
  fChain->SetBranchAddress("mc_etaCM1", mc_etaCM1, &b_mc_etaCM1);
  fChain->SetBranchAddress("mc_etaBF0p", mc_etaBF0p, &b_mc_etaBF0p);
  fChain->SetBranchAddress("mc_etaBF1p", mc_etaBF1p, &b_mc_etaBF1p);
  fChain->SetBranchAddress("mc_etaBF0m", mc_etaBF0m, &b_mc_etaBF0m);
  fChain->SetBranchAddress("mc_etaBF1m", mc_etaBF1m, &b_mc_etaBF1m);
  fChain->SetBranchAddress("mc_etaBF0", mc_etaBF0, &b_mc_etaBF0);
  fChain->SetBranchAddress("mc_etaBF1", mc_etaBF1, &b_mc_etaBF1);
  fChain->SetBranchAddress("mc_phiR_ha", mc_phiR_ha, &b_mc_phiR_ha);
  fChain->SetBranchAddress("mc_plcm0_r", mc_plcm0_r, &b_mc_plcm0_r);
  fChain->SetBranchAddress("mc_plcm1_r", mc_plcm1_r, &b_mc_plcm1_r);
  fChain->SetBranchAddress("mc_phiR_covH", mc_phiR_covH, &b_mc_phiR_covH);
  fChain->SetBranchAddress("mc_E0_phcm", mc_E0_phcm, &b_mc_E0_phcm);
  fChain->SetBranchAddress("mc_E1_phcm", mc_E1_phcm, &b_mc_E1_phcm);
  fChain->SetBranchAddress("mc_R", mc_R, &b_mc_R);
  fChain->SetBranchAddress("mc_KF", mc_KF, &b_mc_KF);
  fChain->SetBranchAddress("Q2", &Q2, &b_Q2);
  fChain->SetBranchAddress("W", &W, &b_W);
  fChain->SetBranchAddress("Nu", &Nu, &b_Nu);
  fChain->SetBranchAddress("Xb", &Xb, &b_Xb);
  fChain->SetBranchAddress("vxec", &vxec, &b_vxec);
  fChain->SetBranchAddress("vyec", &vyec, &b_vyec);
  fChain->SetBranchAddress("vzec", &vzec, &b_vzec);
  fChain->SetBranchAddress("vxe", &vxe, &b_vxe);
  fChain->SetBranchAddress("vye", &vye, &b_vye);
  fChain->SetBranchAddress("vze", &vze, &b_vze);
  fChain->SetBranchAddress("Pex", &Pex, &b_Pex);
  fChain->SetBranchAddress("Pey", &Pey, &b_Pey);
  fChain->SetBranchAddress("Pez", &Pez, &b_Pez);
  fChain->SetBranchAddress("event", &event, &b_event);
  fChain->SetBranchAddress("Pe", &Pe, &b_Pe);
  fChain->SetBranchAddress("Ee", &Ee, &b_Ee);
  fChain->SetBranchAddress("e_Ein", &e_Ein, &b_e_Ein);
  fChain->SetBranchAddress("e_Eout", &e_Eout, &b_e_Eout);
  fChain->SetBranchAddress("e_Epcal", &e_Epcal, &b_e_Epcal);
  fChain->SetBranchAddress("e_npheltcc", &e_npheltcc, &b_e_npheltcc);
  fChain->SetBranchAddress("e_nphehtcc", &e_nphehtcc, &b_e_nphehtcc);
  fChain->SetBranchAddress("helic", &helic, &b_helic);
  fChain->SetBranchAddress("e_chi2pid", &e_chi2pid, &b_e_chi2pid);
  fChain->SetBranchAddress("e_pcal_lu", &e_pcal_lu, &b_e_pcal_lu);
  fChain->SetBranchAddress("e_pcal_lv", &e_pcal_lv, &b_e_pcal_lv);
  fChain->SetBranchAddress("e_pcal_lw", &e_pcal_lw, &b_e_pcal_lw);
  fChain->SetBranchAddress("e_ecin_lu", &e_ecin_lu, &b_e_ecin_lu);
  fChain->SetBranchAddress("e_ecin_lv", &e_ecin_lv, &b_e_ecin_lv);
  fChain->SetBranchAddress("e_ecin_lw", &e_ecin_lw, &b_e_ecin_lw);
  fChain->SetBranchAddress("e_ecout_lu", &e_ecout_lu, &b_e_ecout_lu);
  fChain->SetBranchAddress("e_ecout_lv", &e_ecout_lv, &b_e_ecout_lv);
  fChain->SetBranchAddress("e_ecout_lw", &e_ecout_lw, &b_e_ecout_lw);
  fChain->SetBranchAddress("e_pcal_hx", &e_pcal_hx, &b_e_pcal_hx);
  fChain->SetBranchAddress("e_pcal_hy", &e_pcal_hy, &b_e_pcal_hy);
  fChain->SetBranchAddress("e_pcal_hz", &e_pcal_hz, &b_e_pcal_hz);
  fChain->SetBranchAddress("e_ecin_hx", &e_ecin_hx, &b_e_ecin_hx);
  fChain->SetBranchAddress("e_ecin_hy", &e_ecin_hy, &b_e_ecin_hy);
  fChain->SetBranchAddress("e_ecin_hz", &e_ecin_hz, &b_e_ecin_hz);
  fChain->SetBranchAddress("e_ecout_hx", &e_ecout_hx, &b_e_ecout_hx);
  fChain->SetBranchAddress("e_ecout_hy", &e_ecout_hy, &b_e_ecout_hy);
  fChain->SetBranchAddress("e_ecout_hz", &e_ecout_hz, &b_e_ecout_hz);
  fChain->SetBranchAddress("e_trajx_sl0", &e_trajx_sl0, &b_e_trajx_sl0);
  fChain->SetBranchAddress("e_trajx_sl1", &e_trajx_sl1, &b_e_trajx_sl1);
  fChain->SetBranchAddress("e_trajx_sl2", &e_trajx_sl2, &b_e_trajx_sl2);
  fChain->SetBranchAddress("e_trajx_sl3", &e_trajx_sl3, &b_e_trajx_sl3);
  fChain->SetBranchAddress("e_trajx_sl4", &e_trajx_sl4, &b_e_trajx_sl4);
  fChain->SetBranchAddress("e_trajx_sl5", &e_trajx_sl5, &b_e_trajx_sl5);
  fChain->SetBranchAddress("e_trajy_sl0", &e_trajy_sl0, &b_e_trajy_sl0);
  fChain->SetBranchAddress("e_trajy_sl1", &e_trajy_sl1, &b_e_trajy_sl1);
  fChain->SetBranchAddress("e_trajy_sl2", &e_trajy_sl2, &b_e_trajy_sl2);
  fChain->SetBranchAddress("e_trajy_sl3", &e_trajy_sl3, &b_e_trajy_sl3);
  fChain->SetBranchAddress("e_trajy_sl4", &e_trajy_sl4, &b_e_trajy_sl4);
  fChain->SetBranchAddress("e_trajy_sl5", &e_trajy_sl5, &b_e_trajy_sl5);
  fChain->SetBranchAddress("e_trajz_sl0", &e_trajz_sl0, &b_e_trajz_sl0);
  fChain->SetBranchAddress("e_trajz_sl1", &e_trajz_sl1, &b_e_trajz_sl1);
  fChain->SetBranchAddress("e_trajz_sl2", &e_trajz_sl2, &b_e_trajz_sl2);
  fChain->SetBranchAddress("e_trajz_sl3", &e_trajz_sl3, &b_e_trajz_sl3);
  fChain->SetBranchAddress("e_trajz_sl4", &e_trajz_sl4, &b_e_trajz_sl4);
  fChain->SetBranchAddress("e_trajz_sl5", &e_trajz_sl5, &b_e_trajz_sl5);
  fChain->SetBranchAddress("e_pathtof", &e_pathtof, &b_e_pathtof);
  fChain->SetBranchAddress("e_timetof", &e_timetof, &b_e_timetof);
  fChain->SetBranchAddress("e_sector_tof", &e_sector_tof, &b_e_sector_tof);
  fChain->SetBranchAddress("e_Beta", &e_Beta, &b_e_Beta);
  fChain->SetBranchAddress("STTime", &STTime, &b_STTime);
  fChain->SetBranchAddress("RFTime", &RFTime, &b_RFTime);
  fChain->SetBranchAddress("e_dcx_rot_0", &e_dcx_rot_0, &b_e_dcx_rot_0);
  fChain->SetBranchAddress("e_dcy_rot_0", &e_dcy_rot_0, &b_e_dcy_rot_0);
  fChain->SetBranchAddress("e_dcx_rot_1", &e_dcx_rot_1, &b_e_dcx_rot_1);
  fChain->SetBranchAddress("e_dcy_rot_1", &e_dcy_rot_1, &b_e_dcy_rot_1);
  fChain->SetBranchAddress("e_dcx_rot_2", &e_dcx_rot_2, &b_e_dcx_rot_2);
  fChain->SetBranchAddress("e_dcy_rot_2", &e_dcy_rot_2, &b_e_dcy_rot_2);
  fChain->SetBranchAddress("e_sector_ltcc", &e_sector_ltcc, &b_e_sector_ltcc);
  fChain->SetBranchAddress("e_sector_htcc", &e_sector_htcc, &b_e_sector_htcc);
  fChain->SetBranchAddress("e_sector_ecal", &e_sector_ecal, &b_e_sector_ecal);
  fChain->SetBranchAddress("revent", &revent, &b_revent);
  fChain->SetBranchAddress("e_dc_chi2", &e_dc_chi2, &b_e_dc_chi2);
  fChain->SetBranchAddress("e_ftof1ax", &e_ftof1ax, &b_e_ftof1ax);
  fChain->SetBranchAddress("e_ftof1ay", &e_ftof1ay, &b_e_ftof1ay);
  fChain->SetBranchAddress("e_ftof1az", &e_ftof1az, &b_e_ftof1az);
  fChain->SetBranchAddress("e_pcalx", &e_pcalx, &b_e_pcalx);
  fChain->SetBranchAddress("e_pcaly", &e_pcaly, &b_e_pcaly);
  fChain->SetBranchAddress("e_pcalz", &e_pcalz, &b_e_pcalz);
  fChain->SetBranchAddress("e_ecalx", &e_ecalx, &b_e_ecalx);
  fChain->SetBranchAddress("e_ecaly", &e_ecaly, &b_e_ecaly);
  fChain->SetBranchAddress("e_ecalz", &e_ecalz, &b_e_ecalz);
  fChain->SetBranchAddress("e_ltccx", &e_ltccx, &b_e_ltccx);
  fChain->SetBranchAddress("e_ltccy", &e_ltccy, &b_e_ltccy);
  fChain->SetBranchAddress("e_ltccz", &e_ltccz, &b_e_ltccz);
  fChain->SetBranchAddress("e_htccx", &e_htccx, &b_e_htccx);
  fChain->SetBranchAddress("e_htccy", &e_htccy, &b_e_htccy);
  fChain->SetBranchAddress("e_htccz", &e_htccz, &b_e_htccz);
  fChain->SetBranchAddress("e_ftof1bx", &e_ftof1bx, &b_e_ftof1bx);
  fChain->SetBranchAddress("e_ftof1by", &e_ftof1by, &b_e_ftof1by);
  fChain->SetBranchAddress("e_ftof1bz", &e_ftof1bz, &b_e_ftof1bz);
  fChain->SetBranchAddress("e_ftof2x", &e_ftof2x, &b_e_ftof2x);
  fChain->SetBranchAddress("e_ftof2y", &e_ftof2y, &b_e_ftof2y);
  fChain->SetBranchAddress("e_ftof2z", &e_ftof2z, &b_e_ftof2z);
  fChain->SetBranchAddress("helonline_hel", &helonline_hel, &b_helonline_hel);
  fChain->SetBranchAddress("helonline_helRaw", &helonline_helRaw, &b_helonline_helRaw);
  fChain->SetBranchAddress("helflip_hel", &helflip_hel, &b_helflip_hel);
  fChain->SetBranchAddress("helflip_helRaw", &helflip_helRaw, &b_helflip_helRaw);
  fChain->SetBranchAddress("helflip_event", &helflip_event, &b_helflip_event);
  fChain->SetBranchAddress("e_dc_status", &e_dc_status, &b_e_dc_status);
  fChain->SetBranchAddress("e_dc_ndf", &e_dc_ndf, &b_e_dc_ndf);
  fChain->SetBranchAddress("e_sector_dc", &e_sector_dc, &b_e_sector_dc);
  fChain->SetBranchAddress("e_statPart", &e_statPart, &b_e_statPart);
  fChain->SetBranchAddress("e_DCPx", &e_DCPx, &b_e_DCPx);
  fChain->SetBranchAddress("e_DCPy", &e_DCPy, &b_e_DCPy);
  fChain->SetBranchAddress("e_DCPz", &e_DCPz, &b_e_DCPz);
  fChain->SetBranchAddress("y", &y, &b_y);
  fChain->SetBranchAddress("th_e", &th_e, &b_th_e);
  fChain->SetBranchAddress("phi_e", &phi_e, &b_phi_e);
  fChain->SetBranchAddress("helicRaw", &helicRaw, &b_helicRaw);
  fChain->SetBranchAddress("fA", &fA, &b_fA);
  fChain->SetBranchAddress("fB", &fB, &b_fB);
  fChain->SetBranchAddress("fC", &fC, &b_fC);
  fChain->SetBranchAddress("fV", &fV, &b_fV);
  fChain->SetBranchAddress("fW", &fW, &b_fW);
  fChain->SetBranchAddress("mc_Q2", &mc_Q2, &b_mc_Q2);
  fChain->SetBranchAddress("mc_W", &mc_W, &b_mc_W);
  fChain->SetBranchAddress("mc_Nu", &mc_Nu, &b_mc_Nu);
  fChain->SetBranchAddress("mc_Xb", &mc_Xb, &b_mc_Xb);
  fChain->SetBranchAddress("mc_vxe", &mc_vxe, &b_mc_vxe);
  fChain->SetBranchAddress("mc_vye", &mc_vye, &b_mc_vye);
  fChain->SetBranchAddress("mc_vze", &mc_vze, &b_mc_vze);
  fChain->SetBranchAddress("mc_Pex", &mc_Pex, &b_mc_Pex);
  fChain->SetBranchAddress("mc_Pey", &mc_Pey, &b_mc_Pey);
  fChain->SetBranchAddress("mc_Pez", &mc_Pez, &b_mc_Pez);
  fChain->SetBranchAddress("mc_event", &mc_event, &b_mc_event);
  fChain->SetBranchAddress("e_mcmass", &e_mcmass, &b_e_mcmass);
  fChain->SetBranchAddress("mc_Pe", &mc_Pe, &b_mc_Pe);
  fChain->SetBranchAddress("mc_Ee", &mc_Ee, &b_mc_Ee);
  fChain->SetBranchAddress("mc_revent", &mc_revent, &b_mc_revent);
  fChain->SetBranchAddress("mc_y", &mc_y, &b_mc_y);
  fChain->SetBranchAddress("mc_th_e", &mc_th_e, &b_mc_th_e);
  fChain->SetBranchAddress("mc_phi_e", &mc_phi_e, &b_mc_phi_e);
  fChain->SetBranchAddress("mc_e_Beta", &mc_e_Beta, &b_mc_e_Beta);
  fChain->SetBranchAddress("mc_helic", &mc_helic, &b_mc_helic);
  fChain->SetBranchAddress("mc_fA", &mc_fA, &b_mc_fA);
  fChain->SetBranchAddress("mc_fB", &mc_fB, &b_mc_fB);
  fChain->SetBranchAddress("mc_fC", &mc_fC, &b_mc_fC);
  fChain->SetBranchAddress("mc_fV", &mc_fV, &b_mc_fV);
  fChain->SetBranchAddress("mc_fW", &mc_fW, &b_mc_fW);

  OUTDIR = "";
  
  ///// Setting bins and output dir ////
  
  std::ifstream bf(binfo);
  TString line;
  char cln[500];
  bf.getline(cln,500);
  bf>>pltv>>ttlv;
  fChain->GetEntries();
  std::cout<<"#### Entries to be processed "<<  fChain->GetEntries()<<" ###"<<std::endl;
  std::cout<<"#### binning configuration ####"<<std::endl;
  std::cout<<"plot var: " + pltv + ", title_var: " + ttlv<<std::endl;
  fChain->LoadTree(0);
  while (bf.getline(cln,500)){
    line = cln;
    line.ReplaceAll("\t"," ");
    line = line.Strip();
    if ( line[0]=='#' || line.Length() == 0) continue;
    TString tok="";
    Ssiz_t from=0;
    line.Tokenize(tok,from);
    TString bn = tok;
    line.Tokenize(tok,from);
    TString bfrm = tok;
    std::cout<<bn<<" : "<<bfrm<<"\t";
    brv[bn] = new TTreeFormula(bn,bfrm,fChain);
    fChain->SetNotify(brv[bn]);
    //    brv[bn] = (Float_t *)(fChain->GetBranch(bn))->GetAddress();

    OUTDIR += bn.ReplaceAll("(","_").ReplaceAll(")","_");
    while(line.Tokenize(tok,from)){
      bedg[bn].push_back(atof(tok));
      std::cout<<atof(tok)<<"\t";
    }
    std::cout<<"\n";
  }
  OUTDIR += "_bins";

  time_t timestmp;
  time(&timestmp);
  
  struct tm *lt = localtime(&timestmp);
  TString date = Form("%02d%02d%02d_%02d%02d%02d",lt->tm_year%100,lt->tm_mon+1,lt->tm_mday,lt->tm_hour,lt->tm_min,lt->tm_sec);

  OUTDIR += Form("_%02d%02d%02d",lt->tm_year%100,lt->tm_mon+1,lt->tm_mday);
  system("mkdir " + OUTDIR + " 2>/dev/null");

  ///// end Setting bins and output dir ////
  
  ofile = new TFile(OUTDIR + "/BSA_results.root","recreate");
  Int_t Nb_phi = 12;
  Float_t min_phi = 0, max_phi = 360;
  
  ofile->Add(new TH1D("hM","#pi^{+}#pi^{-} : IM",300,0,3));

  ofile->Add(new TH1D("hp_phiR","#pi^{+}#pi^{-} : #phi_{R#perp} (#lambda = +1)",Nb_phi, min_phi, max_phi));
  ofile->Add(new TH1D("hn_phiR","#pi^{+}#pi^{-} : #phi_{R#perp} (#lambda = -1)",Nb_phi, min_phi, max_phi));

  ofile->Add(new TH1D("hp_phiH","#pi^{+} : #phi_{H} (#lambda = +1)",Nb_phi, min_phi, max_phi));
  ofile->Add(new TH1D("hn_phiH","#pi^{+} : #phi_{H} (#lambda = -1)",Nb_phi, min_phi, max_phi));

  std::cout<<"#######\nALU histograms : "<<std::endl;
      
  for (auto& x: bedg){
    TString bn = x.first;
    for (int k=0;k<x.second.size()-1;k++){
      TString hname = pltv + "_" + bn + Form("_b%d",k);
      TString hnamepip = "phiH_" + bn + Form("_b%d",k);
      TString ttlsuf =  Form("%.2f<%s<%.2f",x.second[k], bn.Data(), x.second[k+1]);
      ofile->Add(new TH1D("hp_"+ hname, "#pi^{+}#pi^{-} : " + ttlv + " (#lambda = +1, " + ttlsuf + ")",Nb_phi, min_phi, max_phi));
      ofile->Add(new TH1D("hn_"+ hname, "#pi^{+}#pi^{-} : " + ttlv + " (#lambda = -1, " + ttlsuf + ")",Nb_phi, min_phi, max_phi));

      ofile->Add(new TH1D("hp_" + hnamepip,"#pi^{+} : #phi_{H} (#lambda = +1, " + ttlsuf + ")",Nb_phi, min_phi, max_phi));
      ofile->Add(new TH1D("hn_" + hnamepip,"#pi^{+} : #phi_{H} (#lambda = -1, " + ttlsuf + ")",Nb_phi, min_phi, max_phi));
    }

    /// 2D histos ///

    ofile->Add(new TH2D("hsin"+pltv + "_" + bn + "_p","sin(" + ttlv +") vs " + bn,x.second.size()-1,x.second.data(),300,-1,1));
    ofile->Add(new TH2D("hsin"+pltv + "_" + bn + "_n","sin(" + ttlv +") vs " + bn,x.second.size()-1,x.second.data(),300,-1,1));
    ofile->Add(new TH2D("hsinphiH_" + bn + "_p","sin(#phi_{H}) vs " + bn,x.second.size()-1,x.second.data(),300,-1,1));
    ofile->Add(new TH2D("hsinphiH_" + bn + "_n","sin(#phi_{H}) vs " + bn,x.second.size()-1,x.second.data(),300,-1,1));
    
    /// end 2D histos ////

    /// ALU histos
    std::cout<<bn<<" : "<<brv[bn]->GetExpFormula()<<std::endl;

    ofile->Add(new TH1D("hALU_"+pltv + "_" + bn,"ALU^{sin(" + ttlv+ ")}",x.second.size()-1,x.second.data()));
    ofile->Add(new TH1D("hALU_phiH_" + bn,"ALU^{sin(#phi_{H})}",x.second.size()-1,x.second.data()));

    /// end ALU histos
    
  }
  
  ///// end REC data ////
  
  Notify();
}

Bool_t BSA_survey_cls::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void BSA_survey_cls::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t BSA_survey_cls::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef BSA_survey_cls_cxx
