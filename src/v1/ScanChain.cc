// A cms3 looper to compare delta-beta and effective area pileup corrections

#include <iostream>
#include "ScanChain.h"

using namespace std;
using namespace tas;

void PileupCorrection::ScanChain (TTree * tree, const char* outname, int MaxEvents) {
  const int nAbsIsoBins = 100;
  const float AbsIsoMin = 0.001;
  const float AbsIsoMax = 2;
  const int nCorrBins = 100;
  const float CorrMin = 0.001;
  const float CorrMax = 2;

  const int nRelIsoBins = 100;
  const float RelIsoMin = 0.001;
  const float RelIsoMax = 1;
  const int nCorrRelBins = 100;
  const float CorrRelMin = 0.001;
  const float CorrRelMax = 1;

  const int nCorrsBins = 24;
  const float CorrsMin = 0;
  const float CorrsMax = 8;
  const int nRelCorrsBins = 24;
  const float RelCorrsMin = 0;
  const float RelCorrsMax = 3;

  const int nDiffBins = 20;
  const int DiffMin = -4;
  const int DiffMax = 4;
  const int nRelDiffBins = 20;
  const float RelDiffMin = -2;
  const float RelDiffMax = 2;

  //////////
  // Mini //
  //////////

  // Absolute Iso
  TH1F h_all_EAcorrs("h_all_EAcorrs","All Effective Area Corrections",nCorrsBins,CorrsMin,CorrsMax);
  TH1F h_all_dBcorrs("h_all_dBcorrs","All #Delta#beta Corrections",nCorrsBins,CorrsMin,CorrsMax);

  TH1F h_all_EAs("h_all_EAs","All Effective Areas",20,0,10);
  TH1F h_all_rhos("h_all_rhos","All Rhos",20,0,10);

  vector<float> temp;
  temp.reserve(1000);
  int nEntriesByNpv[13] = {0};
  vector<vector<float> > EA_corrs_byNpv(13,temp);
  vector<vector<float> > dB_corrs_byNpv(13,temp);
  TH1F h_median_EA("h_median_EA","Median Effective Area Correction by Number of Pileup Vertices",12,0,60);
  TH1F h_median_dB("h_median_dB","Median #Delta#beta Correction by Number of Pileup Vertices",12,0,60);

  // Delta-Beta
  // < 20 vertices
  TH2F h_CI_dB_low_bar("h_CI_dB_low_bar","#Delta#beta Correction vs Uncorrected Absolute Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_CI_dB_low_cap("h_CI_dB_low_cap","#Delta#beta Correction vs Uncorrected Absolute Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax, nCorrBins, CorrMin, CorrMax);
  // 20 < vertices < 40
  TH2F h_CI_dB_med_bar("h_CI_dB_med_bar","#Delta#beta Correction vs Uncorrected Absolute Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_CI_dB_med_cap("h_CI_dB_med_cap","#Delta#beta Correction vs Uncorrected Absolute Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax, nCorrBins, CorrMin, CorrMax);
  // > 40 vertices
  TH2F h_CI_dB_high_bar("h_CI_dB_high_bar","#Delta#beta Correction vs Uncorrected Absolute Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_CI_dB_high_cap("h_CI_dB_high_cap","#Delta#beta Correction vs Uncorrected Absolute Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax, nCorrBins, CorrMin, CorrMax);

  // Effective Area
  // < 20 vertices
  TH2F h_CI_EA_low_bar("h_CI_EA_low_bar","Effective Area Correction vs Uncorrected Absolute Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_CI_EA_low_cap("h_CI_EA_low_cap","Effective Area Correction vs Uncorrected Absolute Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax, nCorrBins, CorrMin, CorrMax);
  // 20 < vertices < 40
  TH2F h_CI_EA_med_bar("h_CI_EA_med_bar","Effective Area Correction vs Uncorrected Absolute Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_CI_EA_med_cap("h_CI_EA_med_cap","Effective Area Correction vs Uncorrected Absolute Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax, nCorrBins, CorrMin, CorrMax);
  // > 40 vertices
  TH2F h_CI_EA_high_bar("h_CI_EA_high_bar","Effective Area Correction vs Uncorrected Absolute Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_CI_EA_high_cap("h_CI_EA_high_cap","Effective Area Correction vs Uncorrected Absolute Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax, nCorrBins, CorrMin, CorrMax);

  // Difference
  TH1F h_diffCorr_low_bar("h_diffCorr_low_bar", "#Delta#beta Correction - Effective Area Correction, Absolute",nDiffBins,DiffMin,DiffMax);
  TH1F h_diffCorr_low_cap("h_diffCorr_low_cap", "#Delta#beta Correction - Effective Area Correction, Absolute",nDiffBins,DiffMin,DiffMax);
  TH1F h_diffCorr_med_bar("h_diffCorr_med_bar", "#Delta#beta Correction - Effective Area Correction, Absolute",nDiffBins,DiffMin,DiffMax);
  TH1F h_diffCorr_med_cap("h_diffCorr_med_cap", "#Delta#beta Correction - Effective Area Correction, Absolute",nDiffBins,DiffMin,DiffMax);
  TH1F h_diffCorr_high_bar("h_diffCorr_high_bar", "#Delta#beta Correction - Effective Area Correction, Absolute",nDiffBins,DiffMin,DiffMax);
  TH1F h_diffCorr_high_cap("h_diffCorr_high_cap", "#Delta#beta Correction - Effective Area Correction, Absolute",nDiffBins,DiffMin,DiffMax);

  // Neutral iso
  TH1F h_neut_low_bar("h_neut_low_bar", "Mini Absolute Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_low_cap("h_neut_low_cap", "Mini Absolute Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_med_bar("h_neut_med_bar", "Mini Absolute Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_med_cap("h_neut_med_cap", "Mini Absolute Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_high_bar("h_neut_high_bar", "Mini Absolute Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_high_cap("h_neut_high_cap", "Mini Absolute Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);


  // 2D Comparison
  TH2F h_dB_vs_EA_low_bar("h_dB_vs_EA_low_bar","#Delta#beta Correction vs Effective Area Correction", nCorrBins, CorrMin, CorrMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_dB_vs_EA_low_cap("h_dB_vs_EA_low_cap","#Delta#beta Correction vs Effective Area Correction", nCorrBins, CorrMin, CorrMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_dB_vs_EA_med_bar("h_dB_vs_EA_med_bar","#Delta#beta Correction vs Effective Area Correction", nCorrBins, CorrMin, CorrMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_dB_vs_EA_med_cap("h_dB_vs_EA_med_cap","#Delta#beta Correction vs Effective Area Correction", nCorrBins, CorrMin, CorrMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_dB_vs_EA_high_bar("h_dB_vs_EA_high_bar","#Delta#beta Correction vs Effective Area Correction", nCorrBins, CorrMin, CorrMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_dB_vs_EA_high_cap("h_dB_vs_EA_high_cap","#Delta#beta Correction vs Effective Area Correction", nCorrBins, CorrMin, CorrMax, nCorrBins, CorrMin, CorrMax);

  // Relative Iso

  TH1F h_all_EAcorrs_rel("h_all_EAcorrs_rel","All Effective Area Corrections",nRelCorrsBins,RelCorrsMin,RelCorrsMax);
  TH1F h_all_dBcorrs_rel("h_all_dBcorrs_rel","All #Delta#beta Corrections",nRelCorrsBins,RelCorrsMin,RelCorrsMax);

  vector<vector<float> > EA_corrs_byNpv_rel(13,temp);
  vector<vector<float> > dB_corrs_byNpv_rel(13,temp);
  TH1F h_median_EA_rel("h_median_EA_rel","Median Effective Area Correction by Number of Pileup Vertices",12,0,60);
  TH1F h_median_dB_rel("h_median_dB_rel","Median #Delta#beta Correction by Number of Pileup Vertices",12,0,60);

  TH2F h_CIrel_dB_low_bar("h_CIrel_dB_low_bar","#Delta#beta Correction vs Uncorrected Relative Isolation", nRelIsoBins, RelIsoMin, RelIsoMax, nCorrRelBins, CorrRelMin, CorrRelMax);
  TH2F h_CIrel_dB_low_cap("h_CIrel_dB_low_cap","#Delta#beta Correction vs Uncorrected Relative Isolation", nRelIsoBins, RelIsoMin, RelIsoMax, nCorrRelBins, CorrRelMin, CorrRelMax);
  // 20 < vertices < 40
  TH2F h_CIrel_dB_med_bar("h_CIrel_dB_med_bar","#Delta#beta Correction vs Uncorrected Relative Isolation", nRelIsoBins, RelIsoMin, RelIsoMax, nCorrRelBins, CorrRelMin, CorrRelMax);
  TH2F h_CIrel_dB_med_cap("h_CIrel_dB_med_cap","#Delta#beta Correction vs Uncorrected Relative Isolation", nRelIsoBins, RelIsoMin, RelIsoMax, nCorrRelBins, CorrRelMin, CorrRelMax);
  // > 40 vertices
  TH2F h_CIrel_dB_high_bar("h_CIrel_dB_high_bar","#Delta#beta Correction vs Uncorrected Relative Isolation", nRelIsoBins, RelIsoMin, RelIsoMax, nCorrRelBins, CorrRelMin, CorrRelMax);
  TH2F h_CIrel_dB_high_cap("h_CIrel_dB_high_cap","#Delta#beta Correction vs Uncorrected Relative Isolation", nRelIsoBins, RelIsoMin, RelIsoMax, nCorrRelBins, CorrRelMin, CorrRelMax);

  // Effective Area
  // < 20 vertices
  TH2F h_CIrel_EA_low_bar("h_CIrel_EA_low_bar","Effective Area Correction vs Uncorrected Relative Isolation", nRelIsoBins, RelIsoMin, RelIsoMax, nCorrRelBins, CorrRelMin, CorrRelMax);
  TH2F h_CIrel_EA_low_cap("h_CIrel_EA_low_cap","Effective Area Correction vs Uncorrected Relative Isolation", nRelIsoBins, RelIsoMin, RelIsoMax, nCorrRelBins, CorrRelMin, CorrRelMax);
  // 20 < vertices < 40
  TH2F h_CIrel_EA_med_bar("h_CIrel_EA_med_bar","Effective Area Correction vs Uncorrected Relative Isolation", nRelIsoBins, RelIsoMin, RelIsoMax, nCorrRelBins, CorrRelMin, CorrRelMax);
  TH2F h_CIrel_EA_med_cap("h_CIrel_EA_med_cap","Effective Area Correction vs Uncorrected Relative Isolation", nRelIsoBins, RelIsoMin, RelIsoMax, nCorrRelBins, CorrRelMin, CorrRelMax);
  // > 40 vertices
  TH2F h_CIrel_EA_high_bar("h_CIrel_EA_high_bar","Effective Area Correction vs Uncorrected Relative Isolation", nRelIsoBins, RelIsoMin, RelIsoMax, nCorrRelBins, CorrRelMin, CorrRelMax);
  TH2F h_CIrel_EA_high_cap("h_CIrel_EA_high_cap","Effective Area Correction vs Uncorrected Relative Isolation", nRelIsoBins, RelIsoMin, RelIsoMax, nCorrRelBins, CorrRelMin, CorrRelMax);

  TH1F h_diffCorrRel_low_bar("h_diffCorrRel_low_bar", "#Delta#beta Correction - Effective Area Correction, Relative",nRelDiffBins,RelDiffMin,RelDiffMax);
  TH1F h_diffCorrRel_low_cap("h_diffCorrRel_low_cap", "#Delta#beta Correction - Effective Area Correction, Relative",nRelDiffBins,RelDiffMin,RelDiffMax);
  TH1F h_diffCorrRel_med_bar("h_diffCorrRel_med_bar", "#Delta#beta Correction - Effective Area Correction, Relative",nRelDiffBins,RelDiffMin,RelDiffMax);
  TH1F h_diffCorrRel_med_cap("h_diffCorrRel_med_cap", "#Delta#beta Correction - Effective Area Correction, Relative",nRelDiffBins,RelDiffMin,RelDiffMax);
  TH1F h_diffCorrRel_high_bar("h_diffCorrRel_high_bar", "#Delta#beta Correction - Effective Area Correction, Relative",nRelDiffBins,RelDiffMin,RelDiffMax);
  TH1F h_diffCorrRel_high_cap("h_diffCorrRel_high_cap", "#Delta#beta Correction - Effective Area Correction, Relative",nRelDiffBins,RelDiffMin,RelDiffMax);

  // Neutral iso
  TH1F h_neut_low_bar("h_neut_low_bar", "Mini Relative Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_low_cap("h_neut_low_cap", "Mini Relative Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_med_bar("h_neut_med_bar", "Mini Relative Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_med_cap("h_neut_med_cap", "Mini Relative Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_high_bar("h_neut_high_bar", "Mini Relative Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_high_cap("h_neut_high_cap", "Mini Relative Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);


  // 2D Comparison
  TH2F h_dB_vs_EA_low_bar_rel("h_dB_vs_EA_low_bar_rel","#Delta#beta Correction vs Effective Area Correction", nCorrBins, CorrMin, CorrMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_dB_vs_EA_low_cap_rel("h_dB_vs_EA_low_cap_rel","#Delta#beta Correction vs Effective Area Correction", nCorrBins, CorrMin, CorrMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_dB_vs_EA_med_bar_rel("h_dB_vs_EA_med_bar_rel","#Delta#beta Correction vs Effective Area Correction", nCorrBins, CorrMin, CorrMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_dB_vs_EA_med_cap_rel("h_dB_vs_EA_med_cap_rel","#Delta#beta Correction vs Effective Area Correction", nCorrBins, CorrMin, CorrMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_dB_vs_EA_high_bar_rel("h_dB_vs_EA_high_bar_rel","#Delta#beta Correction vs Effective Area Correction", nCorrBins, CorrMin, CorrMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_dB_vs_EA_high_cap_rel("h_dB_vs_EA_high_cap_rel","#Delta#beta Correction vs Effective Area Correction", nCorrBins, CorrMin, CorrMax, nCorrBins, CorrMin, CorrMax);

  //////////////////
  // Standard Iso //
  //////////////////
  TH1F h_all_EAcorrs_st("h_all_EAcorrs_st","All Effective Area Corrections",nCorrsBins,CorrsMin,CorrsMax);
  TH1F h_all_dBcorrs_st("h_all_dBcorrs_st","All #Delta#beta Corrections",nCorrsBins,CorrsMin,CorrsMax);

  TH1F h_all_EAs_st("h_all_EAs_st","All Effective Areas",20,0,10);
  TH1F h_all_rhos_st("h_all_rhos_st","All Rhos",20,0,10);

  vector<vector<float> > EA_corrs_byNpv_st(13,temp);
  vector<vector<float> > dB_corrs_byNpv_st(13,temp);
  TH1F h_median_EA_st("h_median_EA_st","Median Effective Area Correction by Number of Pileup Vertices",12,0,60);
  TH1F h_median_dB_st("h_median_dB_st","Median #Delta#beta Correction by Number of Pileup Vertices",12,0,60);

  // Delta-Beta
  // < 20 vertices
  TH2F h_CI_dB_low_bar_st("h_CI_dB_low_bar_st","#Delta#beta Correction vs Uncorrected Absolute Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_CI_dB_low_cap_st("h_CI_dB_low_cap_st","#Delta#beta Correction vs Uncorrected Absolute Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax, nCorrBins, CorrMin, CorrMax);
  // 20 < vertices < 40
  TH2F h_CI_dB_med_bar_st("h_CI_dB_med_bar_st","#Delta#beta Correction vs Uncorrected Absolute Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_CI_dB_med_cap_st("h_CI_dB_med_cap_st","#Delta#beta Correction vs Uncorrected Absolute Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax, nCorrBins, CorrMin, CorrMax);
  // > 40 vertices
  TH2F h_CI_dB_high_bar_st("h_CI_dB_high_bar_st","#Delta#beta Correction vs Uncorrected Absolute Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_CI_dB_high_cap_st("h_CI_dB_high_cap_st","#Delta#beta Correction vs Uncorrected Absolute Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax, nCorrBins, CorrMin, CorrMax);

  // Effective Area
  // < 20 vertices
  TH2F h_CI_EA_low_bar_st("h_CI_EA_low_bar_st","Effective Area Correction vs Uncorrected Absolute Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_CI_EA_low_cap_st("h_CI_EA_low_cap_st","Effective Area Correction vs Uncorrected Absolute Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax, nCorrBins, CorrMin, CorrMax);
  // 20 < vertices < 40
  TH2F h_CI_EA_med_bar_st("h_CI_EA_med_bar_st","Effective Area Correction vs Uncorrected Absolute Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_CI_EA_med_cap_st("h_CI_EA_med_cap_st","Effective Area Correction vs Uncorrected Absolute Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax, nCorrBins, CorrMin, CorrMax);
  // > 40 vertices
  TH2F h_CI_EA_high_bar_st("h_CI_EA_high_bar_st","Effective Area Correction vs Uncorrected Absolute Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_CI_EA_high_cap_st("h_CI_EA_high_cap_st","Effective Area Correction vs Uncorrected Absolute Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax, nCorrBins, CorrMin, CorrMax);

  TH1F h_diffCorr_low_bar_st("h_diffCorr_low_bar_st", "#Delta#beta Correction - Effective Area Correction, Absolute",nDiffBins,DiffMin,DiffMax);
  TH1F h_diffCorr_low_cap_st("h_diffCorr_low_cap_st", "#Delta#beta Correction - Effective Area Correction, Absolute",nDiffBins,DiffMin,DiffMax);
  TH1F h_diffCorr_med_bar_st("h_diffCorr_med_bar_st", "#Delta#beta Correction - Effective Area Correction, Absolute",nDiffBins,DiffMin,DiffMax);
  TH1F h_diffCorr_med_cap_st("h_diffCorr_med_cap_st", "#Delta#beta Correction - Effective Area Correction, Absolute",nDiffBins,DiffMin,DiffMax);
  TH1F h_diffCorr_high_bar_st("h_diffCorr_high_bar_st", "#Delta#beta Correction - Effective Area Correction, Absolute",nDiffBins,DiffMin,DiffMax);
  TH1F h_diffCorr_high_cap_st("h_diffCorr_high_cap_st", "#Delta#beta Correction - Effective Area Correction, Absolute",nDiffBins,DiffMin,DiffMax);

  // 2D Comparison
  TH2F h_dB_vs_EA_low_bar_st("h_dB_vs_EA_low_bar_st","#Delta#beta Correction vs Effective Area Correction", nCorrBins, CorrMin, CorrMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_dB_vs_EA_low_cap_st("h_dB_vs_EA_low_cap_st","#Delta#beta Correction vs Effective Area Correction", nCorrBins, CorrMin, CorrMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_dB_vs_EA_med_bar_st("h_dB_vs_EA_med_bar_st","#Delta#beta Correction vs Effective Area Correction", nCorrBins, CorrMin, CorrMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_dB_vs_EA_med_cap_st("h_dB_vs_EA_med_cap_st","#Delta#beta Correction vs Effective Area Correction", nCorrBins, CorrMin, CorrMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_dB_vs_EA_high_bar_st("h_dB_vs_EA_high_bar_st","#Delta#beta Correction vs Effective Area Correction", nCorrBins, CorrMin, CorrMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_dB_vs_EA_high_cap_st("h_dB_vs_EA_high_cap_st","#Delta#beta Correction vs Effective Area Correction", nCorrBins, CorrMin, CorrMax, nCorrBins, CorrMin, CorrMax);

  // Relative

  TH1F h_all_EAcorrs_rel_st("h_all_EAcorrs_rel_st","All Effective Area Corrections",nRelCorrsBins,RelCorrsMin,RelCorrsMax);
  TH1F h_all_dBcorrs_rel_st("h_all_dBcorrs_rel_st","All #Delta#beta Corrections",nRelCorrsBins,RelCorrsMin,RelCorrsMax);

  vector<vector<float> > EA_corrs_byNpv_rel_st(13,temp);
  vector<vector<float> > dB_corrs_byNpv_rel_st(13,temp);
  TH1F h_median_EA_rel_st("h_median_EA_rel_st","Median Effective Area Correction by Number of Pileup Vertices",12,0,60);
  TH1F h_median_dB_rel_st("h_median_dB_rel_st","Median #Delta#beta Correction by Number of Pileup Vertices",12,0,60);


  TH2F h_CIrel_dB_low_bar_st("h_CIrel_dB_low_bar_st","#Delta#beta Correction vs Uncorrected Relative Isolation", nRelIsoBins, RelIsoMin, RelIsoMax, nCorrRelBins, CorrRelMin, CorrRelMax);
  TH2F h_CIrel_dB_low_cap_st("h_CIrel_dB_low_cap_st","#Delta#beta Correction vs Uncorrected Relative Isolation", nRelIsoBins, RelIsoMin, RelIsoMax, nCorrRelBins, CorrRelMin, CorrRelMax);
  // 20 < vertices < 40
  TH2F h_CIrel_dB_med_bar_st("h_CIrel_dB_med_bar_st","#Delta#beta Correction vs Uncorrected Relative Isolation", nRelIsoBins, RelIsoMin, RelIsoMax, nCorrRelBins, CorrRelMin, CorrRelMax);
  TH2F h_CIrel_dB_med_cap_st("h_CIrel_dB_med_cap_st","#Delta#beta Correction vs Uncorrected Relative Isolation", nRelIsoBins, RelIsoMin, RelIsoMax, nCorrRelBins, CorrRelMin, CorrRelMax);
  // > 40 vertices
  TH2F h_CIrel_dB_high_bar_st("h_CIrel_dB_high_bar_st","#Delta#beta Correction vs Uncorrected Relative Isolation", nRelIsoBins, RelIsoMin, RelIsoMax, nCorrRelBins, CorrRelMin, CorrRelMax);
  TH2F h_CIrel_dB_high_cap_st("h_CIrel_dB_high_cap_st","#Delta#beta Correction vs Uncorrected Relative Isolation", nRelIsoBins, RelIsoMin, RelIsoMax, nCorrRelBins, CorrRelMin, CorrRelMax);

  // Effective Area
  // < 20 vertices
  TH2F h_CIrel_EA_low_bar_st("h_CIrel_EA_low_bar_st","Effective Area Correction vs Uncorrected Relative Isolation", nRelIsoBins, RelIsoMin, RelIsoMax, nCorrRelBins, CorrRelMin, CorrRelMax);
  TH2F h_CIrel_EA_low_cap_st("h_CIrel_EA_low_cap_st","Effective Area Correction vs Uncorrected Relative Isolation", nRelIsoBins, RelIsoMin, RelIsoMax, nCorrRelBins, CorrRelMin, CorrRelMax);
  // 20 < vertices < 40
  TH2F h_CIrel_EA_med_bar_st("h_CIrel_EA_med_bar_st","Effective Area Correction vs Uncorrected Relative Isolation", nRelIsoBins, RelIsoMin, RelIsoMax, nCorrRelBins, CorrRelMin, CorrRelMax);
  TH2F h_CIrel_EA_med_cap_st("h_CIrel_EA_med_cap_st","Effective Area Correction vs Uncorrected Relative Isolation", nRelIsoBins, RelIsoMin, RelIsoMax, nCorrRelBins, CorrRelMin, CorrRelMax);
  // > 40 vertices
  TH2F h_CIrel_EA_high_bar_st("h_CIrel_EA_high_bar_st","Effective Area Correction vs Uncorrected Relative Isolation", nRelIsoBins, RelIsoMin, RelIsoMax, nCorrRelBins, CorrRelMin, CorrRelMax);
  TH2F h_CIrel_EA_high_cap_st("h_CIrel_EA_high_cap_st","Effective Area Correction vs Uncorrected Relative Isolation", nRelIsoBins, RelIsoMin, RelIsoMax, nCorrRelBins, CorrRelMin, CorrRelMax);

  TH1F h_diffCorrRel_low_bar_st("h_diffCorrRel_low_bar_st", "#Delta#beta Correction - Effective Area Correction, Relative",nRelDiffBins,RelDiffMin,RelDiffMax);
  TH1F h_diffCorrRel_low_cap_st("h_diffCorrRel_low_cap_st", "#Delta#beta Correction - Effective Area Correction, Relative",nRelDiffBins,RelDiffMin,RelDiffMax);
  TH1F h_diffCorrRel_med_bar_st("h_diffCorrRel_med_bar_st", "#Delta#beta Correction - Effective Area Correction, Relative",nRelDiffBins,RelDiffMin,RelDiffMax);
  TH1F h_diffCorrRel_med_cap_st("h_diffCorrRel_med_cap_st", "#Delta#beta Correction - Effective Area Correction, Relative",nRelDiffBins,RelDiffMin,RelDiffMax);
  TH1F h_diffCorrRel_high_bar_st("h_diffCorrRel_high_bar_st", "#Delta#beta Correction - Effective Area Correction, Relative",nRelDiffBins,RelDiffMin,RelDiffMax);
  TH1F h_diffCorrRel_high_cap_st("h_diffCorrRel_high_cap_st", "#Delta#beta Correction - Effective Area Correction, Relative",nRelDiffBins,RelDiffMin,RelDiffMax);

  // 2D Comparison
  TH2F h_dB_vs_EA_low_bar_rel_st("h_dB_vs_EA_low_bar_rel_st","#Delta#beta Correction vs Effective Area Correction", nCorrBins, CorrMin, CorrMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_dB_vs_EA_low_cap_rel_st("h_dB_vs_EA_low_cap_rel_st","#Delta#beta Correction vs Effective Area Correction", nCorrBins, CorrMin, CorrMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_dB_vs_EA_med_bar_rel_st("h_dB_vs_EA_med_bar_rel_st","#Delta#beta Correction vs Effective Area Correction", nCorrBins, CorrMin, CorrMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_dB_vs_EA_med_cap_rel_st("h_dB_vs_EA_med_cap_rel_st","#Delta#beta Correction vs Effective Area Correction", nCorrBins, CorrMin, CorrMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_dB_vs_EA_high_bar_rel_st("h_dB_vs_EA_high_bar_rel_st","#Delta#beta Correction vs Effective Area Correction", nCorrBins, CorrMin, CorrMax, nCorrBins, CorrMin, CorrMax);
  TH2F h_dB_vs_EA_high_cap_rel_st("h_dB_vs_EA_high_cap_rel_st","#Delta#beta Correction vs Effective Area Correction", nCorrBins, CorrMin, CorrMax, nCorrBins, CorrMin, CorrMax);

  const int nEvents = tree->GetEntries();
  const int event_max = MaxEvents < 0 ? nEvents : min(MaxEvents,nEvents);

  cout << "Running on " << event_max << " events out of " << nEvents << endl;

  //  TTreeCache::SetLearnEntries(10);
  //tree->SetCacheSize(128*1024*1024);
  cms3.Init(tree);
  for (int event = 0; event < event_max; event++) {
    tree->LoadTree(event);
    cms3.GetEntry(event);
    //    CMS3::progress(event, nEvents);

    const int nPUvertices = cms3.puInfo_nPUvertices().at(0);
    const float rho_ctr = cms3.evt_fixgridfastjet_centralneutral_rho();
    const float rho_all = cms3.evt_fixgridfastjet_all_rho();
    //const float w_ = cms3.evt_scale1fb();
    const float w_ = 1;

    // Loop over Reco muons
    for (unsigned int i = 0; i < cms3.mus_p4().size(); i++){
      const float cand_pt = cms3.mus_p4().at(i).pt();
      if (cand_pt < 5) continue; // only use cands with pt > 5 GeV
      const float cand_eta = cms3.mus_p4().at(i).eta();
      if (fabs(cand_eta) > 2.4) continue;
      const float cand_phi = cms3.mus_p4().at(i).phi();

      // make sure there's a gen match
      bool matched = false;
      for (unsigned int j = 0; j < cms3.genps_p4().size() && !matched; j++){
	if (abs(cms3.genps_id().at(j)) != 13) continue;
	const float gen_eta = cms3.genps_p4().at(j).eta();
	const float gen_phi = cms3.genps_p4().at(j).phi();
	matched = DeltaR(cand_eta,gen_eta,cand_phi,gen_phi) < 0.3;
      }
      if (!matched) continue;

      // See IsolationTools.cc for definitions of all these functions

      // at this point we have a Reco muon with pT > 5 and |eta| < 2.4
      const float absiso_uncorr = mus_miniIso_ch().at(i) + mus_miniIso_nh().at(i) + mus_miniIso_em().at(i);
      const float absiso_uncorr_st = mus_isoR03_pf_ChargedHadronPt().at(i) + mus_isoR03_pf_NeutralHadronEt().at(i) + mus_isoR03_pf_PhotonEt().at(i);

      const float dr = getMiniDR(cand_pt);
      const float effective_area = muEA03(i,1) * (dr/0.3) * (dr/0.3);
      const float correction_EA = rho_ctr * effective_area;
      const float effective_area_st = muEA03(i,1);
      const float correction_EA_st = rho_all * effective_area_st;

      h_all_EAcorrs.Fill(correction_EA, w_);
      h_all_EAcorrs_rel.Fill(correction_EA/cand_pt, w_);
      h_all_EAs.Fill(effective_area,w_);
      h_all_rhos.Fill(rho_ctr,w_);
      h_all_EAcorrs_st.Fill(correction_EA_st, w_);
      h_all_EAs_st.Fill(effective_area_st,w_);
      h_all_EAcorrs_rel_st.Fill(correction_EA_st/cand_pt, w_);
      h_all_rhos.Fill(rho_all,w_);


      const float correction_dB = 0.5 * mus_miniIso_db().at(i);
      const float correction_dB_st = 0.5 * mus_isoR03_pf_PUPt().at(i);

      h_all_dBcorrs.Fill(correction_dB, w_);
      h_all_dBcorrs_st.Fill(correction_dB_st, w_);
      h_all_dBcorrs_rel.Fill(correction_dB/cand_pt, w_);
      h_all_dBcorrs_rel_st.Fill(correction_dB_st/cand_pt, w_);

      int PU_index = nPUvertices / 5;
      if (PU_index > 12) PU_index = 12;
      nEntriesByNpv[PU_index]++;
      EA_corrs_byNpv.at(PU_index).push_back(correction_EA);
      dB_corrs_byNpv.at(PU_index).push_back(correction_dB);
      EA_corrs_byNpv_st.at(PU_index).push_back(correction_EA_st);
      dB_corrs_byNpv_st.at(PU_index).push_back(correction_dB_st);
      EA_corrs_byNpv_rel.at(PU_index).push_back(correction_EA/cand_pt);
      dB_corrs_byNpv_rel.at(PU_index).push_back(correction_dB/cand_pt);
      EA_corrs_byNpv_rel_st.at(PU_index).push_back(correction_EA_st/cand_pt);
      dB_corrs_byNpv_rel_st.at(PU_index).push_back(correction_dB_st/cand_pt);

      const float diffCorr = correction_dB - correction_EA;
      const float diffCorr_st = correction_dB_st - correction_EA_st;

      if (nPUvertices < 20){
	if (fabs(cand_eta) < 1.4) {
	  h_CI_dB_low_bar.Fill(absiso_uncorr, correction_dB, w_);
	  h_CI_EA_low_bar.Fill(absiso_uncorr, correction_EA, w_);
	  h_diffCorr_low_bar.Fill(diffCorr);
	  h_CIrel_dB_low_bar.Fill(absiso_uncorr/cand_pt, correction_dB/cand_pt, w_);
	  h_CIrel_EA_low_bar.Fill(absiso_uncorr/cand_pt, correction_EA/cand_pt, w_);
	  h_diffCorrRel_low_bar.Fill(diffCorr/cand_pt);
	  h_dB_vs_EA_low_bar.Fill(correction_EA, correction_dB, w_);
	  h_dB_vs_EA_low_bar_rel.Fill(correction_EA/cand_pt, correction_dB/cand_pt, w_);

	  h_CI_dB_low_bar_st.Fill(absiso_uncorr_st, correction_dB_st, w_);
	  h_CI_EA_low_bar_st.Fill(absiso_uncorr_st, correction_EA_st, w_);
	  h_diffCorr_low_bar_st.Fill(diffCorr_st);
	  h_CIrel_dB_low_bar_st.Fill(absiso_uncorr_st/cand_pt, correction_dB_st/cand_pt, w_);
	  h_CIrel_EA_low_bar_st.Fill(absiso_uncorr_st/cand_pt, correction_EA_st/cand_pt, w_);
	  h_diffCorrRel_low_bar_st.Fill(diffCorr_st/cand_pt);
	  h_dB_vs_EA_low_bar_st.Fill(correction_EA_st, correction_dB_st, w_);
	  h_dB_vs_EA_low_bar_rel_st.Fill(correction_EA_st/cand_pt, correction_dB_st/cand_pt, w_);
	} else {
	  h_CI_dB_low_cap.Fill(absiso_uncorr, correction_dB, w_);
	  h_CI_EA_low_cap.Fill(absiso_uncorr, correction_EA, w_);
	  h_diffCorr_low_cap.Fill(diffCorr, w_);
	  h_CIrel_dB_low_cap.Fill(absiso_uncorr/cand_pt, correction_dB/cand_pt, w_);
	  h_CIrel_EA_low_cap.Fill(absiso_uncorr/cand_pt, correction_EA/cand_pt, w_);
	  h_diffCorrRel_low_cap.Fill(diffCorr/cand_pt, w_);
	  h_dB_vs_EA_low_cap.Fill(correction_EA, correction_dB, w_);
	  h_dB_vs_EA_low_cap_rel.Fill(correction_EA/cand_pt, correction_dB/cand_pt, w_);

	  h_CI_dB_low_cap_st.Fill(absiso_uncorr_st, correction_dB_st, w_);
	  h_CI_EA_low_cap_st.Fill(absiso_uncorr_st, correction_EA_st, w_);
	  h_diffCorr_low_cap_st.Fill(diffCorr_st, w_);
	  h_CIrel_dB_low_cap_st.Fill(absiso_uncorr_st/cand_pt, correction_dB_st/cand_pt, w_);
	  h_CIrel_EA_low_cap_st.Fill(absiso_uncorr_st/cand_pt, correction_EA_st/cand_pt, w_);
	  h_diffCorrRel_low_cap_st.Fill(diffCorr_st/cand_pt, w_);
	  h_dB_vs_EA_low_cap_st.Fill(correction_EA_st, correction_dB_st, w_);
	  h_dB_vs_EA_low_cap_rel_st.Fill(correction_EA_st/cand_pt, correction_dB_st/cand_pt, w_);
	}
      } else if (nPUvertices < 40) {
	if (fabs(cand_eta) < 1.4) {
	  h_CI_dB_med_bar.Fill(absiso_uncorr, correction_dB, w_);
	  h_CI_EA_med_bar.Fill(absiso_uncorr, correction_EA, w_);
	  h_diffCorr_med_bar.Fill(diffCorr, w_);
	  h_CIrel_dB_med_bar.Fill(absiso_uncorr/cand_pt, correction_dB/cand_pt, w_);
	  h_CIrel_EA_med_bar.Fill(absiso_uncorr/cand_pt, correction_EA/cand_pt, w_);
	  h_diffCorrRel_med_bar.Fill(diffCorr/cand_pt, w_);
	  h_dB_vs_EA_med_bar.Fill(correction_EA, correction_dB, w_);
	  h_dB_vs_EA_med_bar_rel.Fill(correction_EA/cand_pt, correction_dB/cand_pt, w_);

	  h_CI_dB_med_bar_st.Fill(absiso_uncorr_st, correction_dB_st, w_);
	  h_CI_EA_med_bar_st.Fill(absiso_uncorr_st, correction_EA_st, w_);
	  h_diffCorr_med_bar_st.Fill(diffCorr_st, w_);
	  h_CIrel_dB_med_bar_st.Fill(absiso_uncorr_st/cand_pt, correction_dB_st/cand_pt, w_);
	  h_CIrel_EA_med_bar_st.Fill(absiso_uncorr_st/cand_pt, correction_EA_st/cand_pt, w_);
	  h_diffCorrRel_med_bar_st.Fill(diffCorr_st/cand_pt, w_);
	  h_dB_vs_EA_med_bar_st.Fill(correction_EA_st, correction_dB_st, w_);
	  h_dB_vs_EA_med_bar_rel_st.Fill(correction_EA/cand_pt, correction_dB_st/cand_pt, w_);
	} else {
	  h_CI_dB_med_cap.Fill(absiso_uncorr, correction_dB, w_);
	  h_CI_EA_med_cap.Fill(absiso_uncorr, correction_EA, w_);
	  h_diffCorr_med_cap.Fill(diffCorr, w_);
	  h_CIrel_dB_med_cap.Fill(absiso_uncorr/cand_pt, correction_dB/cand_pt, w_);
	  h_CIrel_EA_med_cap.Fill(absiso_uncorr/cand_pt, correction_EA/cand_pt, w_);
	  h_diffCorrRel_med_cap.Fill(diffCorr/cand_pt, w_);
	  h_dB_vs_EA_med_cap.Fill(correction_EA, correction_dB, w_);
	  h_dB_vs_EA_med_cap_rel.Fill(correction_EA/cand_pt, correction_dB/cand_pt, w_);

	  h_CI_dB_med_cap_st.Fill(absiso_uncorr_st, correction_dB_st, w_);
	  h_CI_EA_med_cap_st.Fill(absiso_uncorr_st, correction_EA_st, w_);
	  h_diffCorr_med_cap_st.Fill(diffCorr_st, w_);
	  h_CIrel_dB_med_cap_st.Fill(absiso_uncorr_st/cand_pt, correction_dB_st/cand_pt, w_);
	  h_CIrel_EA_med_cap_st.Fill(absiso_uncorr_st/cand_pt, correction_EA_st/cand_pt, w_);
	  h_diffCorrRel_med_cap_st.Fill(diffCorr_st/cand_pt, w_);
	  h_dB_vs_EA_med_cap_st.Fill(correction_EA_st, correction_dB_st, w_);
	  h_dB_vs_EA_med_cap_rel_st.Fill(correction_EA_st/cand_pt, correction_dB_st/cand_pt, w_);
	}
      }	else {
	if (fabs(cand_eta) < 1.4) {
	  h_CI_dB_high_bar.Fill(absiso_uncorr, correction_dB, w_);
	  h_CI_EA_high_bar.Fill(absiso_uncorr, correction_EA, w_);
	  h_diffCorr_high_bar.Fill(diffCorr, w_);
	  h_CIrel_dB_high_bar.Fill(absiso_uncorr/cand_pt, correction_dB/cand_pt, w_);
	  h_CIrel_EA_high_bar.Fill(absiso_uncorr/cand_pt, correction_EA/cand_pt, w_);
	  h_diffCorrRel_high_bar.Fill(diffCorr/cand_pt, w_);
	  h_dB_vs_EA_high_bar.Fill(correction_EA, correction_dB, w_);
	  h_dB_vs_EA_high_bar_rel.Fill(correction_EA/cand_pt, correction_dB/cand_pt, w_);

	  h_CI_dB_high_bar_st.Fill(absiso_uncorr_st, correction_dB_st, w_);
	  h_CI_EA_high_bar_st.Fill(absiso_uncorr_st, correction_EA_st, w_);
	  h_diffCorr_high_bar_st.Fill(diffCorr_st, w_);
	  h_CIrel_dB_high_bar_st.Fill(absiso_uncorr_st/cand_pt, correction_dB_st/cand_pt, w_);
	  h_CIrel_EA_high_bar_st.Fill(absiso_uncorr_st/cand_pt, correction_EA_st/cand_pt, w_);
	  h_diffCorrRel_high_bar_st.Fill(diffCorr_st/cand_pt, w_);
	  h_dB_vs_EA_high_bar_st.Fill(correction_EA_st, correction_dB_st, w_);
	  h_dB_vs_EA_high_bar_rel_st.Fill(correction_EA_st/cand_pt, correction_dB_st/cand_pt, w_);
	} else {
	  h_CI_dB_high_cap.Fill(absiso_uncorr, correction_dB);
	  h_CI_EA_high_cap.Fill(absiso_uncorr, correction_EA);
	  h_diffCorr_high_cap.Fill(diffCorr, w_);
	  h_CIrel_dB_high_cap.Fill(absiso_uncorr/cand_pt, correction_dB/cand_pt);
	  h_CIrel_EA_high_cap.Fill(absiso_uncorr/cand_pt, correction_EA/cand_pt);
	  h_diffCorrRel_high_cap.Fill(diffCorr/cand_pt, w_);
	  h_dB_vs_EA_high_cap.Fill(correction_EA, correction_dB, w_);
	  h_dB_vs_EA_high_cap_rel.Fill(correction_EA/cand_pt, correction_dB/cand_pt, w_);

	  h_CI_dB_high_cap_st.Fill(absiso_uncorr_st, correction_dB_st);
	  h_CI_EA_high_cap_st.Fill(absiso_uncorr_st, correction_EA_st);
	  h_diffCorr_high_cap_st.Fill(diffCorr_st, w_);
	  h_CIrel_dB_high_cap_st.Fill(absiso_uncorr_st/cand_pt, correction_dB_st/cand_pt);
	  h_CIrel_EA_high_cap_st.Fill(absiso_uncorr_st/cand_pt, correction_EA_st/cand_pt);
	  h_diffCorrRel_high_cap_st.Fill(diffCorr_st/cand_pt, w_);
	  h_dB_vs_EA_high_cap_st.Fill(correction_EA_st, correction_dB_st, w_);
	  h_dB_vs_EA_high_cap_rel_st.Fill(correction_EA_st/cand_pt, correction_dB_st/cand_pt, w_);

	}
      } // Histogram filling block
    } // Muon loop

  } // Event Loop

  // get median corrections as function of Npv
  for (int i = 0; i < 13; i++) {
    sort(EA_corrs_byNpv.at(i).begin(),EA_corrs_byNpv.at(i).begin()+nEntriesByNpv[i]);
    sort(dB_corrs_byNpv.at(i).begin(),dB_corrs_byNpv.at(i).begin()+nEntriesByNpv[i]);
    sort(EA_corrs_byNpv_st.at(i).begin(),EA_corrs_byNpv_st.at(i).begin()+nEntriesByNpv[i]);
    sort(dB_corrs_byNpv_st.at(i).begin(),dB_corrs_byNpv_st.at(i).begin()+nEntriesByNpv[i]);
    sort(EA_corrs_byNpv_rel.at(i).begin(),EA_corrs_byNpv_rel.at(i).begin()+nEntriesByNpv[i]);
    sort(dB_corrs_byNpv_rel.at(i).begin(),dB_corrs_byNpv_rel.at(i).begin()+nEntriesByNpv[i]);
    sort(EA_corrs_byNpv_rel_st.at(i).begin(),EA_corrs_byNpv_rel_st.at(i).begin()+nEntriesByNpv[i]);
    sort(dB_corrs_byNpv_rel_st.at(i).begin(),dB_corrs_byNpv_rel_st.at(i).begin()+nEntriesByNpv[i]);
    float EAmedian = 0, dBmedian = 0;
    float EAmedian_st = 0, dBmedian_st = 0;
    float EAmedian_rel = 0, dBmedian_rel = 0;
    float EAmedian_rel_st = 0, dBmedian_rel_st = 0;
    if (nEntriesByNpv[i] > 0) {
      if (nEntriesByNpv[i] % 2) {
	// Odd
	EAmedian = EA_corrs_byNpv.at(i).at( (nEntriesByNpv[i] / 2) );
	dBmedian = dB_corrs_byNpv.at(i).at( (nEntriesByNpv[i] / 2) );
	EAmedian_st = EA_corrs_byNpv_st.at(i).at( (nEntriesByNpv[i] / 2) );
	dBmedian_st = dB_corrs_byNpv_st.at(i).at( (nEntriesByNpv[i] / 2) );
	EAmedian_rel = EA_corrs_byNpv_rel.at(i).at( (nEntriesByNpv[i] / 2) );
	dBmedian_rel = dB_corrs_byNpv_rel.at(i).at( (nEntriesByNpv[i] / 2) );
	EAmedian_rel_st = EA_corrs_byNpv_rel_st.at(i).at( (nEntriesByNpv[i] / 2) );
	dBmedian_rel_st = dB_corrs_byNpv_rel_st.at(i).at( (nEntriesByNpv[i] / 2) );
      }
      else {
	// Even
	EAmedian = 0.5 * (EA_corrs_byNpv.at(i).at( (nEntriesByNpv[i] / 2) ) + EA_corrs_byNpv.at(i).at( (nEntriesByNpv[i] / 2) - 1));
	dBmedian = 0.5 * (dB_corrs_byNpv.at(i).at( (nEntriesByNpv[i] / 2) ) + dB_corrs_byNpv.at(i).at( (nEntriesByNpv[i] / 2) - 1));
	EAmedian_st = 0.5 * (EA_corrs_byNpv_st.at(i).at( (nEntriesByNpv[i] / 2) ) + EA_corrs_byNpv_st.at(i).at( (nEntriesByNpv[i] / 2) - 1));
	dBmedian_st = 0.5 * (dB_corrs_byNpv_st.at(i).at( (nEntriesByNpv[i] / 2) ) + dB_corrs_byNpv_st.at(i).at( (nEntriesByNpv[i] / 2) - 1));
	EAmedian_rel = 0.5 * (EA_corrs_byNpv_rel.at(i).at( (nEntriesByNpv[i] / 2) ) + EA_corrs_byNpv_rel.at(i).at( (nEntriesByNpv[i] / 2) - 1));
	dBmedian_rel = 0.5 * (dB_corrs_byNpv_rel.at(i).at( (nEntriesByNpv[i] / 2) ) + dB_corrs_byNpv_rel.at(i).at( (nEntriesByNpv[i] / 2) - 1));
	EAmedian_rel_st = 0.5 * (EA_corrs_byNpv_rel_st.at(i).at( (nEntriesByNpv[i] / 2) ) + EA_corrs_byNpv_rel_st.at(i).at( (nEntriesByNpv[i] / 2) - 1));
	dBmedian_rel_st = 0.5 * (dB_corrs_byNpv_rel_st.at(i).at( (nEntriesByNpv[i] / 2) ) + dB_corrs_byNpv_rel_st.at(i).at( (nEntriesByNpv[i] / 2) - 1));
      }
    }
    h_median_EA.SetBinContent(i,EAmedian);
    h_median_EA_st.SetBinContent(i,EAmedian_st);
    h_median_dB.SetBinContent(i,dBmedian);
    h_median_dB_st.SetBinContent(i,dBmedian_st);
    h_median_EA_rel.SetBinContent(i,EAmedian_rel);
    h_median_EA_rel_st.SetBinContent(i,EAmedian_rel_st);
    h_median_dB_rel.SetBinContent(i,dBmedian_rel);
    h_median_dB_rel_st.SetBinContent(i,dBmedian_rel_st);
    if (nEntriesByNpv[i] > 0) {
      h_median_EA.SetBinError(i,EAmedian/sqrt(nEntriesByNpv[i]));
      h_median_EA_st.SetBinError(i,EAmedian_st/sqrt(nEntriesByNpv[i]));
      h_median_dB.SetBinError(i,dBmedian/sqrt(nEntriesByNpv[i]));
      h_median_dB_st.SetBinError(i,dBmedian_st/sqrt(nEntriesByNpv[i]));
      h_median_EA_rel.SetBinError(i,EAmedian_rel/sqrt(nEntriesByNpv[i]));
      h_median_EA_rel_st.SetBinError(i,EAmedian_rel_st/sqrt(nEntriesByNpv[i]));
      h_median_dB_rel.SetBinError(i,dBmedian_rel/sqrt(nEntriesByNpv[i]));
      h_median_dB_rel_st.SetBinError(i,dBmedian_rel_st/sqrt(nEntriesByNpv[i]));
    }
  }


  TFile * OutFile_ = new TFile(Form("%s.root", outname), "RECREATE");
  OutFile_->cd();
  h_all_EAs.Write();
  h_all_rhos.Write();
  h_median_EA.Write();
  h_median_dB.Write();
  h_all_EAcorrs.Write();
  h_all_dBcorrs.Write();
  h_median_EA_rel.Write();
  h_median_dB_rel.Write();
  h_all_EAcorrs_rel.Write();
  h_all_dBcorrs_rel.Write();

  h_CI_dB_low_bar.Write();
  h_CI_dB_low_cap.Write();
  h_CI_dB_med_bar.Write();
  h_CI_dB_med_cap.Write();
  h_CI_dB_high_bar.Write();
  h_CI_dB_high_cap.Write();
  h_CI_EA_low_bar.Write();
  h_CI_EA_low_cap.Write();
  h_CI_EA_med_bar.Write();
  h_CI_EA_med_cap.Write();
  h_CI_EA_high_bar.Write();
  h_CI_EA_high_cap.Write();
  h_diffCorr_low_bar.Write();
  h_diffCorr_low_cap.Write();
  h_diffCorr_med_bar.Write();
  h_diffCorr_med_cap.Write();
  h_diffCorr_high_bar.Write();
  h_diffCorr_high_cap.Write();
  h_dB_vs_EA_low_bar.Write();
  h_dB_vs_EA_low_cap.Write();
  h_dB_vs_EA_med_bar.Write();
  h_dB_vs_EA_med_cap.Write();
  h_dB_vs_EA_high_bar.Write();
  h_dB_vs_EA_high_cap.Write();

  h_CIrel_dB_low_bar.Write();
  h_CIrel_dB_low_cap.Write();
  h_CIrel_dB_med_bar.Write();
  h_CIrel_dB_med_cap.Write();
  h_CIrel_dB_high_bar.Write();
  h_CIrel_dB_high_cap.Write();
  h_CIrel_EA_low_bar.Write();
  h_CIrel_EA_low_cap.Write();
  h_CIrel_EA_med_bar.Write();
  h_CIrel_EA_med_cap.Write();
  h_CIrel_EA_high_bar.Write();
  h_CIrel_EA_high_cap.Write();
  h_diffCorrRel_low_bar.Write();
  h_diffCorrRel_low_cap.Write();
  h_diffCorrRel_med_bar.Write();
  h_diffCorrRel_med_cap.Write();
  h_diffCorrRel_high_bar.Write();
  h_diffCorrRel_high_cap.Write();
  h_dB_vs_EA_low_bar_rel.Write();
  h_dB_vs_EA_low_cap_rel.Write();
  h_dB_vs_EA_med_bar_rel.Write();
  h_dB_vs_EA_med_cap_rel.Write();
  h_dB_vs_EA_high_bar_rel.Write();
  h_dB_vs_EA_high_cap_rel.Write();

  h_all_EAs_st.Write();
  h_all_rhos_st.Write();
  h_median_EA_st.Write();
  h_median_dB_st.Write();
  h_all_EAcorrs_st.Write();
  h_all_dBcorrs_st.Write();
  h_median_EA_rel_st.Write();
  h_median_dB_rel_st.Write();
  h_all_EAcorrs_rel_st.Write();
  h_all_dBcorrs_rel_st.Write();

  h_CI_dB_low_bar_st.Write();
  h_CI_dB_low_cap_st.Write();
  h_CI_dB_med_bar_st.Write();
  h_CI_dB_med_cap_st.Write();
  h_CI_dB_high_bar_st.Write();
  h_CI_dB_high_cap_st.Write();
  h_CI_EA_low_bar_st.Write();
  h_CI_EA_low_cap_st.Write();
  h_CI_EA_med_bar_st.Write();
  h_CI_EA_med_cap_st.Write();
  h_CI_EA_high_bar_st.Write();
  h_CI_EA_high_cap_st.Write();
  h_diffCorr_low_bar_st.Write();
  h_diffCorr_low_cap_st.Write();
  h_diffCorr_med_bar_st.Write();
  h_diffCorr_med_cap_st.Write();
  h_diffCorr_high_bar_st.Write();
  h_diffCorr_high_cap_st.Write();
  h_dB_vs_EA_low_bar_st.Write();
  h_dB_vs_EA_low_cap_st.Write();
  h_dB_vs_EA_med_bar_st.Write();
  h_dB_vs_EA_med_cap_st.Write();
  h_dB_vs_EA_high_bar_st.Write();
  h_dB_vs_EA_high_cap_st.Write();

  h_CIrel_dB_low_bar_st.Write();
  h_CIrel_dB_low_cap_st.Write();
  h_CIrel_dB_med_bar_st.Write();
  h_CIrel_dB_med_cap_st.Write();
  h_CIrel_dB_high_bar_st.Write();
  h_CIrel_dB_high_cap_st.Write();
  h_CIrel_EA_low_bar_st.Write();
  h_CIrel_EA_low_cap_st.Write();
  h_CIrel_EA_med_bar_st.Write();
  h_CIrel_EA_med_cap_st.Write();
  h_CIrel_EA_high_bar_st.Write();
  h_CIrel_EA_high_cap_st.Write();
  h_diffCorrRel_low_bar_st.Write();
  h_diffCorrRel_low_cap_st.Write();
  h_diffCorrRel_med_bar_st.Write();
  h_diffCorrRel_med_cap_st.Write();
  h_diffCorrRel_high_bar_st.Write();
  h_diffCorrRel_high_cap_st.Write();
  h_dB_vs_EA_low_bar_rel_st.Write();
  h_dB_vs_EA_low_cap_rel_st.Write();
  h_dB_vs_EA_med_bar_rel_st.Write();
  h_dB_vs_EA_med_cap_rel_st.Write();
  h_dB_vs_EA_high_bar_rel_st.Write();
  h_dB_vs_EA_high_cap_rel_st.Write();
  OutFile_->Close();
  
}

int main (int argc, char ** argv) {
  if (argc < 3) {
    std::cout << "USAGE: ./ScanChain.exe <tag> <filename> [<max_num_events>]" << std::endl;
    return 1;
  }

  TString outfileid(argv[1]); 
  TString sample(argv[2]); 

  int max_events = -1;
  if (argc >= 4) max_events = atoi(argv[3]);
  std::cout << "set max number of events to: " << max_events << std::endl;
  
  TChain *chain = new TChain("Events");
  TString inputs = Form("%s", sample.Data());
  chain->Add(inputs);
  TTree *tree = (TTree*)chain->Clone("Events");
  delete chain;

  PileupCorrection puc;
  puc.ScanChain(tree,outfileid.Data(),max_events);
  return 0;
}
