// A cms3 looper to compare delta-beta and effective area pileup corrections

#include <iostream>
#include "ScanChain.h"

using namespace std;
using namespace tas;

void PileupCorrection::ScanChain (TTree * tree, const char* outname, int MaxEvents) {
  const int nAbsIsoBins = 25;
  const float AbsIsoMin = 0;
  const float AbsIsoMax = 2;
  // const int nCorrBins = 100;
  // const float CorrMin = 0.001;
  // const float CorrMax = 2;

  const int nRelIsoBins = 25;
  const float RelIsoMin = 0;
  const float RelIsoMax = 1;
  // const int nCorrRelBins = 100;
  // const float CorrRelMin = 0.001;
  // const float CorrRelMax = 1;

  const int nCorrsBins = 24;
  const float CorrsMin = 0;
  const float CorrsMax = 5;
  const int nRelCorrsBins = 24;
  const float RelCorrsMin = 0;
  const float RelCorrsMax = 1;

  const int nDiffBins = 20;
  const int DiffMin = -4;
  const int DiffMax = 4;
  const int nRelDiffBins = 20;
  const float RelDiffMin = -2;
  const float RelDiffMax = 2;


  // TProfiles (standard only)
  TProfile p_pt_EA("p_pt_EA","EA Corrections by p_{T}",20,0,100);
  TProfile p_npv_EA("p_npv_EA","EA Corrections by N_{PV}",12,0,60);
  TProfile p_pt_dB("p_pt_dB","#Delta#beta Corrections by p_{T}",20,0,100);
  TProfile p_npv_dB("p_npv_dB","#Delta#beta Corrections by N_{PV}",12,0,60);
  TProfile p_npv_rho_all("p_npv_rho_all","#rho_{all} by N_{PV}",12,0,60);
  TProfile p_npv_rho_ctr("p_npv_rho_ctr","#rho_{ctr} by N_{PV}",12,0,60);

  TH1F h_npv_over_EA("h_npv_over_EA","Fraction of EA Overcorrections by N_{PV}",12,0,60);
  TH1F h_npv_over_dB("h_npv_over_dB","Fraction of #Delta#beta Overcorrections by N_{PV}",12,0,60);
  // use nEntriesByNpv array defined below as denom
  int nOverEA[13] = {0};
  int nOverdB[13] = {0};

  //////////
  // Mini //
  //////////

  TH1::SetDefaultSumw2(true);

  // Absolute Iso
  TH1F h_all_EAcorrs("h_all_EAcorrs","All Effective Area Corrections",nCorrsBins,CorrsMin,CorrsMax);
  TH1F h_all_dBcorrs("h_all_dBcorrs","All #Delta#beta Corrections",nCorrsBins,CorrsMin,CorrsMax);

  TH1F h_all_EAs("h_all_EAs","All Effective Areas",20,0,10);
  TH1F h_all_rhos_ctr("h_all_rhos_ctr","All Rhos",20,0,10);

  vector<float> temp;
  temp.reserve(1000);
  int nEntriesByNpv[13] = {0};
  vector<vector<float> > EA_corrs_byNpv(13,temp);
  vector<vector<float> > dB_corrs_byNpv(13,temp);
  TH1F h_median_EA("h_median_EA","Median Effective Area Correction by Number of Pileup Vertices",12,0,60);
  TH1F h_median_dB("h_median_dB","Median #Delta#beta Correction by Number of Pileup Vertices",12,0,60);

  // Difference
  TH1F h_diffCorr_low_bar("h_diffCorr_low_bar", "#Delta#beta Correction - Effective Area Correction, Absolute",nDiffBins,DiffMin,DiffMax);
  TH1F h_diffCorr_low_cap("h_diffCorr_low_cap", "#Delta#beta Correction - Effective Area Correction, Absolute",nDiffBins,DiffMin,DiffMax);
  TH1F h_diffCorr_med_bar("h_diffCorr_med_bar", "#Delta#beta Correction - Effective Area Correction, Absolute",nDiffBins,DiffMin,DiffMax);
  TH1F h_diffCorr_med_cap("h_diffCorr_med_cap", "#Delta#beta Correction - Effective Area Correction, Absolute",nDiffBins,DiffMin,DiffMax);
  TH1F h_diffCorr_high_bar("h_diffCorr_high_bar", "#Delta#beta Correction - Effective Area Correction, Absolute",nDiffBins,DiffMin,DiffMax);
  TH1F h_diffCorr_high_cap("h_diffCorr_high_cap", "#Delta#beta Correction - Effective Area Correction, Absolute",nDiffBins,DiffMin,DiffMax);

  // Neutral iso
  // By Npv and Eta
  TH1F h_neut_low_bar("h_neut_low_bar", "Mini Absolute Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_low_cap("h_neut_low_cap", "Mini Absolute Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_med_bar("h_neut_med_bar", "Mini Absolute Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_med_cap("h_neut_med_cap", "Mini Absolute Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_high_bar("h_neut_high_bar", "Mini Absolute Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_high_cap("h_neut_high_cap", "Mini Absolute Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);

  // By Npv and Cand Pt
  TH1F h_neut_soft("h_neut_soft", "Mini Absolute Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_med("h_neut_med", "Mini Absolute Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_hard("h_neut_hard", "Mini Absolute Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);

  // Relative

  TH1F h_all_EAcorrs_rel("h_all_EAcorrs_rel","All Effective Area Corrections",nRelCorrsBins,RelCorrsMin,RelCorrsMax);
  TH1F h_all_dBcorrs_rel("h_all_dBcorrs_rel","All #Delta#beta Corrections",nRelCorrsBins,RelCorrsMin,RelCorrsMax);

  vector<vector<float> > EA_corrs_byNpv_rel(13,temp);
  vector<vector<float> > dB_corrs_byNpv_rel(13,temp);
  TH1F h_median_EA_rel("h_median_EA_rel","Median Effective Area Correction by Number of Pileup Vertices",12,0,60);
  TH1F h_median_dB_rel("h_median_dB_rel","Median #Delta#beta Correction by Number of Pileup Vertices",12,0,60);

  TH1F h_diffCorr_low_bar_rel("h_diffCorr_low_bar_rel", "#Delta#beta Correction - Effective Area Correction, Relative",nRelDiffBins,RelDiffMin,RelDiffMax);
  TH1F h_diffCorr_low_cap_rel("h_diffCorr_low_cap_rel", "#Delta#beta Correction - Effective Area Correction, Relative",nRelDiffBins,RelDiffMin,RelDiffMax);
  TH1F h_diffCorr_med_bar_rel("h_diffCorr_med_bar_rel", "#Delta#beta Correction - Effective Area Correction, Relative",nRelDiffBins,RelDiffMin,RelDiffMax);
  TH1F h_diffCorr_med_cap_rel("h_diffCorr_med_cap_rel", "#Delta#beta Correction - Effective Area Correction, Relative",nRelDiffBins,RelDiffMin,RelDiffMax);
  TH1F h_diffCorr_high_bar_rel("h_diffCorr_high_bar_rel", "#Delta#beta Correction - Effective Area Correction, Relative",nRelDiffBins,RelDiffMin,RelDiffMax);
  TH1F h_diffCorr_high_cap_rel("h_diffCorr_high_cap_rel", "#Delta#beta Correction - Effective Area Correction, Relative",nRelDiffBins,RelDiffMin,RelDiffMax);

  // Neutral iso
  TH1F h_neut_low_bar_rel("h_neut_low_bar_rel", "Mini Relative Neutral Isolation", nRelIsoBins, RelIsoMin, RelIsoMax);
  TH1F h_neut_low_cap_rel("h_neut_low_cap_rel", "Mini Relative Neutral Isolation", nRelIsoBins, RelIsoMin, RelIsoMax);
  TH1F h_neut_med_bar_rel("h_neut_med_bar_rel", "Mini Relative Neutral Isolation", nRelIsoBins, RelIsoMin, RelIsoMax);
  TH1F h_neut_med_cap_rel("h_neut_med_cap_rel", "Mini Relative Neutral Isolation", nRelIsoBins, RelIsoMin, RelIsoMax);
  TH1F h_neut_high_bar_rel("h_neut_high_bar_rel", "Mini Relative Neutral Isolation", nRelIsoBins, RelIsoMin, RelIsoMax);
  TH1F h_neut_high_cap_rel("h_neut_high_cap_rel", "Mini Relative Neutral Isolation", nRelIsoBins, RelIsoMin, RelIsoMax);

  // By Npv and Cand Pt
  TH1F h_neut_soft_rel("h_neut_soft_rel", "Mini Relative Neutral Isolation", nRelIsoBins, RelIsoMin, RelIsoMax);
  TH1F h_neut_med_rel("h_neut_med_rel", "Mini Relative Neutral Isolation", nRelIsoBins, RelIsoMin, RelIsoMax);
  TH1F h_neut_hard_rel("h_neut_hard_rel", "Mini Relative Neutral Isolation", nRelIsoBins, RelIsoMin, RelIsoMax);

  //////////////////
  // Standard Iso //
  //////////////////
  TH1F h_all_EAcorrs_st("h_all_EAcorrs_st","All Effective Area Corrections",nCorrsBins,CorrsMin,CorrsMax);
  TH1F h_all_dBcorrs_st("h_all_dBcorrs_st","All #Delta#beta Corrections",nCorrsBins,CorrsMin,CorrsMax);

  TH1F h_all_EAs_st("h_all_EAs_st","All Effective Areas",20,0,10);
  TH1F h_all_rhos_all("h_all_rhos_all","All Rhos",20,0,20);

  vector<vector<float> > EA_corrs_byNpv_st(13,temp);
  vector<vector<float> > dB_corrs_byNpv_st(13,temp);
  TH1F h_median_EA_st("h_median_EA_st","Median Effective Area Correction by Number of Pileup Vertices",12,0,60);
  TH1F h_median_dB_st("h_median_dB_st","Median #Delta#beta Correction by Number of Pileup Vertices",12,0,60);

  TH1F h_diffCorr_low_bar_st("h_diffCorr_low_bar_st", "#Delta#beta Correction - Effective Area Correction, Absolute",nDiffBins,DiffMin,DiffMax);
  TH1F h_diffCorr_low_cap_st("h_diffCorr_low_cap_st", "#Delta#beta Correction - Effective Area Correction, Absolute",nDiffBins,DiffMin,DiffMax);
  TH1F h_diffCorr_med_bar_st("h_diffCorr_med_bar_st", "#Delta#beta Correction - Effective Area Correction, Absolute",nDiffBins,DiffMin,DiffMax);
  TH1F h_diffCorr_med_cap_st("h_diffCorr_med_cap_st", "#Delta#beta Correction - Effective Area Correction, Absolute",nDiffBins,DiffMin,DiffMax);
  TH1F h_diffCorr_high_bar_st("h_diffCorr_high_bar_st", "#Delta#beta Correction - Effective Area Correction, Absolute",nDiffBins,DiffMin,DiffMax);
  TH1F h_diffCorr_high_cap_st("h_diffCorr_high_cap_st", "#Delta#beta Correction - Effective Area Correction, Absolute",nDiffBins,DiffMin,DiffMax);

  // Neutral iso
  TH1F h_neut_low_bar_st("h_neut_low_bar_st", "#DeltaR = 0.3 Absolute Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_low_cap_st("h_neut_low_cap_st", "#DeltaR = 0.3 Absolute Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_med_bar_st("h_neut_med_bar_st", "#DeltaR = 0.3 Absolute Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_med_cap_st("h_neut_med_cap_st", "#DeltaR = 0.3 Absolute Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_high_bar_st("h_neut_high_bar_st", "#DeltaR = 0.3 Absolute Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_high_cap_st("h_neut_high_cap_st", "#DeltaR = 0.3 Absolute Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);

  // By Npv and Cand Pt
  TH1F h_neut_soft_st("h_neut_soft_st", "#DeltaR=0.3 Absolute Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_med_st("h_neut_med_st", "#DeltaR=0.3 Absolute Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_hard_st("h_neut_hard_st", "#DeltaR=0.3  Absolute Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);


  // Relative

  TH1F h_all_EAcorrs_rel_st("h_all_EAcorrs_rel_st","All Effective Area Corrections",nRelCorrsBins,RelCorrsMin,RelCorrsMax);
  TH1F h_all_dBcorrs_rel_st("h_all_dBcorrs_rel_st","All #Delta#beta Corrections",nRelCorrsBins,RelCorrsMin,RelCorrsMax);

  vector<vector<float> > EA_corrs_byNpv_rel_st(13,temp);
  vector<vector<float> > dB_corrs_byNpv_rel_st(13,temp);
  TH1F h_median_EA_rel_st("h_median_EA_rel_st","Median Effective Area Correction by Number of Pileup Vertices",12,0,60);
  TH1F h_median_dB_rel_st("h_median_dB_rel_st","Median #Delta#beta Correction by Number of Pileup Vertices",12,0,60);

  TH1F h_diffCorr_low_bar_rel_st("h_diffCorr_low_bar_rel_st", "#Delta#beta Correction - Effective Area Correction, Relative",nRelDiffBins,RelDiffMin,RelDiffMax);
  TH1F h_diffCorr_low_cap_rel_st("h_diffCorr_low_cap_rel_st", "#Delta#beta Correction - Effective Area Correction, Relative",nRelDiffBins,RelDiffMin,RelDiffMax);
  TH1F h_diffCorr_med_bar_rel_st("h_diffCorr_med_bar_rel_st", "#Delta#beta Correction - Effective Area Correction, Relative",nRelDiffBins,RelDiffMin,RelDiffMax);
  TH1F h_diffCorr_med_cap_rel_st("h_diffCorr_med_cap_rel_st", "#Delta#beta Correction - Effective Area Correction, Relative",nRelDiffBins,RelDiffMin,RelDiffMax);
  TH1F h_diffCorr_high_bar_rel_st("h_diffCorr_high_bar_rel_st", "#Delta#beta Correction - Effective Area Correction, Relative",nRelDiffBins,RelDiffMin,RelDiffMax);
  TH1F h_diffCorr_high_cap_rel_st("h_diffCorr_high_cap_rel_st", "#Delta#beta Correction - Effective Area Correction, Relative",nRelDiffBins,RelDiffMin,RelDiffMax);

  // Neutral iso
  // By Npv and Eta
  TH1F h_neut_low_bar_rel_st("h_neut_low_bar_rel_st", "#DeltaR = 0.3 Relative Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_low_cap_rel_st("h_neut_low_cap_rel_st", "#DeltaR = 0.3 Relative Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_med_bar_rel_st("h_neut_med_bar_rel_st", "#DeltaR = 0.3 Relative Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_med_cap_rel_st("h_neut_med_cap_rel_st", "#DeltaR = 0.3 Relative Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_high_bar_rel_st("h_neut_high_bar_rel_st", "#DeltaR = 0.3 Relative Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);
  TH1F h_neut_high_cap_rel_st("h_neut_high_cap_rel_st", "#DeltaR = 0.3 Relative Neutral Isolation", nAbsIsoBins, AbsIsoMin, AbsIsoMax);

  // By Cand Pt
  TH1F h_neut_soft_rel_st("h_neut_soft_rel_st", "#DeltaR=0.3 Relative Neutral Isolation", nRelIsoBins, RelIsoMin, RelIsoMax);
  TH1F h_neut_med_rel_st("h_neut_med_rel_st", "#DeltaR=0.3 Relative Neutral Isolation", nRelIsoBins, RelIsoMin, RelIsoMax);
  TH1F h_neut_hard_rel_st("h_neut_hard_rel_st", "#DeltaR=0.3 Relative Neutral Isolation", nRelIsoBins, RelIsoMin, RelIsoMax);

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
      for (unsigned int j = 0; j < cms3.mus_mc_p4().size() && !matched; j++){
	if ( cms3.mus_mc_motherid().at(j) != 23 ) continue;
	const float gen_eta = cms3.mus_mc_p4().at(j).eta();
	const float gen_phi = cms3.mus_mc_p4().at(j).phi();
	matched = DeltaR(cand_eta,gen_eta,cand_phi,gen_phi) < 0.3;
      }
      if (!matched) continue;

      // See IsolationTools.cc for definitions of all these functions

      // at this point we have a Reco muon with pT > 5 and |eta| < 2.4
      const float neut_iso = mus_miniIso_nh().at(i) + mus_miniIso_em().at(i);
      //      const float absiso_uncorr = mus_miniIso_ch().at(i) + neut_iso; //mus_miniIso_nh().at(i) + mus_miniIso_em().at(i);
      const float neut_iso_st = mus_isoR03_pf_NeutralHadronEt().at(i) + mus_isoR03_pf_PhotonEt().at(i);
      //      const float absiso_uncorr_st = mus_isoR03_pf_ChargedHadronPt().at(i) + neut_iso_st; //mus_isoR03_pf_NeutralHadronEt().at(i) + mus_isoR03_pf_PhotonEt().at(i);

      const float dr = getMiniDR(cand_pt);
      const float effective_area = muEA03(i,1) * (dr/0.3) * (dr/0.3);
      const float correction_EA = rho_ctr * effective_area;
      const float effective_area_st = muEA03(i,1);
      const float correction_EA_st = rho_all * effective_area_st;

      h_all_EAcorrs.Fill(correction_EA, w_);
      h_all_EAcorrs_rel.Fill(correction_EA/cand_pt, w_);
      h_all_EAs.Fill(effective_area,w_);
      h_all_rhos_ctr.Fill(rho_ctr,w_);
      h_all_EAcorrs_st.Fill(correction_EA_st, w_);
      h_all_EAs_st.Fill(effective_area_st,w_);
      h_all_EAcorrs_rel_st.Fill(correction_EA_st/cand_pt, w_);
      h_all_rhos_all.Fill(rho_all,w_);

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
      // Overcorrection
      if ( neut_iso_st < correction_EA_st ) nOverEA[PU_index]++;
      if ( neut_iso_st < correction_dB_st ) nOverdB[PU_index]++;

      // TProfiles
      p_pt_EA.Fill(cand_pt,correction_EA_st,w_);
      p_pt_dB.Fill(cand_pt,correction_dB_st,w_);
      p_npv_EA.Fill(nPUvertices,correction_EA_st,w_);
      p_npv_dB.Fill(nPUvertices,correction_dB_st,w_);
      p_npv_rho_all.Fill(nPUvertices,rho_all,w_);
      p_npv_rho_ctr.Fill(nPUvertices,rho_ctr,w_);

      const float diffCorr = correction_dB - correction_EA;
      const float diffCorr_st = correction_dB_st - correction_EA_st;

      if (cand_pt < 50) {
	h_neut_soft.Fill(neut_iso, w_);
	h_neut_soft_st.Fill(neut_iso_st, w_);

	h_neut_soft_rel.Fill(neut_iso/cand_pt, w_);
	h_neut_soft_rel_st.Fill(neut_iso_st/cand_pt, w_);
      } else if (cand_pt < 200) {
	h_neut_med.Fill(neut_iso, w_);
	h_neut_med_st.Fill(neut_iso_st, w_);

	h_neut_med_rel.Fill(neut_iso/cand_pt, w_);
	h_neut_med_rel_st.Fill(neut_iso_st/cand_pt, w_);
      } else {
	h_neut_hard.Fill(neut_iso, w_);
	h_neut_hard_st.Fill(neut_iso_st, w_);

	h_neut_hard_rel.Fill(neut_iso/cand_pt, w_);
	h_neut_hard_rel_st.Fill(neut_iso_st/cand_pt, w_);
      }

      if (nPUvertices < 20){
	if (fabs(cand_eta) < 1.4) {
	  h_neut_low_bar.Fill(neut_iso, w_);
	  h_neut_low_bar_st.Fill(neut_iso_st, w_);

	  h_neut_low_bar_rel.Fill(neut_iso/cand_pt, w_);
	  h_neut_low_bar_rel_st.Fill(neut_iso_st/cand_pt, w_);

	  h_diffCorr_low_bar.Fill(diffCorr, w_);
	  h_diffCorr_low_bar_st.Fill(diffCorr_st, w_);

	  h_diffCorr_low_bar_rel.Fill(diffCorr/cand_pt, w_);
	  h_diffCorr_low_bar_rel_st.Fill(diffCorr_st/cand_pt, w_);
	} else {
	  h_neut_low_cap.Fill(neut_iso, w_);
	  h_neut_low_cap_st.Fill(neut_iso_st, w_);

	  h_neut_low_cap_rel.Fill(neut_iso/cand_pt, w_);
	  h_neut_low_cap_rel_st.Fill(neut_iso_st/cand_pt, w_);

	  h_diffCorr_low_cap.Fill(diffCorr, w_);
	  h_diffCorr_low_cap_st.Fill(diffCorr_st, w_);

	  h_diffCorr_low_cap_rel.Fill(diffCorr/cand_pt, w_);
	  h_diffCorr_low_cap_rel_st.Fill(diffCorr_st/cand_pt, w_);
	}
      } else if (nPUvertices < 40) {
	if (fabs(cand_eta) < 1.4) {
	  h_neut_med_bar.Fill(neut_iso, w_);
	  h_neut_med_bar_st.Fill(neut_iso, w_);

	  h_neut_med_bar_rel.Fill(neut_iso/cand_pt, w_);
	  h_neut_med_bar_rel_st.Fill(neut_iso/cand_pt, w_);

	  h_diffCorr_med_bar.Fill(diffCorr, w_);
	  h_diffCorr_med_bar_st.Fill(diffCorr_st, w_);

	  h_diffCorr_med_bar_rel.Fill(diffCorr/cand_pt, w_);
	  h_diffCorr_med_bar_rel_st.Fill(diffCorr_st/cand_pt, w_);
	} else {
	  h_neut_med_cap.Fill(neut_iso, w_);
	  h_neut_med_cap_st.Fill(neut_iso, w_);

	  h_neut_med_cap_rel.Fill(neut_iso/cand_pt, w_);
	  h_neut_med_cap_rel_st.Fill(neut_iso/cand_pt, w_);

	  h_diffCorr_med_cap.Fill(diffCorr, w_);
	  h_diffCorr_med_cap_st.Fill(diffCorr_st, w_);

	  h_diffCorr_med_cap_rel.Fill(diffCorr/cand_pt, w_);
	  h_diffCorr_med_cap_rel_st.Fill(diffCorr_st/cand_pt, w_);
	}
      }	else {
	if (fabs(cand_eta) < 1.4) {
	  h_neut_high_bar.Fill(neut_iso, w_);
	  h_neut_high_bar_st.Fill(neut_iso, w_);

	  h_neut_high_bar_rel.Fill(neut_iso/cand_pt, w_);
	  h_neut_high_bar_rel_st.Fill(neut_iso/cand_pt, w_);

	  h_diffCorr_high_bar.Fill(diffCorr, w_);
	  h_diffCorr_high_bar_st.Fill(diffCorr_st, w_);

	  h_diffCorr_high_bar_rel.Fill(diffCorr/cand_pt, w_);
	  h_diffCorr_high_bar_rel_st.Fill(diffCorr_st/cand_pt, w_);
	} else {
	  h_neut_high_cap.Fill(neut_iso, w_);
	  h_neut_high_cap_st.Fill(neut_iso, w_);

	  h_neut_high_cap_rel.Fill(neut_iso/cand_pt, w_);
	  h_neut_high_cap_rel_st.Fill(neut_iso/cand_pt, w_);

	  h_diffCorr_high_cap.Fill(diffCorr, w_);
	  h_diffCorr_high_cap_st.Fill(diffCorr_st, w_);

	  h_diffCorr_high_cap_rel.Fill(diffCorr/cand_pt, w_);
	  h_diffCorr_high_cap_rel_st.Fill(diffCorr_st/cand_pt, w_);
	}
      } //Histogram filling block

      
    } // Muon loop

  } // Event Loop

  cout << "Finished event loop." << endl;

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

    const float overEA = (float) nOverEA[i] / nEntriesByNpv[i];
    const float overdB = (float) nOverdB[i] / nEntriesByNpv[i];
    h_npv_over_EA.SetBinContent(i, overEA);
    h_npv_over_dB.SetBinContent(i, overdB);

    if (nEntriesByNpv[i] > 0) {
      h_median_EA.SetBinError(i,EAmedian/sqrt(nEntriesByNpv[i]));
      h_median_EA_st.SetBinError(i,EAmedian_st/sqrt(nEntriesByNpv[i]));
      h_median_dB.SetBinError(i,dBmedian/sqrt(nEntriesByNpv[i]));
      h_median_dB_st.SetBinError(i,dBmedian_st/sqrt(nEntriesByNpv[i]));
      h_median_EA_rel.SetBinError(i,EAmedian_rel/sqrt(nEntriesByNpv[i]));
      h_median_EA_rel_st.SetBinError(i,EAmedian_rel_st/sqrt(nEntriesByNpv[i]));
      h_median_dB_rel.SetBinError(i,dBmedian_rel/sqrt(nEntriesByNpv[i]));
      h_median_dB_rel_st.SetBinError(i,dBmedian_rel_st/sqrt(nEntriesByNpv[i]));

      h_npv_over_EA.SetBinError(i,overEA/sqrt(nEntriesByNpv[i]));
      h_npv_over_dB.SetBinError(i,overdB/sqrt(nEntriesByNpv[i]));
    }
  }

  cout << "Finished median calculations. About to save histograms." << endl;


  TFile * OutFile_ = new TFile(Form("%s.root", outname), "RECREATE");
  OutFile_->cd();

  // TProfiles
  p_pt_EA.Write();
  p_pt_dB.Write();
  p_npv_EA.Write();
  p_npv_dB.Write();
  p_npv_rho_all.Write();
  p_npv_rho_ctr.Write();

  h_npv_over_EA.Write();
  h_npv_over_dB.Write();

  // Mini

  h_all_EAs.Write();
  h_all_rhos_ctr.Write();

  h_median_EA.Write();
  h_median_dB.Write();
  h_all_EAcorrs.Write();
  h_all_dBcorrs.Write();
  h_median_EA_rel.Write();
  h_median_dB_rel.Write();
  h_all_EAcorrs_rel.Write();
  h_all_dBcorrs_rel.Write();

  h_diffCorr_low_bar.Write();
  h_diffCorr_low_cap.Write();
  h_diffCorr_med_bar.Write();
  h_diffCorr_med_cap.Write();
  h_diffCorr_high_bar.Write();
  h_diffCorr_high_cap.Write();

  h_diffCorr_low_bar_rel.Write();
  h_diffCorr_low_cap_rel.Write();
  h_diffCorr_med_bar_rel.Write();
  h_diffCorr_med_cap_rel.Write();
  h_diffCorr_high_bar_rel.Write();
  h_diffCorr_high_cap_rel.Write();

  h_neut_low_bar.Write();
  h_neut_low_cap.Write();
  h_neut_med_bar.Write();
  h_neut_med_cap.Write();
  h_neut_high_bar.Write();
  h_neut_high_cap.Write();

  h_neut_low_bar_rel.Write();
  h_neut_low_cap_rel.Write();
  h_neut_med_bar_rel.Write();
  h_neut_med_cap_rel.Write();
  h_neut_high_bar_rel.Write();
  h_neut_high_cap_rel.Write();

  h_neut_soft.Write();
  h_neut_med.Write();
  h_neut_hard.Write();

  h_neut_soft_rel.Write();
  h_neut_med_rel.Write();
  h_neut_hard_rel.Write();

  // DeltaR = 0.3

  h_all_EAs_st.Write();
  h_all_rhos_all.Write();

  h_median_EA_st.Write();
  h_median_dB_st.Write();
  h_all_EAcorrs_st.Write();
  h_all_dBcorrs_st.Write();
  h_median_EA_rel_st.Write();
  h_median_dB_rel_st.Write();
  h_all_EAcorrs_rel_st.Write();
  h_all_dBcorrs_rel_st.Write();

  h_diffCorr_low_bar_st.Write();
  h_diffCorr_low_cap_st.Write();
  h_diffCorr_med_bar_st.Write();
  h_diffCorr_med_cap_st.Write();
  h_diffCorr_high_bar_st.Write();
  h_diffCorr_high_cap_st.Write();

  h_diffCorr_low_bar_rel_st.Write();
  h_diffCorr_low_cap_rel_st.Write();
  h_diffCorr_med_bar_rel_st.Write();
  h_diffCorr_med_cap_rel_st.Write();
  h_diffCorr_high_bar_rel_st.Write();
  h_diffCorr_high_cap_rel_st.Write();

  h_neut_low_bar_st.Write();
  h_neut_low_cap_st.Write();
  h_neut_med_bar_st.Write();
  h_neut_med_cap_st.Write();
  h_neut_high_bar_st.Write();
  h_neut_high_cap_st.Write();

  h_neut_low_bar_rel_st.Write();
  h_neut_low_cap_rel_st.Write();
  h_neut_med_bar_rel_st.Write();
  h_neut_med_cap_rel_st.Write();
  h_neut_high_bar_rel_st.Write();
  h_neut_high_cap_rel_st.Write();

  h_neut_soft_st.Write();
  h_neut_med_st.Write();
  h_neut_hard_st.Write();

  h_neut_soft_rel_st.Write();
  h_neut_med_rel_st.Write();
  h_neut_hard_rel_st.Write();

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
