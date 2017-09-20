// A cms3 looper to compare delta-beta and effective area pileup corrections

#include <iostream>
#include "ScanChain.h"

using namespace std;
using namespace tas;

void PileupCorrection::ScanChain (TTree * tree, const char* outname, int MaxEvents) {
  const int nIsoDiffBins = 40;
  const float IsoDiffMin = -1;
  const float IsoDiffMax = 1;

  const int nCorrsBins = 24;
  const float CorrsMin = 0;
  const float CorrsMax = 5;
  const int nRelCorrsBins = 24;
  const float RelCorrsMin = 0;
  const float RelCorrsMax = 1;

  const int nPUbins = 12;
  const float PUmin = 0;
  const float PUmax = 60;

  TH1::SetDefaultSumw2(true);

  // TProfiles (standard only)
  TProfile p_pt_EA("p_pt_EA","EA Corrections by p_{T}",20,0,100);
  TProfile p_pt_dB("p_pt_dB","#Delta#beta Corrections by p_{T}",20,0,100);
  TProfile p_npv_rho_all("p_npv_rho_all","#rho_{all} by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_rho_ctr("p_npv_rho_ctr","#rho_{ctr} by N_{PV}",nPUbins,PUmin,PUmax);

  TProfile p_npv_incl_EA("p_npv_incl_EA","EA Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_incl_dB("p_npv_incl_dB","#Delta#beta Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_incl_neut("p_npv_incl_neut","Neut Isolation by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_incl_EA_st("p_npv_incl_EA_st","EA Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_incl_dB_st("p_npv_incl_dB_st","#Delta#beta Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_incl_neut_st("p_npv_incl_neut_st","Neut Isolation by N_{PV}",nPUbins,PUmin,PUmax);

  TProfile p_npv_0p8_EA("p_npv_0p8_EA","EA Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_1p3_EA("p_npv_1p3_EA","EA Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p0_EA("p_npv_2p0_EA","EA Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p2_EA("p_npv_2p2_EA","EA Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p5_EA("p_npv_2p5_EA","EA Corrections by N_{PV}",nPUbins,PUmin,PUmax);

  TProfile p_npv_0p8_dB("p_npv_0p8_dB","#Delta#beta Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_1p3_dB("p_npv_1p3_dB","#Delta#beta Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p0_dB("p_npv_2p0_dB","#Delta#beta Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p2_dB("p_npv_2p2_dB","#Delta#beta Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p5_dB("p_npv_2p5_dB","#Delta#beta Corrections by N_{PV}",nPUbins,PUmin,PUmax);

  TProfile p_npv_0p8_EA_st("p_npv_0p8_EA_st","EA Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_1p3_EA_st("p_npv_1p3_EA_st","EA Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p0_EA_st("p_npv_2p0_EA_st","EA Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p2_EA_st("p_npv_2p2_EA_st","EA Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p5_EA_st("p_npv_2p5_EA_st","EA Corrections by N_{PV}",nPUbins,PUmin,PUmax);

  TProfile p_npv_0p8_dB_st("p_npv_0p8_dB_st","#Delta#beta Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_1p3_dB_st("p_npv_1p3_dB_st","#Delta#beta Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p0_dB_st("p_npv_2p0_dB_st","#Delta#beta Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p2_dB_st("p_npv_2p2_dB_st","#Delta#beta Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p5_dB_st("p_npv_2p5_dB_st","#Delta#beta Corrections by N_{PV}",nPUbins,PUmin,PUmax);

  TProfile p_npv_0p8_neut("p_npv_0p8_neut","Neutral Isolation by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_1p3_neut("p_npv_1p3_neut","Neutral Isolation by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p0_neut("p_npv_2p0_neut","Neutral Isolation by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p2_neut("p_npv_2p2_neut","Neutral Isolation by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p5_neut("p_npv_2p5_neut","Neutral Isolation by N_{PV}",nPUbins,PUmin,PUmax);
  
  TProfile p_npv_0p8_neut_st("p_npv_0p8_neut_st","Neutral Isolation by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_1p3_neut_st("p_npv_1p3_neut_st","Neutral Isolation by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p0_neut_st("p_npv_2p0_neut_st","Neutral Isolation by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p2_neut_st("p_npv_2p2_neut_st","Neutral Isolation by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p5_neut_st("p_npv_2p5_neut_st","Neutral Isolation by N_{PV}",nPUbins,PUmin,PUmax);

  TH1F h_npv_over_EA("h_npv_over_EA","Fraction of EA Overcorrections by N_{PV}",nPUbins,PUmin,PUmax);
  TH1F h_npv_over_dB("h_npv_over_dB","Fraction of #Delta#beta Overcorrections by N_{PV}",nPUbins,PUmin,PUmax);
  // use nEntriesByNpv array defined below as denom
  int nOverEA[13] = {0};
  int nOverdB[13] = {0};

  //////////
  // Mini //
  //////////
  // Absolute Iso
  TH1F h_all_EAcorrs("h_all_EAcorrs","All Effective Area Corrections",nCorrsBins,CorrsMin,CorrsMax);
  TH1F h_all_dBcorrs("h_all_dBcorrs","All #Delta#beta Corrections",nCorrsBins,CorrsMin,CorrsMax);

  TH1F h_all_EAs("h_all_EAs","All Effective Areas",20,0,10);
  TH1F h_all_rhos_ctr("h_all_rhos_ctr","All Rhos",30,0,30);

  vector<float> temp;
  temp.reserve(1000);
  int nEntriesByNpv[13] = {0};
  vector<vector<float> > EA_corrs_byNpv(13,temp);
  vector<vector<float> > dB_corrs_byNpv(13,temp);
  TH1F h_median_EA("h_median_EA","Median Effective Area Correction by Number of Pileup Vertices",nPUbins,PUmin,PUmax);
  TH1F h_median_dB("h_median_dB","Median #Delta#beta Correction by Number of Pileup Vertices",nPUbins,PUmin,PUmax);

  TH1F h_diff_0p8_EA("h_diff_0p8_EA","Neutral MiniIso - Effective Area for |#eta| < 0.800", nIsoDiffBins,IsoDiffMin,IsoDiffMax);
  TH1F h_diff_1p3_EA("h_diff_1p3_EA","Neutral MiniIso - Effective Area for 0.800 < |#eta| < 1.300", nIsoDiffBins,IsoDiffMin,IsoDiffMax);
  TH1F h_diff_2p0_EA("h_diff_2p0_EA","Neutral MiniIso - Effective Area for 1.300 < |#eta| < 2.000", nIsoDiffBins,IsoDiffMin,IsoDiffMax);
  TH1F h_diff_2p2_EA("h_diff_2p2_EA","Neutral MiniIso - Effective Area for 2.000 < |#eta| < 2.200", nIsoDiffBins,IsoDiffMin,IsoDiffMax);
  TH1F h_diff_2p5_EA("h_diff_2p5_EA","Neutral MiniIso - Effective Area for 2.200 < |#eta| < 2.500", nIsoDiffBins,IsoDiffMin,IsoDiffMax);

  TH1F h_diff_0p8_dB("h_diff_0p8_dB","Neutral MiniIso - #Delta#beta for |#eta| < 0.800", nIsoDiffBins,IsoDiffMin,IsoDiffMax);
  TH1F h_diff_1p3_dB("h_diff_1p3_dB","Neutral MiniIso - #Delta#beta for 0.800 < |#eta| < 1.300", nIsoDiffBins,IsoDiffMin,IsoDiffMax);
  TH1F h_diff_2p0_dB("h_diff_2p0_dB","Neutral MiniIso - #Delta#beta for 1.300 < |#eta| < 2.000", nIsoDiffBins,IsoDiffMin,IsoDiffMax);
  TH1F h_diff_2p2_dB("h_diff_2p2_dB","Neutral MiniIso - #Delta#beta for 2.000 < |#eta| < 2.200", nIsoDiffBins,IsoDiffMin,IsoDiffMax);
  TH1F h_diff_2p5_dB("h_diff_2p5_dB","Neutral MiniIso - #Delta#beta for 2.200 < |#eta| < 2.500", nIsoDiffBins,IsoDiffMin,IsoDiffMax);

  // Relative

  TH1F h_all_EAcorrs_rel("h_all_EAcorrs_rel","All Effective Area Corrections",nRelCorrsBins,RelCorrsMin,RelCorrsMax);
  TH1F h_all_dBcorrs_rel("h_all_dBcorrs_rel","All #Delta#beta Corrections",nRelCorrsBins,RelCorrsMin,RelCorrsMax);

  vector<vector<float> > EA_corrs_byNpv_rel(13,temp);
  vector<vector<float> > dB_corrs_byNpv_rel(13,temp);
  TH1F h_median_EA_rel("h_median_EA_rel","Median Effective Area Correction by Number of Pileup Vertices",nPUbins,PUmin,PUmax);
  TH1F h_median_dB_rel("h_median_dB_rel","Median #Delta#beta Correction by Number of Pileup Vertices",nPUbins,PUmin,PUmax);

  //////////////////
  // Standard Iso //
  //////////////////
  TH1F h_all_EAcorrs_st("h_all_EAcorrs_st","All Effective Area Corrections",nCorrsBins,CorrsMin,CorrsMax);
  TH1F h_all_dBcorrs_st("h_all_dBcorrs_st","All #Delta#beta Corrections",nCorrsBins,CorrsMin,CorrsMax);

  TH1F h_all_EAs_st("h_all_EAs_st","All Effective Areas",20,0,10);
  TH1F h_all_rhos_all("h_all_rhos_all","All Rhos",30,0,30);

  vector<vector<float> > EA_corrs_byNpv_st(13,temp);
  vector<vector<float> > dB_corrs_byNpv_st(13,temp);
  TH1F h_median_EA_st("h_median_EA_st","Median Effective Area Correction by Number of Pileup Vertices",nPUbins,PUmin,PUmax);
  TH1F h_median_dB_st("h_median_dB_st","Median #Delta#beta Correction by Number of Pileup Vertices",nPUbins,PUmin,PUmax);

  TH1F h_diff_0p8_EA_st("h_diff_0p8_EA_st","Neutral #DeltaR=0.3 Iso - Effective Area for |#eta| < 0.800", nIsoDiffBins,IsoDiffMin,IsoDiffMax);
  TH1F h_diff_1p3_EA_st("h_diff_1p3_EA_st","Neutral #DeltaR=0.3 Iso - Effective Area for 0.800 < |#eta| < 1.300", nIsoDiffBins,IsoDiffMin,IsoDiffMax);
  TH1F h_diff_2p0_EA_st("h_diff_2p0_EA_st","Neutral #DeltaR=0.3 Iso - Effective Area for 1.300 < |#eta| < 2.000", nIsoDiffBins,IsoDiffMin,IsoDiffMax);
  TH1F h_diff_2p2_EA_st("h_diff_2p2_EA_st","Neutral #DeltaR=0.3 Iso - Effective Area for 2.000 < |#eta| < 2.200", nIsoDiffBins,IsoDiffMin,IsoDiffMax);
  TH1F h_diff_2p5_EA_st("h_diff_2p5_EA_st","Neutral #DeltaR=0.3 Iso - Effective Area for 2.200 < |#eta| < 2.500", nIsoDiffBins,IsoDiffMin,IsoDiffMax);

  TH1F h_diff_0p8_dB_st("h_diff_0p8_dB_st","Neutral #DeltaR=0.3 Iso - #Delta#beta for |#eta| < 0.800", nIsoDiffBins,IsoDiffMin,IsoDiffMax);
  TH1F h_diff_1p3_dB_st("h_diff_1p3_dB_st","Neutral #DeltaR=0.3 Iso - #Delta#beta for 0.800 < |#eta| < 1.300", nIsoDiffBins,IsoDiffMin,IsoDiffMax);
  TH1F h_diff_2p0_dB_st("h_diff_2p0_dB_st","Neutral #DeltaR=0.3 Iso - #Delta#beta for 1.300 < |#eta| < 2.000", nIsoDiffBins,IsoDiffMin,IsoDiffMax);
  TH1F h_diff_2p2_dB_st("h_diff_2p2_dB_st","Neutral #DeltaR=0.3 Iso - #Delta#beta for 2.000< |#eta| < 2.200", nIsoDiffBins,IsoDiffMin,IsoDiffMax);
  TH1F h_diff_2p5_dB_st("h_diff_2p5_dB_st","Neutral #DeltaR=0.3 Iso - #Delta#beta for 2.200 < |#eta| < 2.500", nIsoDiffBins,IsoDiffMin,IsoDiffMax);

  // Relative

  TH1F h_all_EAcorrs_rel_st("h_all_EAcorrs_rel_st","All Effective Area Corrections",nRelCorrsBins,RelCorrsMin,RelCorrsMax);
  TH1F h_all_dBcorrs_rel_st("h_all_dBcorrs_rel_st","All #Delta#beta Corrections",nRelCorrsBins,RelCorrsMin,RelCorrsMax);

  vector<vector<float> > EA_corrs_byNpv_rel_st(13,temp);
  vector<vector<float> > dB_corrs_byNpv_rel_st(13,temp);
  TH1F h_median_EA_rel_st("h_median_EA_rel_st","Median Effective Area Correction by Number of Pileup Vertices",nPUbins,PUmin,PUmax);
  TH1F h_median_dB_rel_st("h_median_dB_rel_st","Median #Delta#beta Correction by Number of Pileup Vertices",nPUbins,PUmin,PUmax);

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

    h_all_rhos_ctr.Fill(rho_ctr,w_);
    h_all_rhos_all.Fill(rho_all,w_);

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
      h_all_EAcorrs_st.Fill(correction_EA_st, w_);
      h_all_EAs_st.Fill(effective_area_st,w_);
      h_all_EAcorrs_rel_st.Fill(correction_EA_st/cand_pt, w_);

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
      p_npv_rho_all.Fill(nPUvertices,rho_all,w_);
      p_npv_rho_ctr.Fill(nPUvertices,rho_ctr,w_);

      p_npv_incl_EA.Fill(nPUvertices, correction_EA, w_);
      p_npv_incl_dB.Fill(nPUvertices, correction_dB, w_);
      p_npv_incl_EA_st.Fill(nPUvertices, correction_EA_st, w_);
      p_npv_incl_dB_st.Fill(nPUvertices, correction_dB_st, w_);
      p_npv_incl_neut.Fill(nPUvertices, neut_iso, w_);
      p_npv_incl_neut_st.Fill(nPUvertices, neut_iso_st, w_);

      if (cand_eta < 0.8) {
	h_diff_0p8_EA.Fill(neut_iso - correction_EA, w_);
	h_diff_0p8_dB.Fill(neut_iso - correction_dB, w_);
	h_diff_0p8_EA_st.Fill(neut_iso_st - correction_EA_st, w_);
	h_diff_0p8_dB_st.Fill(neut_iso_st - correction_dB_st, w_);
	p_npv_0p8_EA.Fill(nPUvertices, correction_EA, w_);
	p_npv_0p8_dB.Fill(nPUvertices, correction_dB, w_);
	p_npv_0p8_EA_st.Fill(nPUvertices, correction_EA_st, w_);
	p_npv_0p8_dB_st.Fill(nPUvertices, correction_dB_st, w_);
	p_npv_0p8_neut.Fill(nPUvertices, neut_iso, w_);
	p_npv_0p8_neut_st.Fill(nPUvertices, neut_iso_st, w_);
      } else if (cand_eta < 1.3) {
	h_diff_1p3_EA.Fill(neut_iso - correction_EA, w_);
	h_diff_1p3_dB.Fill(neut_iso - correction_dB, w_);
	h_diff_1p3_EA_st.Fill(neut_iso_st - correction_EA_st, w_);
	h_diff_1p3_dB_st.Fill(neut_iso_st - correction_dB_st, w_);
	p_npv_1p3_EA.Fill(nPUvertices, correction_EA, w_);
	p_npv_1p3_dB.Fill(nPUvertices, correction_dB, w_);
	p_npv_1p3_EA_st.Fill(nPUvertices, correction_EA_st, w_);
	p_npv_1p3_dB_st.Fill(nPUvertices, correction_dB_st, w_);
	p_npv_1p3_neut.Fill(nPUvertices, neut_iso, w_);
	p_npv_1p3_neut_st.Fill(nPUvertices, neut_iso_st, w_);
      } else if (cand_eta < 2.0) {
	h_diff_2p0_EA.Fill(neut_iso - correction_EA, w_);
	h_diff_2p0_dB.Fill(neut_iso - correction_dB, w_);
	h_diff_2p0_EA_st.Fill(neut_iso_st - correction_EA_st, w_);
	h_diff_2p0_dB_st.Fill(neut_iso_st - correction_dB_st, w_);
	p_npv_2p0_EA.Fill(nPUvertices, correction_EA, w_);
	p_npv_2p0_dB.Fill(nPUvertices, correction_dB, w_);
	p_npv_2p0_EA_st.Fill(nPUvertices, correction_EA_st, w_);
	p_npv_2p0_dB_st.Fill(nPUvertices, correction_dB_st, w_);
	p_npv_2p0_neut.Fill(nPUvertices, neut_iso, w_);
	p_npv_2p0_neut_st.Fill(nPUvertices, neut_iso_st, w_);
      } else if (cand_eta < 2.2) {
	h_diff_2p2_EA.Fill(neut_iso - correction_EA, w_);
	h_diff_2p2_dB.Fill(neut_iso - correction_dB, w_);
	h_diff_2p2_EA_st.Fill(neut_iso_st - correction_EA_st, w_);
	h_diff_2p2_dB_st.Fill(neut_iso_st - correction_dB_st, w_);
	p_npv_2p2_EA.Fill(nPUvertices, correction_EA, w_);
	p_npv_2p2_dB.Fill(nPUvertices, correction_dB, w_);
	p_npv_2p2_EA_st.Fill(nPUvertices, correction_EA_st, w_);
	p_npv_2p2_dB_st.Fill(nPUvertices, correction_dB_st, w_);
	p_npv_2p2_neut.Fill(nPUvertices, neut_iso, w_);
	p_npv_2p2_neut_st.Fill(nPUvertices, neut_iso_st, w_);
      } else if (cand_eta < 2.5) {
	h_diff_2p5_EA.Fill(neut_iso - correction_EA, w_);
	h_diff_2p5_dB.Fill(neut_iso - correction_dB, w_);
	h_diff_2p5_EA_st.Fill(neut_iso_st - correction_EA_st, w_);
	h_diff_2p5_dB_st.Fill(neut_iso_st - correction_dB_st, w_);
	p_npv_2p5_EA.Fill(nPUvertices, correction_EA, w_);
	p_npv_2p5_dB.Fill(nPUvertices, correction_dB, w_);
	p_npv_2p5_EA_st.Fill(nPUvertices, correction_EA_st, w_);
	p_npv_2p5_dB_st.Fill(nPUvertices, correction_dB_st, w_);
	p_npv_2p5_neut.Fill(nPUvertices, neut_iso, w_);
	p_npv_2p5_neut_st.Fill(nPUvertices, neut_iso_st, w_);
      } // Eta-slicing histogram filling block      
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
  p_npv_rho_all.Write();
  p_npv_rho_ctr.Write();

  p_npv_incl_EA.Write();
  p_npv_incl_dB.Write();
  p_npv_incl_neut.Write();
  p_npv_incl_EA_st.Write();
  p_npv_incl_dB_st.Write();
  p_npv_incl_neut_st.Write();

  p_npv_0p8_EA.Write();
  p_npv_1p3_EA.Write();
  p_npv_2p0_EA.Write();
  p_npv_2p2_EA.Write();
  p_npv_2p5_EA.Write();

  p_npv_0p8_dB.Write();
  p_npv_1p3_dB.Write();
  p_npv_2p0_dB.Write();
  p_npv_2p2_dB.Write();
  p_npv_2p5_dB.Write();

  p_npv_0p8_EA_st.Write();
  p_npv_1p3_EA_st.Write();
  p_npv_2p0_EA_st.Write();
  p_npv_2p2_EA_st.Write();
  p_npv_2p5_EA_st.Write();

  p_npv_0p8_dB_st.Write();
  p_npv_1p3_dB_st.Write();
  p_npv_2p0_dB_st.Write();
  p_npv_2p2_dB_st.Write();
  p_npv_2p5_dB_st.Write();

  p_npv_0p8_neut.Write();
  p_npv_1p3_neut.Write();
  p_npv_2p0_neut.Write();
  p_npv_2p2_neut.Write();
  p_npv_2p5_neut.Write();

  p_npv_0p8_neut_st.Write();
  p_npv_1p3_neut_st.Write();
  p_npv_2p0_neut_st.Write();
  p_npv_2p2_neut_st.Write();
  p_npv_2p5_neut_st.Write();  

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

  h_diff_0p8_EA.Write();
  h_diff_1p3_EA.Write();
  h_diff_2p0_EA.Write();
  h_diff_2p2_EA.Write();
  h_diff_2p5_EA.Write();

  h_diff_0p8_dB.Write();
  h_diff_1p3_dB.Write();
  h_diff_2p0_dB.Write();
  h_diff_2p2_dB.Write();
  h_diff_2p5_dB.Write();


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

  h_diff_0p8_EA_st.Write();
  h_diff_1p3_EA_st.Write();
  h_diff_2p0_EA_st.Write();
  h_diff_2p2_EA_st.Write();
  h_diff_2p5_EA_st.Write();

  h_diff_0p8_dB_st.Write();
  h_diff_1p3_dB_st.Write();
  h_diff_2p0_dB_st.Write();
  h_diff_2p2_dB_st.Write();
  h_diff_2p5_dB_st.Write();

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
