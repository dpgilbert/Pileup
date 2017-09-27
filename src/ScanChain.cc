// A cms3 looper to compare delta-beta and effective area pileup corrections

#include <iostream>
#include <assert.h>
#include "ScanChain.h"

using namespace std;
using namespace tas;

void PileupCorrection::ScanChain (TTree * tree, const char* outname, int MaxEvents) {
  const int nCorrsBins = 40;
  const float CorrsMin = 0;
  const float CorrsMax = 2;

  const int nPUbins = 12;
  const float PUmin = 0;
  const float PUmax = 60;

  TH1::SetDefaultSumw2(true);

  int nEntriesByNpv[13] = {0};
  vector<float> temp;
  temp.reserve(1000);

  TH1F h_all_rhos_ctr("h_all_rhos_ctr","Rho-CN",30,0,30);
  TH1F h_all_rhos_all("h_all_rhos_all","Rho-All",30,0,30);

  TProfile p_npv_rho_ctr("p_npv_rho_ctr","Rho-CN",nPUbins,PUmin,PUmax);
  TProfile p_npv_rho_all("p_npv_rho_all","Rho-All",nPUbins,PUmin,PUmax);

  TProfile p_npv_eff_ctr("p_npv_eff_ctr","MiniRelIso < 0.2, Rho-CN",nPUbins,PUmin,PUmax);
  TProfile p_npv_eff_all("p_npv_eff_all","MiniRelIso < 0.2, Rho-All",nPUbins,PUmin,PUmax);

  //////////
  // Mini //
  //////////

  TH1F h_all_bug_mini("h_all_bug_mini","All Bugged Corrections",nCorrsBins,CorrsMin,CorrsMax);
  TH1F h_all_nobug_mini("h_all_nobug_mini","All Un-Bugged Corrections",nCorrsBins,CorrsMin,CorrsMax);

  TH1F h_all_EAs_mini("h_all_EAs_mini","All Effective Areas",20,0,1);

  vector<vector<float> > bugged_corrs_byNpv_mini(13,temp);
  vector<vector<float> > unbugged_corrs_byNpv_mini(13,temp);
  TH1F h_median_bug_mini("h_median_bug_mini","Median Bugged Correction by N_{PV}",nPUbins,PUmin,PUmax);
  TH1F h_median_unbug_mini("h_median_unbug_mini","Median Un-bugged Correction by N_{PV}",nPUbins,PUmin,PUmax);

  TProfile p_pt_bug_mini("p_pt_bug_mini","Bugged Corrections by p_{T}",20,0,100);
  TProfile p_pt_nobug_mini("p_pt_nobug_mini","Un-bugged Corrections by p_{T}",20,0,100);

  TProfile p_npv_incl_bug_mini("p_npv_incl_bug_mini","Bugged EA Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_incl_nobug_mini("p_npv_incl_nobug_mini","Un-bugged Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_incl_neut_mini("p_npv_incl_neut_mini","Uncorrected Neutral Isolation by N_{PV}",nPUbins,PUmin,PUmax);

  TProfile p_npv_0p8_bug_mini("p_npv_0p8_bug_mini","Bugged Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_1p3_bug_mini("p_npv_1p3_bug_mini","Bugged Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p0_bug_mini("p_npv_2p0_bug_mini","Bugged Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p2_bug_mini("p_npv_2p2_bug_mini","Bugged Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p5_bug_mini("p_npv_2p5_bug_mini","Bugged Corrections by N_{PV}",nPUbins,PUmin,PUmax);

  TProfile p_npv_0p8_nobug_mini("p_npv_0p8_nobug_mini","Un-bugged Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_1p3_nobug_mini("p_npv_1p3_nobug_mini","Un-bugged Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p0_nobug_mini("p_npv_2p0_nobug_mini","Un-bugged Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p2_nobug_mini("p_npv_2p2_nobug_mini","Un-bugged Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p5_nobug_mini("p_npv_2p5_nobug_mini","Un-bugged Corrections by N_{PV}",nPUbins,PUmin,PUmax);

  TProfile p_npv_0p8_neut_mini("p_npv_0p8_neut_mini","Uncorrected Neutral Isolation by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_1p3_neut_mini("p_npv_1p3_neut_mini","Uncorrected Neutral Isolation by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p0_neut_mini("p_npv_2p0_neut_mini","Uncorrected Neutral Isolation by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p2_neut_mini("p_npv_2p2_neut_mini","Uncorrected Neutral Isolation by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p5_neut_mini("p_npv_2p5_neut_mini","Uncorrected Neutral Isolation by N_{PV}",nPUbins,PUmin,PUmax);
  
  TProfile p_npv_incl_corr_bug_mini("p_npv_incl_corr_bug_mini","Average Corrected Neut MiniIso by N_{PV} (bugged)",nPUbins,PUmin,PUmax);
  TProfile p_npv_incl_corr_nobug_mini("p_npv_incl_corr_nobug_mini","Average Corrected Neut MiniIso by N_{PV} (un-bugged)",nPUbins,PUmin,PUmax);

  ////////////////
  // Fixed Cone //
  ////////////////

  TH1F h_all_bug_fc("h_all_bug_fc","All Bugged Corrections",nCorrsBins,CorrsMin,CorrsMax);
  TH1F h_all_nobug_fc("h_all_nobug_fc","All Un-Bugged Corrections",nCorrsBins,CorrsMin,CorrsMax);

  TH1F h_all_EAs_fc("h_all_EAs_fc","All Effective Areas",20,0,1);
  TH1F h_all_rhos_ctr_fc("h_all_rhos_ctr_fc","All Rhos",30,0,30);
  TH1F h_all_rhos_all_fc("h_all_rhos_all_fc","All Rhos",30,0,30);

  vector<vector<float> > bugged_corrs_byNpv_fc(13,temp);
  vector<vector<float> > unbugged_corrs_byNpv_fc(13,temp);
  TH1F h_median_bug_fc("h_median_bug_fc","Median Bugged Correction by N_{PV}",nPUbins,PUmin,PUmax);
  TH1F h_median_unbug_fc("h_median_unbug_fc","Median Un-bugged Correction by N_{PV}",nPUbins,PUmin,PUmax);

  TProfile p_pt_bug_fc("p_pt_bug_fc","Bugged Corrections by p_{T}",20,0,100);
  TProfile p_pt_nobug_fc("p_pt_nobug_fc","Un-bugged Corrections by p_{T}",20,0,100);

  TProfile p_npv_incl_bug_fc("p_npv_incl_bug_fc","Bugged EA Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_incl_nobug_fc("p_npv_incl_nobug_fc","Un-bugged Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_incl_neut_fc("p_npv_incl_neut_fc","Uncorrected Neutral Isolation by N_{PV}",nPUbins,PUmin,PUmax);

  TProfile p_npv_0p8_bug_fc("p_npv_0p8_bug_fc","Bugged Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_1p3_bug_fc("p_npv_1p3_bug_fc","Bugged Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p0_bug_fc("p_npv_2p0_bug_fc","Bugged Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p2_bug_fc("p_npv_2p2_bug_fc","Bugged Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p5_bug_fc("p_npv_2p5_bug_fc","Bugged Corrections by N_{PV}",nPUbins,PUmin,PUmax);

  TProfile p_npv_0p8_nobug_fc("p_npv_0p8_nobug_fc","Un-bugged Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_1p3_nobug_fc("p_npv_1p3_nobug_fc","Un-bugged Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p0_nobug_fc("p_npv_2p0_nobug_fc","Un-bugged Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p2_nobug_fc("p_npv_2p2_nobug_fc","Un-bugged Corrections by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p5_nobug_fc("p_npv_2p5_nobug_fc","Un-bugged Corrections by N_{PV}",nPUbins,PUmin,PUmax);

  TProfile p_npv_0p8_neut_fc("p_npv_0p8_neut_fc","Uncorrected Neutral Isolation by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_1p3_neut_fc("p_npv_1p3_neut_fc","Uncorrected Neutral Isolation by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p0_neut_fc("p_npv_2p0_neut_fc","Uncorrected Neutral Isolation by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p2_neut_fc("p_npv_2p2_neut_fc","Uncorrected Neutral Isolation by N_{PV}",nPUbins,PUmin,PUmax);
  TProfile p_npv_2p5_neut_fc("p_npv_2p5_neut_fc","Uncorrected Neutral Isolation by N_{PV}",nPUbins,PUmin,PUmax);
  
  TProfile p_npv_incl_corr_bug_fc("p_npv_incl_corr_bug_fc","Average Corrected Neut Iso by N_{PV} (bugged)",nPUbins,PUmin,PUmax);
  TProfile p_npv_incl_corr_nobug_fc("p_npv_incl_corr_nobug_fc","Average Corrected Neut Iso by N_{PV} (un-bugged)",nPUbins,PUmin,PUmax);

  cms3.Init(tree);
  const int nEvents = tree->GetEntries();
  const int event_max = MaxEvents < 0 ? nEvents : min(MaxEvents,nEvents);
  for (int event = 0; event < event_max; event++) {
    tree->LoadTree(event);
    cms3.GetEntry(event);

    // Truth Npv
    //const int nPUvertices = cms3.puInfo_nPUvertices().at(12);    

    // Reco vertex counting copied from 
    // https://github.com/cmstas/MT2Analysis/blob/master/babymaker/ScanChain.cc#L509
    int nVert = 0;
    for(unsigned int ivtx=0; ivtx < cms3.evt_nvtxs(); ivtx++){
      if(cms3.vtxs_isFake().at(ivtx)) continue;
      if(cms3.vtxs_ndof().at(ivtx) <= 4) continue;
      if(fabs(cms3.vtxs_position().at(ivtx).z()) > 24) continue;
      if(cms3.vtxs_position().at(ivtx).Rho() > 2) continue;

      nVert++;
    }
    if (nVert == 0) continue;
    const int nPUvertices = nVert - 1;

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

      /*
      // Reject muons that overlap jets
      bool overlap = false;
      for (unsigned int j = 0; j < cms3.pfjets_p4().size() && !overlap; j++) {
	const float jet_eta = cms3.pfjets_p4().at(j).eta();
	const float jet_phi = cms3.pfjets_p4().at(j).phi();
	overlap = DeltaR(cand_eta,jet_eta,cand_phi,jet_phi) < 0.3;
      }
      if (overlap) continue;
      */

      // See IsolationTools.cc for definitions of all these functions

      // at this point we have a Reco muon with pT > 5 and |eta| < 2.4
      const float neut_iso_mini = mus_miniIso_nh().at(i) + mus_miniIso_em().at(i);
      const float neut_iso_fc = mus_isoR03_pf_NeutralHadronEt().at(i) + mus_isoR03_pf_PhotonEt().at(i);

      const float dr = getMiniDR(cand_pt);
      const float effective_area_mini = muEA03(i,1) * (dr/0.3) * (dr/0.3);
      const float effective_area_fc = muEA03(i,1);
      const float correction_bug_mini = rho_ctr * effective_area_mini;
      const float correction_unbug_mini = rho_all * effective_area_mini;
      const float correction_bug_fc = rho_ctr * effective_area_fc;
      const float correction_unbug_fc = rho_all * effective_area_fc;

      h_all_EAs_mini.Fill(effective_area_mini,w_);
      h_all_EAs_fc.Fill(effective_area_fc,w_);
      h_all_bug_mini.Fill(correction_bug_mini, w_);
      h_all_nobug_mini.Fill(correction_unbug_mini, w_);
      h_all_bug_fc.Fill(correction_bug_fc, w_);
      h_all_nobug_fc.Fill(correction_unbug_fc, w_);

      p_pt_bug_mini.Fill(cand_pt,correction_bug_mini,w_);
      p_pt_nobug_mini.Fill(cand_pt,correction_unbug_mini,w_);
      p_pt_bug_fc.Fill(cand_pt,correction_bug_fc,w_);
      p_pt_nobug_fc.Fill(cand_pt,correction_unbug_fc,w_);

      int PU_index = nPUvertices / 5;
      if (PU_index > 12) PU_index = 12;
      nEntriesByNpv[PU_index]++;
      bugged_corrs_byNpv_mini.at(PU_index).push_back(correction_bug_mini);
      unbugged_corrs_byNpv_mini.at(PU_index).push_back(correction_unbug_mini);
      bugged_corrs_byNpv_fc.at(PU_index).push_back(correction_bug_fc);
      unbugged_corrs_byNpv_fc.at(PU_index).push_back(correction_unbug_fc);

      p_npv_rho_ctr.Fill(nPUvertices, rho_ctr, w_);
      p_npv_rho_all.Fill(nPUvertices, rho_all, w_);

      const float ch_iso_mini = mus_miniIso_ch().at(i);

      if ( (ch_iso_mini + neut_iso_mini - correction_bug_mini) / cand_pt < 0.2 ) p_npv_eff_ctr.Fill(nPUvertices, 1, w_);
      else p_npv_eff_ctr.Fill(nPUvertices, 0, w_);

      if ( (ch_iso_mini + neut_iso_mini - correction_unbug_mini) / cand_pt < 0.2 ) p_npv_eff_all.Fill(nPUvertices, 1, w_);
      else p_npv_eff_all.Fill(nPUvertices, 0, w_);

      p_npv_incl_bug_mini.Fill(nPUvertices, correction_bug_mini, w_);
      p_npv_incl_nobug_mini.Fill(nPUvertices, correction_unbug_mini, w_);
      p_npv_incl_neut_mini.Fill(nPUvertices, neut_iso_mini, w_);
      p_npv_incl_bug_fc.Fill(nPUvertices, correction_bug_fc, w_);
      p_npv_incl_nobug_fc.Fill(nPUvertices, correction_unbug_fc, w_);
      p_npv_incl_neut_fc.Fill(nPUvertices, neut_iso_fc, w_);

      p_npv_incl_corr_bug_mini.Fill(nPUvertices, neut_iso_mini - correction_bug_mini, w_);
      p_npv_incl_corr_nobug_mini.Fill(nPUvertices, neut_iso_mini - correction_unbug_mini, w_);
      p_npv_incl_corr_bug_fc.Fill(nPUvertices, neut_iso_fc - correction_bug_fc, w_);
      p_npv_incl_corr_nobug_fc.Fill(nPUvertices, neut_iso_fc - correction_unbug_fc, w_);

      if (cand_eta < 0.8) {
	p_npv_0p8_bug_mini.Fill(nPUvertices, correction_bug_mini, w_);
	p_npv_0p8_nobug_mini.Fill(nPUvertices, correction_unbug_mini, w_);
	p_npv_0p8_neut_mini.Fill(nPUvertices, neut_iso_mini, w_);

	p_npv_0p8_bug_fc.Fill(nPUvertices, correction_bug_fc, w_);
	p_npv_0p8_nobug_fc.Fill(nPUvertices, correction_unbug_fc, w_);
	p_npv_0p8_neut_fc.Fill(nPUvertices, neut_iso_fc, w_);
      } else if (cand_eta < 1.3) {
	p_npv_1p3_bug_mini.Fill(nPUvertices, correction_bug_mini, w_);
	p_npv_1p3_nobug_mini.Fill(nPUvertices, correction_unbug_mini, w_);
	p_npv_1p3_neut_mini.Fill(nPUvertices, neut_iso_mini, w_);

	p_npv_1p3_bug_fc.Fill(nPUvertices, correction_bug_fc, w_);
	p_npv_1p3_nobug_fc.Fill(nPUvertices, correction_unbug_fc, w_);
	p_npv_1p3_neut_fc.Fill(nPUvertices, neut_iso_fc, w_);
      } else if (cand_eta < 2.0) {
	p_npv_2p0_bug_mini.Fill(nPUvertices, correction_bug_mini, w_);
	p_npv_2p0_nobug_mini.Fill(nPUvertices, correction_unbug_mini, w_);
	p_npv_2p0_neut_mini.Fill(nPUvertices, neut_iso_mini, w_);

	p_npv_2p0_bug_fc.Fill(nPUvertices, correction_bug_fc, w_);
	p_npv_2p0_nobug_fc.Fill(nPUvertices, correction_unbug_fc, w_);
	p_npv_2p0_neut_fc.Fill(nPUvertices, neut_iso_fc, w_);
      } else if (cand_eta < 2.2) {
	p_npv_2p2_bug_mini.Fill(nPUvertices, correction_bug_mini, w_);
	p_npv_2p2_nobug_mini.Fill(nPUvertices, correction_unbug_mini, w_);
	p_npv_2p2_neut_mini.Fill(nPUvertices, neut_iso_mini, w_);

	p_npv_2p2_bug_fc.Fill(nPUvertices, correction_bug_fc, w_);
	p_npv_2p2_nobug_fc.Fill(nPUvertices, correction_unbug_fc, w_);
	p_npv_2p2_neut_fc.Fill(nPUvertices, neut_iso_fc, w_);
      } else if (cand_eta < 2.5) {
	p_npv_2p5_bug_mini.Fill(nPUvertices, correction_bug_mini, w_);
	p_npv_2p5_nobug_mini.Fill(nPUvertices, correction_unbug_mini, w_);
	p_npv_2p5_neut_mini.Fill(nPUvertices, neut_iso_mini, w_);

	p_npv_2p5_bug_fc.Fill(nPUvertices, correction_bug_fc, w_);
	p_npv_2p5_nobug_fc.Fill(nPUvertices, correction_unbug_fc, w_);
	p_npv_2p5_neut_fc.Fill(nPUvertices, neut_iso_fc, w_);
      } // Eta-slicing histogram filling block      
    } // Muon loop    
  } // Event Loop

  cout << "Finished event loop." << endl;

  // get median corrections as function of Npv
  for (int i = 0; i < 13; i++) {
    sort(bugged_corrs_byNpv_mini.at(i).begin(),bugged_corrs_byNpv_mini.at(i).begin()+nEntriesByNpv[i]);
    sort(unbugged_corrs_byNpv_mini.at(i).begin(),unbugged_corrs_byNpv_mini.at(i).begin()+nEntriesByNpv[i]);
    sort(bugged_corrs_byNpv_fc.at(i).begin(),bugged_corrs_byNpv_fc.at(i).begin()+nEntriesByNpv[i]);
    sort(unbugged_corrs_byNpv_fc.at(i).begin(),unbugged_corrs_byNpv_fc.at(i).begin()+nEntriesByNpv[i]);
    float bugged_median_fc = 0, bugged_median_mini = 0;
    float unbugged_median_fc = 0, unbugged_median_mini = 0;
    if (nEntriesByNpv[i] > 0) {
      if (nEntriesByNpv[i] % 2) {
	// Odd
	bugged_median_mini = bugged_corrs_byNpv_mini.at(i).at( (nEntriesByNpv[i] / 2) );
	unbugged_median_mini = unbugged_corrs_byNpv_mini.at(i).at( (nEntriesByNpv[i] / 2) );
	bugged_median_fc = bugged_corrs_byNpv_fc.at(i).at( (nEntriesByNpv[i] / 2) );
	unbugged_median_fc = unbugged_corrs_byNpv_fc.at(i).at( (nEntriesByNpv[i] / 2) );
      }
      else {
	// Even
	bugged_median_mini = 0.5 * (bugged_corrs_byNpv_mini.at(i).at( (nEntriesByNpv[i] / 2) ) + bugged_corrs_byNpv_mini.at(i).at( (nEntriesByNpv[i] / 2) - 1));
	unbugged_median_mini = 0.5 * (unbugged_corrs_byNpv_mini.at(i).at( (nEntriesByNpv[i] / 2) ) + unbugged_corrs_byNpv_mini.at(i).at( (nEntriesByNpv[i] / 2) - 1));
	bugged_median_fc = 0.5 * (bugged_corrs_byNpv_fc.at(i).at( (nEntriesByNpv[i] / 2) ) + bugged_corrs_byNpv_fc.at(i).at( (nEntriesByNpv[i] / 2) - 1));
	unbugged_median_fc = 0.5 * (unbugged_corrs_byNpv_fc.at(i).at( (nEntriesByNpv[i] / 2) ) + unbugged_corrs_byNpv_fc.at(i).at( (nEntriesByNpv[i] / 2) - 1));

      }
    }
    h_median_bug_mini.SetBinContent(i,bugged_median_mini);
    h_median_unbug_mini.SetBinContent(i,unbugged_median_mini);
    h_median_bug_fc.SetBinContent(i,bugged_median_fc);
    h_median_unbug_fc.SetBinContent(i,unbugged_median_fc);

    if (nEntriesByNpv[i] > 0) {
      h_median_bug_mini.SetBinError(i,bugged_median_mini/sqrt(nEntriesByNpv[i]));
      h_median_unbug_mini.SetBinError(i,unbugged_median_fc/sqrt(nEntriesByNpv[i]));
      h_median_bug_fc.SetBinError(i,bugged_median_fc/sqrt(nEntriesByNpv[i]));
      h_median_unbug_fc.SetBinError(i,unbugged_median_fc/sqrt(nEntriesByNpv[i]));
    }
  }

  cout << "Finished median calculations. About to save histograms." << endl;


  TFile * OutFile_ = new TFile(Form("%s.root", outname), "RECREATE");
  OutFile_->cd();

  h_all_rhos_ctr.Write();
  h_all_rhos_all.Write();

  p_npv_rho_ctr.Write();
  p_npv_rho_all.Write();

  p_npv_eff_ctr.Write();
  p_npv_eff_all.Write();

  // Mini
  h_all_bug_mini.Write();
  h_all_nobug_mini.Write();

  h_all_EAs_mini.Write();

  h_median_bug_mini.Write();
  h_median_unbug_mini.Write();

  p_pt_bug_mini.Write();
  p_pt_nobug_mini.Write();

  p_npv_incl_bug_mini.Write();
  p_npv_incl_nobug_mini.Write();
  p_npv_incl_neut_mini.Write();

  p_npv_0p8_bug_mini.Write();
  p_npv_1p3_bug_mini.Write();
  p_npv_2p0_bug_mini.Write();
  p_npv_2p2_bug_mini.Write();
  p_npv_2p5_bug_mini.Write();

  p_npv_0p8_nobug_mini.Write();
  p_npv_1p3_nobug_mini.Write();
  p_npv_2p0_nobug_mini.Write();
  p_npv_2p2_nobug_mini.Write();
  p_npv_2p5_nobug_mini.Write();

  p_npv_0p8_neut_mini.Write();
  p_npv_1p3_neut_mini.Write();
  p_npv_2p0_neut_mini.Write();
  p_npv_2p2_neut_mini.Write();
  p_npv_2p5_neut_mini.Write();
  
  p_npv_incl_corr_bug_mini.Write();
  p_npv_incl_corr_nobug_mini.Write();

  // Fixed Cone
  h_all_bug_fc.Write();
  h_all_nobug_fc.Write();

  h_all_EAs_fc.Write();

  h_median_bug_fc.Write();
  h_median_unbug_fc.Write();

  p_pt_bug_fc.Write();
  p_pt_nobug_fc.Write();

  p_npv_incl_bug_fc.Write();
  p_npv_incl_nobug_fc.Write();
  p_npv_incl_neut_fc.Write();

  p_npv_0p8_bug_fc.Write();
  p_npv_1p3_bug_fc.Write();
  p_npv_2p0_bug_fc.Write();
  p_npv_2p2_bug_fc.Write();
  p_npv_2p5_bug_fc.Write();

  p_npv_0p8_nobug_fc.Write();
  p_npv_1p3_nobug_fc.Write();
  p_npv_2p0_nobug_fc.Write();
  p_npv_2p2_nobug_fc.Write();
  p_npv_2p5_nobug_fc.Write();

  p_npv_0p8_neut_fc.Write();
  p_npv_1p3_neut_fc.Write();
  p_npv_2p0_neut_fc.Write();
  p_npv_2p2_neut_fc.Write();
  p_npv_2p5_neut_fc.Write();
  
  p_npv_incl_corr_bug_fc.Write();
  p_npv_incl_corr_nobug_fc.Write();

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
