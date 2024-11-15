void blip_plot_pair_variables_training(int mass = 100, bool b_validation = false)
{
  TString s_sample = "training";
  if ( b_validation ) {
    s_sample = "validation";
  }
  const char * ch_sample = s_sample.Data();
  TString filename;
  if ( mass != 0 ) {
    filename = Form("%s/blippairs_MillichargeBlipAna_500kev_2hits_%imev.root",ch_sample,mass);
  } else {
    filename = Form("%s/blippairs_MillichargeBlipAna_500kev_2hits_combinedmasses.root",ch_sample);
  }
  TFile *f = new TFile(filename.Data());
  TTree *tsig = (TTree *)f->Get("sig_blip_pairs");
  TTree *tmix = (TTree *)f->Get("mix_blip_pairs");
  TTree *tbkg = (TTree *)f->Get("bkg_blip_pairs");


  int n_branches = tsig->GetListOfBranches()->GetEntries();

  std::map<TString, TString> m_x_title;
  m_x_title["bp_px1"] = "Position of blip 1 in x (cm)";
  m_x_title["bp_py1"] = "Position of blip 1 in y (cm)";
  m_x_title["bp_pz1"] = "Position of blip 1 in z (cm)";
  m_x_title["bp_px2"] = "Position of blip 2 in x (cm)";
  m_x_title["bp_py2"] = "Position of blip 2 in y (cm)";
  m_x_title["bp_pz2"] = "Position of blip 2 in z (cm)";
  m_x_title["bp_dx1"] = "Extent of blip 1 in x (cm)";
  m_x_title["bp_dyz1"] = "Wire extent of blip 1 in yz (cm)";
  m_x_title["bp_dx2"] = "Extent of blip 2 in x (cm)";
  m_x_title["bp_dyz2"] = "Wire extent of blip 2 in yz (cm)";
  m_x_title["bp_emin"] = "Energy of least energetic blip in the pair (MeV)";
  m_x_title["bp_emax"] = "Energy of most energetic blip in the pair (MeV)";
  m_x_title["bp_sizemin"] = "Size of smallest blip in the pair (cm)";
  m_x_title["bp_sizemax"] = "Size of largest blip in the pair (cm)";
  m_x_title["bp_nhits_p2_min"] = "Min number of hits in collection plane";
  m_x_title["bp_nhits_p2_max"] = "Max number of hits in collection plane";
  m_x_title["bp_nhits_p01_min"] = "Min sum of hits in induction planes";
  m_x_title["bp_nhits_p01_max"] = "Max sum of hits in induction planes";
  m_x_title["bp_dist"] = "Distance among blips in pair (cm)";
  m_x_title["bp_th"] = "Angle #theta (deg)";
  m_x_title["bp_ph"] = "Angle #phi (deg)";
  m_x_title["bp_proxtrkdist_min"] = "Minimum distance from blip to nearest track (cm)";
  m_x_title["bp_proxtrkdist_max"] = "Maximum distance from blip to nearest track (cm)";
  m_x_title["bp_dist_closest_overall"]      = "#lbar#vec{#it{R}}#cbar of nearest third blip to blip-pair center (cm)";
  m_x_title["bp_perp_dist_closest_overall"] = "#it{R}_{#perp}  of nearest third blip to blip-pair center (cm)";
  m_x_title["bp_perp_dist_closest_behind"]  = "#it{R}_{#perp}  of nearest third blip behind of blip-pair center (cm)";
  m_x_title["bp_perp_dist_closest_between"] = "#it{R}_{#perp}  of nearest third blip between blip-pair (cm)";
  m_x_title["bp_perp_dist_closest_ahead"]   = "#it{R}_{#perp}  of nearest third blip ahead of blip-pair center (cm)";
  
  
  
  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas();
  gSystem->Exec(Form("mkdir -p img/%s/",ch_sample));
  c1->Print(Form("img/%s/blip_pairs_training_vars_%i.pdf[",ch_sample,mass));
  for ( int i = 0; i < n_branches; i++ ) {
    TString b_name = tsig->GetListOfBranches()->At(i)->GetName();

    tbkg->Draw(Form("%s>>hbkg",b_name.Data()));
    tmix->Draw(Form("%s>>hmix",b_name.Data()),"","same hist");
    tsig->Draw(Form("%s>>hsig",b_name.Data()),"","same hist");

    TH1F* h_sig = (TH1F*)gDirectory->Get("hsig");
    TH1F* h_bkg = (TH1F*)gDirectory->Get("hbkg");
    TH1F* h_mix = (TH1F*)gDirectory->Get("hmix");

    h_sig->SetLineColor(kBlue);
    h_bkg->SetLineColor(kRed);
    h_bkg->SetFillColorAlpha(kRed+1,0.5);
    h_mix->SetLineColor(kGreen);
    h_sig->SetLineWidth(2);
    h_bkg->SetLineWidth(2);
    h_mix->SetLineWidth(2);

    double norm_sig = h_bkg->Integral()/h_sig->Integral();
    double norm_mix = h_bkg->Integral()/h_mix->Integral();
    h_sig->Scale(norm_sig);
    h_mix->Scale(norm_mix);

    double max = std::max(h_bkg->GetMaximum(),h_sig->GetMaximum());
    max = std::max(max,h_mix->GetMaximum());
    h_bkg->SetMaximum(max*1.5);
    if ( m_x_title[b_name] == "" ) {
      h_bkg->SetTitle(Form(";%s",b_name.Data()));
    } else {
      h_bkg->SetTitle(Form(";%s",m_x_title[b_name].Data()));
    }

    TLegend myleg(0.4,0.7,0.9,0.9);
    myleg.AddEntry(h_bkg,Form("Background blip-pairs (%i)",
			      (int)h_bkg->GetEntries()));
    if ( mass != 0 ) {
      myleg.AddEntry(h_sig,
		     Form("Signal blip-pairs (%i MeV, %i, x%.1f)",
			  mass,
			  (int)h_sig->GetEntries(),
			  norm_sig));
    } else {
      myleg.AddEntry(h_sig,
		     Form("Signal blip-pairs (Combined masses, %i, x%.1f)",
			  (int)h_sig->GetEntries(),
			  norm_sig));
    }
    myleg.AddEntry(h_mix,
		   Form("Mixed blip-pairs (%i, x%.1f)",
			(int)h_mix->GetEntries(),
			norm_mix));
    myleg.Draw("same");
    
    c1->Print(Form("img/%s/blip_pairs_training_vars_%i.pdf",ch_sample,mass));
    delete h_sig;
    delete h_bkg;
    delete h_mix;
  }
  c1->Print(Form("img/%s/blip_pairs_training_vars_%i.pdf]",ch_sample,mass));
}
  
