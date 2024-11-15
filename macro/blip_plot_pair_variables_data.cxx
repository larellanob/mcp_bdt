void blip_plot_pair_variables_data(TString g_working_dir, int g_run)
{
  // working dir
  // /exp/uboone/data/users/arellano/bdt_mcp/blipreco/numi/v1/
  



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

  std::vector<TString> v_samples = {
    "dat",
    "ext",
    "fhc_nu_overlay",
    "sig_30mev"
  };

  g_working_dir = Form("%s/run%i",g_working_dir.Data(),g_run);

  TFile *f_data  = new TFile(Form("%s/blippairs/hist_blippairs_numi_run%i_data.root",
				g_working_dir.Data(),g_run));
  TFile *f_genie = new TFile(Form("%s/blippairs/hist_blippairs_numi_run%i_fhc_nu_overlay.root",
				  g_working_dir.Data(),g_run));
  TFile *f_ext   = new TFile(Form("%s/blippairs/hist_blippairs_numi_run%i_ext.root",
				g_working_dir.Data(),g_run));
  TFile *f_sig   = new TFile(Form("%s/blippairs/sig_validation/hist_blippairs_numi_run%i_sig_30mev.root",
				g_working_dir.Data(),g_run));
  TTree *tdat = (TTree *)f_data->Get("dat_blip_pairs");
  TTree *tnu = (TTree *)f_genie->Get("sig_blip_pairs");
  TTree *tcos = (TTree *)f_genie->Get("bkg_blip_pairs");
  TTree *text = (TTree *)f_ext->Get("dat_blip_pairs");
  TTree *tsig = (TTree *)f_sig->Get("sig_blip_pairs");

  double POT = 8.9e+18;
  double trig_on  = 232231.0;
  double trig_off = 489644.000000;

  double pot_genie = 0.0;
  TTree * t_pot_genie = (TTree*)f_genie->Get("total_pot");
  t_pot_genie->SetBranchAddress("tot_pot",&pot_genie);
  t_pot_genie->GetEntry(0);
  double pot_sig = 0.0;
  TTree * t_pot_sig = (TTree*)f_sig->Get("total_pot");
  t_pot_sig->SetBranchAddress("tot_pot",&pot_sig);
  t_pot_sig->GetEntry(0);
    


  TString img_dir = Form("%s/blippairs/img/",g_working_dir.Data());
  gSystem->Exec(Form("mkdir -p %s",img_dir.Data()));
  int n_branches = tsig->GetListOfBranches()->GetEntries();

  c1->Print(Form("%s/blip_pairs_vars_data.pdf[",img_dir.Data()));
  for ( int i = 0; i < n_branches; i++ ) {
    TString b_name = tsig->GetListOfBranches()->At(i)->GetName();

    tdat->Draw(Form("%s>>hbase",b_name.Data()),"","P ");
    tdat->Draw(Form("%s>>hdat",b_name.Data()),"","P same");
    tnu->Draw(Form("%s>>hnu",b_name.Data()),"","P same");
    tcos->Draw(Form("%s>>hcos",b_name.Data()),"","P same");
    text->Draw(Form("%s>>hext",b_name.Data()),"","P same");
    tsig->Draw(Form("%s>>hsig",b_name.Data()),"","P same");

    TH1F* h_base = (TH1F*)gDirectory->Get("hbase");
    TH1F* h_dat = (TH1F*)gDirectory->Get("hdat");
    TH1F* h_nu = (TH1F*)gDirectory->Get("hnu");
    TH1F* h_cos = (TH1F*)gDirectory->Get("hcos");
    TH1F* h_ext = (TH1F*)gDirectory->Get("hext");
    TH1F* h_sig = (TH1F*)gDirectory->Get("hsig");


    //h_bkg->SetLineColor(kRed);
    h_nu->SetFillColorAlpha(kRed-4,0.5);
    h_cos->SetFillColorAlpha(kCyan+2,0.5);
    h_ext->SetFillColorAlpha(kGreen-3,0.5);
    h_sig->SetLineColor(kBlue);
    h_sig->SetLineWidth(2);
    h_dat->SetMarkerStyle(kFullCircle);


    //h_sig->Scale(10000.*POT/pot_sig);
    double sig_norm = h_dat->Integral()/h_sig->Integral();
    h_sig->Scale(sig_norm);
    h_nu->Scale(POT/pot_genie);
    h_cos->Scale(POT/pot_genie);
    h_ext->Scale(trig_on/trig_off);

    THStack * h_stack = new THStack("h_stack","");
    h_stack->Add(h_nu);
    h_stack->Add(h_cos);
    h_stack->Add(h_ext);
    //h_base->Clear();
    h_base->Draw("goff");
    h_stack->Draw("same hist");
    h_sig->Draw("same hist");
    h_dat->Draw("E0 same");

    
    double max = std::max(h_stack->GetMaximum(),h_sig->GetMaximum());
    max = std::max(max,h_dat->GetMaximum());
    h_base->SetMaximum(max*1.5);

    /*
    double max = std::max(h_bkg->GetMaximum(),h_sig->GetMaximum());
    max = std::max(max,h_mix->GetMaximum());
    h_bkg->SetMaximum(max*1.5);
    */
    
    if ( m_x_title[b_name] == "" ) {
      h_stack->SetTitle(Form(";%s",b_name.Data()));
    } else {
      h_base->SetTitle(Form("Blip pairs;%s",m_x_title[b_name].Data()));
    }

    TLegend myleg(0.4,0.7,0.9,0.9);
    myleg.AddEntry(h_nu,
		   Form("Genie #nu (%i pairs)",
			(int)h_nu->GetEntries())
		   );
    myleg.AddEntry(h_cos,
		   Form("Cosmic (%i pairs)",
			(int)h_cos->GetEntries())
		   );
    myleg.AddEntry(h_ext,
		   Form("Beam-off (%i pairs)",
			(int)h_ext->GetEntries())
		   );
    myleg.AddEntry(h_sig,
		   Form("30 MeV Signal (%.2e POT %i pairs)",
			pot_sig*sig_norm,
			(int)h_sig->GetEntries())
		   );
    myleg.AddEntry(h_dat,
		   Form("Data (Run 1, %.2e POT, %i pairs)",
			POT,
			(int)h_dat->GetEntries())
		   );
    myleg.Draw("same");
    
    c1->Print(Form("%s/blip_pairs_vars_data.pdf",img_dir.Data()));
    delete h_dat;
    delete h_nu;
    delete h_cos;
    delete h_ext;
    delete h_sig;
    delete h_stack;
    delete h_base;
  }
  c1->Print(Form("%s/blip_pairs_vars_data.pdf]",img_dir.Data()));
}
  
