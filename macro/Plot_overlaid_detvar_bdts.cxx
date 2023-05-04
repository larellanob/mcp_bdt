void Plot_overlaid_detvar_bdts(TString tag = "def-th_1000kev_2hits")
{

  gStyle->SetOptStat(0);
  std::vector<int> mass_v {100, 150, 200, 300, 350, 400};
  //std::vector<int> mass_v {100};

  std::vector<TString> detvar_v
    {
      "wiremod_x",
      "wiremod_yz",
      "wiremod_anglexz",
      "wiremod_angleyz",
      "wiremod_dEdx",
      "ly_down",
      "ly_rayleigh",
      "ly_atten",
      "sce",
      "recomb"
    };

  for ( int mass: mass_v ) {
    bool draw_detvar_leg = true;
    // open cv file
    TCanvas *c1 = new TCanvas();
    TLegend *leg = new TLegend(0.6,0.65,0.9,0.9);
    TString cv_filename = Form("systematics/detvar_%s_%imev_v08_00_00_61/hist_BDT_scores_CV.root",
			       tag.Data(),mass);
    TFile *f_cv = new TFile(cv_filename);
    TH1F *cv_sig = (TH1F*)f_cv->Get("sig_hist");
    TH1F *cv_bkg = (TH1F*)f_cv->Get("bkg_hist");
    cv_sig->SetLineWidth(4);
    cv_bkg->SetLineWidth(4);
    cv_bkg->SetLineColor(kRed);
    cv_bkg->SetTitle(Form("BDT scores (CV and detvars) - %i MeV",mass));
    cv_bkg->Draw();
    cv_sig->Draw("same");
    leg->AddEntry(cv_sig,"CV signal");
    leg->AddEntry(cv_bkg,"CV background");
    // open detvar file and add to plot

    for ( auto detvar: detvar_v ) {
      TString detvar_filename = Form("systematics/detvar_%s_%imev_v08_00_00_61/hist_BDT_scores_%s.root",
				     tag.Data(),mass, detvar.Data());
      TFile *f_detvar = new TFile(detvar_filename);
      TH1F *detvar_sig = (TH1F*)f_detvar->Get("sig_hist");
      TH1F *detvar_bkg = (TH1F*)f_detvar->Get("bkg_hist");
      detvar_sig->SetLineWidth(1);
      detvar_bkg->SetLineWidth(1);
      detvar_bkg->SetLineColor(kRed+2);
      detvar_sig->SetLineColor(kBlue+2);
      detvar_bkg->Draw("same");
      detvar_sig->Draw("same");

      if ( draw_detvar_leg ) {
	leg->AddEntry(detvar_sig,"detvars signal");
	leg->AddEntry(detvar_bkg,"detvars background");
	draw_detvar_leg = false;
      }
      
    }
    leg->Draw("same");
    c1->SaveAs(Form("img/%s/systematics/overlaid_syst_%i.pdf",tag.Data(),mass));
    c1->SaveAs(Form("img/%s/systematics/overlaid_syst_%i.png",tag.Data(),mass));
  }
  
}
