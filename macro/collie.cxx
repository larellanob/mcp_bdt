void cout_bins(TH1F* h )
{
  int nbins = h->GetNbinsX();
  std::cout << "Number of bins for " << h->GetName() << ": "  << nbins << std::endl;
  std::cout << "[ ";
  for ( int i = 1; i < nbins+1; i++ ) {
    std::cout << h->GetBinLowEdge(i) << ": " << h->GetBinContent(i);
    if ( i < nbins ) {
      std::cout << ", ";
    }
  }
  std::cout << " ]" << std::endl;
}

void cout_bins(TH1F* bdt, TH1F *sys)
{
  int nbins = bdt->GetNbinsX();
  std::cout << "Number of bins for " << bdt->GetName() << ": "  << nbins << std::endl;
  std::cout << "[ ";
  for ( int i = 1; i < nbins+1; i++ ) {
    std::cout << bdt->GetBinLowEdge(i) << ": " << bdt->GetBinContent(i)
	      << "±" << bdt->GetBinError(i)<<"±"<<bdt->GetBinContent(i)*sys->GetBinContent(i);
    if ( i < nbins ) {
      std::cout << ", ";
    }
  }
  std::cout << " ]" << std::endl;
}

void collie()
{
  TString analysis = "def-th_1000kev_2hits";
  std::vector<int> mass_v {100, 150, 200, 300, 350, 400};
  //std::vector<int> mass_v {100};
  for ( int mass: mass_v ) {
    std::cout << "--------------------" << std::endl;
    std::cout << "Mass: " << mass << std::endl;
    // import bdt output file
    TString bdt_filename = Form("root/%s/hist_BDT_scores_%imev.root",analysis.Data(),mass);
    TFile *f_bdt = new TFile(bdt_filename);
    // import systematics file
    TString sys_sig_filename = Form("systematics/detvar_%s_%imev_v08_00_00_61/Quadrature_signal.root",
				    analysis.Data(),mass);
    TString sys_bkg_filename = Form("systematics/detvar_%s_%imev_v08_00_00_61/Quadrature_background.root",
				    analysis.Data(),mass);
    TFile *f_sys_sig = new TFile(sys_sig_filename);
    TFile *f_sys_bkg = new TFile(sys_bkg_filename);
    // get bdt histograms for sig and bkg
    TH1F *bdt_sig = (TH1F*)f_bdt->Get("sig_hist");
    TH1F *bdt_bkg = (TH1F*)f_bdt->Get("bkg_hist");
    // bdt pot/event scaling
    TTree *pot_t = (TTree*)f_bdt->Get("total_pot");
    double pot, evt;
    pot_t->SetBranchAddress("tot_pot",&pot);
    pot_t->SetBranchAddress("tot_evt",&evt);
    pot_t->GetEntry(0);
    std::cout << "POT, events: " << pot << " " << evt << std::endl;
    float pot_runs123 = 1.5e21;
    float eve_runs123 = 3081362.;
    bdt_sig->Scale(pot_runs123/pot);
    bdt_bkg->Scale(eve_runs123/evt);
    
    TH1F *sys_sig = (TH1F*)f_sys_sig->Get("quadrature");
    TH1F *sys_bkg = (TH1F*)f_sys_bkg->Get("quadrature");
    TH1F *bdt_sig_full = (TH1F*)bdt_sig->Clone();
    TH1F *bdt_bkg_full = (TH1F*)bdt_bkg->Clone();
    for ( int i = 1; i < bdt_sig_full->GetNbinsX()+1; i++ ) {
      float sig_err = bdt_sig->GetBinError(i)+bdt_sig->GetBinContent(i)*sys_sig->GetBinContent(i);
      float bkg_err = bdt_bkg->GetBinError(i)+bdt_bkg->GetBinContent(i)*sys_bkg->GetBinContent(i);
      bdt_sig_full->SetBinError(i, sig_err);
      bdt_bkg_full->SetBinError(i, bkg_err);
    }
    

    std::cout << bdt_sig->GetNbinsX() << std::endl;
    std::cout << bdt_bkg->GetNbinsX() << std::endl;
    std::cout << sys_sig->GetNbinsX() << std::endl;
    std::cout << bdt_bkg->GetNbinsX() << std::endl;
    std::cout << "[ ";
    for ( int i = 1; i < bdt_sig->GetNbinsX()+1; i++ ) {
      std::cout << bdt_sig->GetBinLowEdge(i);
      if ( i < bdt_sig->GetNbinsX() ) {
	std::cout << ", ";
      }
    }
    std::cout << " ]" << std::endl;
    cout_bins(bdt_sig);
    cout_bins(bdt_bkg);
    cout_bins(sys_sig);
    cout_bins(sys_bkg);

    std::cout << std::endl;
    cout_bins(bdt_sig, sys_sig);
    cout_bins(bdt_bkg, sys_bkg);

    TCanvas * c1 = new TCanvas();

    gStyle->SetOptStat(0);

    TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);

    bdt_bkg_full->SetLineColor(kRed);
    bdt_bkg_full->SetTitle(Form("BDT scores (stat. + syst.) - %i MeV",mass));
    bdt_bkg->SetLineColor(kRed);
    leg->AddEntry(bdt_bkg_full,"Background");
    
    bdt_bkg_full->Draw("H X0 E1");
    bdt_bkg->Draw("H X0 E1 same");
    bdt_sig->SetMaximum(800);
    bdt_sig_full->SetMaximum(800);
    bdt_sig_full->SetMinimum(0);
    leg->AddEntry(bdt_sig_full,"Signal");
    
    bdt_sig_full->Draw("H X0 E1 same");
    bdt_sig->Draw("HX0E1 same");
    leg->Draw("same");
    gSystem->Exec(Form("mkdir -p img/%s/systematics",analysis.Data()));
    c1->SaveAs(Form("img/%s/systematics/BDT_scores_stat_syst_%imev.pdf",analysis.Data(),mass));
    c1->SaveAs(Form("img/%s/systematics/BDT_scores_stat_syst_%imev.png",analysis.Data(),mass));

    TFile *of = new TFile(Form("limit_cross_check_%imev.root",mass),"recreate");
    TH1F * sys_sig_abs = (TH1F*)bdt_sig_full->Clone();
    TH1F * sys_bkg_abs = (TH1F*)bdt_bkg_full->Clone();
    TH1F * stat_sig_abs = (TH1F*)bdt_sig_full->Clone();
    TH1F * stat_bkg_abs = (TH1F*)bdt_bkg_full->Clone();
    sys_sig_abs->Reset();
    sys_bkg_abs->Reset();
    stat_sig_abs->Reset();
    stat_bkg_abs->Reset();
    for ( int i = 1; i < sys_sig_abs->GetNbinsX()+1; i++ ) {
      sys_sig_abs->SetBinContent(i,bdt_sig_full->GetBinError(i));
      sys_bkg_abs->SetBinContent(i,bdt_bkg_full->GetBinError(i));
      stat_sig_abs->SetBinContent(i,bdt_sig->GetBinError(i));
      stat_bkg_abs->SetBinContent(i,bdt_bkg->GetBinError(i));
    }
    bdt_sig->SetName("signal");
    bdt_sig->Write();
    sys_sig->SetName("signal_DetVar_uncertainty_frac");
    sys_sig->Write();
    sys_sig_abs->SetName("signal_DetVar_unvertainty");
    sys_sig_abs->Write();
    stat_sig_abs->SetName("signal_stat_unvertainty");
    stat_sig_abs->Write();
    //
    bdt_bkg->SetName("background");
    bdt_bkg->Write();
    sys_bkg->SetName("background_DetVar_uncertainty_frac");
    sys_bkg->Write();
    sys_bkg_abs->SetName("background_DetVar_unvertainty");
    sys_bkg_abs->Write();
    stat_bkg_abs->SetName("background_stat_unvertainty");
    stat_bkg_abs->Write();
  }



  
  
}
