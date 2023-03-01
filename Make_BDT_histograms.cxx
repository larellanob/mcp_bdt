void Make_BDT_histograms(TString model = "model_230123_normal_thresholds")
{
  bool fVerbose {true};
  std::vector<int> masses {
    100,
    150,
    200,
    300,
    350,
    400
  };

  for ( auto mass: masses ) {
    if ( fVerbose ) std::cout << "\n\n\n\nMass \n\n\n\n" << mass << std::endl;
    
    TString filename = Form("root/%s/test_1000kev_2hits_%imev_BDT_scores.root",model.Data(),mass);
    TFile *f = new TFile(filename);
    
    TTree* bkg_tree = (TTree*)f->Get("test_bkg");
    TTree* sig_tree = (TTree*)f->Get("test_sig");
            
    TH1F *bkg_hist = new TH1F("bkg_hist","Test BDT bkg;Score;Entries",40,-10,10);
    TH1F *sig_hist = new TH1F("sig_hist","Test BDT sig;Score;Entries",40,-10,10);
    double bin_width = 0.5;
    
    bkg_tree->Draw("bdt>>bkg_hist");
    sig_tree->Draw("bdt>>sig_hist");
    
    TFile outfile(Form("hist/%s/no_bdtcut_s_b_mass_%i.root",model.Data(),mass),"recreate");
    // save the histograms, not the trees
    bkg_hist->Write();
    sig_hist->Write();

    TCanvas c3;
    TH1F *h_base = new TH1F("h_base","Test BDT bkg;Score;Entries",40,-10,10);
    h_base->Draw();
    h_base->SetTitle(Form("Mass %i MeV BDT score;BDT score;Entries",mass));
    h_base->SetMaximum(std::max(sig_hist->GetMaximum(),bkg_hist->GetMaximum())*1.1);
    bkg_hist->SetLineColor(kRed);
    sig_hist->SetLineColor(kBlue);
    sig_hist->Draw("same");
    bkg_hist->Draw("same");
    sig_hist->SetMinimum(0);
    c3.SaveAs(Form("img/m%i_hist.pdf",mass));
  }
}
