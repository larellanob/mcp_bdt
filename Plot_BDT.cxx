//void Plot_BDT(TString model = "test_model_221116_combinedmass_workshop")
void Plot_BDT(TString model = "model_230123_normal_thresholds", TString test_model = "")
{
  TFile *f = new TFile(Form("root/%s_BDT_scores.root",model.Data()));

  int test_mass;
  if ( test_model.IsDigit() ) {
    std::cout << "TEST MODEL IS DIGIT: " << test_model << std::endl;
    test_mass = test_model.Atoi();
    test_model = Form("1000kev_2hits_%imev_normal_thresholds",test_mass);
    std::cout << "TEST MODEL: " << test_model << std::endl;
  }
  bool bkg_nu;
  if ( test_model == "" ) {
    bkg_nu = false;
  } else {
    bkg_nu = true;
  }
  
  // trees
  // tbg_train
  // tbg_test
  // tsig_train
  // tsig_test
  auto c1 = new TCanvas();
  TH1F *h0 = new TH1F("h0","Test BDTs;Score;Entries",40,-6,6);
  TH1F *h1 = new TH1F("h1","Test BDTs;Score;Entries",40,-6,6);
  TH1F *h2 = new TH1F("h2","Test BDTs;Score;Entries",40,-6,6);
  TH1F *h3 = new TH1F("h3","Train bkg;Score;Entries",40,-6,6);
  TH1F *h4 = new TH1F("h4","Train sig;Score;Entries",40,-6,6);
  TH1F *h5 = new TH1F("h5","Train sig;Score;Entries",40,-6,6); // bkg_nu
  TH1F *h6 = new TH1F("h6","Train sig;Score;Entries",40,-6,6); // bkg_nu cosmics

  
  TTree* bkg = (TTree*)f->Get("test_bkg");
  TTree* sig = (TTree*)f->Get("test_sig");
  TTree* bkg_train = (TTree*)f->Get("train_bkg");
  TTree* sig_train = (TTree*)f->Get("train_sig");

  TFile *g;
  TTree *bkg_nu_test;
  TTree *bkg_nu_cosmics;

  if ( bkg_nu ) {
    // test model
    g = new TFile(Form("root/%s/test_%s_BDT_scores.root",model.Data(),test_model.Data()));
    //g = new TFile(Form("root/%s/test_neutrino_BDT_scores.root",model.Data()));
    //g = new TFile(Form("root/%s/test_1nu1cos_BDT_scores.root",model.Data()));
    bkg_nu_test = (TTree*)g->Get("test_sig");
    bkg_nu_cosmics = (TTree*)g->Get("test_bkg");
    f->cd();
  }

  h0->SetMaximum(1400);
  h0->Draw();
  if ( bkg_nu ) {
    std::cout << "doing nu" << std::endl;
    // neutrino
    bkg_nu_test->Draw("bdt>>h5","","hist same");
    double nu_scaling = 4.0;
    //h5->SetTitle(Form("#nu-overlay blip pair (Test, x%.1f)",nu_scaling));
    h5->SetTitle(Form("Normal thresholds mCP pairs (Test, x%.1f)",nu_scaling));
    h5->Scale(nu_scaling);
    h5->SetFillColor(kGreen);

    // cosmics
    //bkg_nu_cosmics->Draw("bdt>>h6","","P same");
    double ov_scaling = 0.65;
    h6->SetTitle(Form("Blips from overlay (Neutrino sample, x%.1f)",ov_scaling));
    h6->Scale(ov_scaling);
    h6->SetMarkerStyle(kFullSquare);
    //h6->SetMarkerColor(kGreen);
    //h6->SetFillColor(kPink);
  }
  bkg->Draw("bdt>>h1","","hist same");
  bkg_train->Draw("bdt>>h3","","P same");
  sig->Draw("bdt>>h2","","hist same");
  sig_train->Draw("bdt>>h4","","P same");
  

  //h1->SetMaximum(1400);
  h3->SetMarkerStyle(kFullTriangleDown);
  h3->SetMarkerColor(kRed);
  h1->SetLineColor(kRed);
  h3->SetTitle("Blips from overlay (Train, scaled 3/7.)");

  h3->Scale(3/7.);
  h4->SetMarkerStyle(kFullTriangleUp);
  h4->SetMarkerColor(kBlue);

    h4->Scale(3/7.);
  h4->SetTitle("Millicharge (Train, scaled 3/7.)");

  h1->SetTitle("Blips from overlay (Test)");
  h2->SetTitle("Millicharge (Test)");
  TLegend *myleg = new TLegend();
  myleg = c1->BuildLegend(0.1,0.6,0.4,0.9,"","");
  //myleg->Draw();
  h1->SetTitle("Test BDTs");
  if ( test_model != "" ) {
    c1->SaveAs(Form("img/%s/%s_BDT_score.pdf",model.Data(),test_model.Data()));
  } else {
    c1->SaveAs(Form("img/%s/%s_BDT_score.pdf",model.Data(),model.Data()));
  }
    
  
}
