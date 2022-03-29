void Plot_BDT(TString model = "my_model_1")
{
  TFile *f = new TFile(Form("root/%s_BDT_scores.root",model.Data()));

  // trees
  // tbg_train
  // tbg_test
  // tsig_train
  // tsig_test
  auto c1 = new TCanvas();
  TH1F *h1 = new TH1F("h1","Test BDTs;Score;Entries",40,-10,10);
  TH1F *h2 = new TH1F("h2","Test BDTs;Score;Entries",40,-10,10);
  TH1F *h3 = new TH1F("h3","Train bkg;Score;Entries",40,-10,10);
  TH1F *h4 = new TH1F("h4","Train sig;Score;Entries",40,-10,10);

  
  TTree* bkg = (TTree*)f->Get("test_bkg");
  TTree* sig = (TTree*)f->Get("test_sig");
  TTree* bkg_train = (TTree*)f->Get("train_bkg");
  TTree* sig_train = (TTree*)f->Get("train_sig");

  
  bkg->Draw("bdt>>h1","","hist same");
  sig->Draw("bdt>>h2","","hist same");

  h1->SetMaximum(80);
  h3->SetMarkerStyle(kFullTriangleDown);
  h3->SetMarkerColor(kRed);
  h1->SetLineColor(kRed);
  bkg_train->Draw("bdt>>h3","","P same");
  h3->SetTitle("Train bkg (scaled 0.4)");

  h3->Scale(0.4);
  h4->SetMarkerStyle(kFullTriangleUp);
  h4->SetMarkerColor(kBlue);

  sig_train->Draw("bdt>>h4","","P same");
  h4->Scale(0.4);
  h4->SetTitle("Train sig (scaled 0.4)");

  h1->SetTitle("Test Background");
  h2->SetTitle("Test Signal");
  TLegend *myleg = new TLegend();
  myleg = c1->BuildLegend(0.1,0.7,0.3,0.9,"","");
  //myleg->Draw();
  h1->SetTitle("Test BDTs");  
  c1->SaveAs(Form("img/%s_BDT_score.pdf",model.Data()));
  
}
