void Make_BDT_histograms(TString model_sample,TString sample = "", bool b_data = true)
{
  
  bool fVerbose {false};
  TFile *f = new TFile(model_sample);
  /*
  TTree *pot;
  pot = (TTree*)f->Get("total_pot");
  */
  TTree* bkg_tree;
  TTree* sig_tree;
  bkg_tree = (TTree*)f->Get("train_bkg");
  sig_tree = (TTree*)f->Get("train_sig");
  
  TH1F *bkg_hist = new TH1F("bkg_hist","Test BDT bkg;Score;Entries",40,-10,10);
  TH1F *sig_hist = new TH1F("sig_hist","Test BDT sig;Score;Entries",40,-10,10);
  TH1F *data_hist = new TH1F("data_hist","Test BDT sig;Score;Entries",40,-10,10);

  bkg_tree->Draw("bdt>>bkg_hist","","goff");
  sig_tree->Draw("bdt>>sig_hist","","goff");
  TTree * data_tree;
  if ( b_data ) {
    TFile * fdata= new TFile(sample);
    //fdata->cd();
    data_tree = (TTree*)fdata->Get("tdata");
    //std::cout << "data entries1:  "  << data_tree->GetEntries() << std::endl;
    f->cd();
    data_tree->Draw("bdt>>data_hist","","goff");
    std::cout << "data entries:  "  << data_hist->GetEntries() << std::endl;
  }

  TString out_filename = sample;
  out_filename.ReplaceAll("BDT_scores","BDT_scores_hist");
  TFile outfile(out_filename,"recreate");
  // save the histograms, not the trees
  bkg_hist->Write();
  sig_hist->Write();
  /*
  TTree * newpottree;
  newpottree = pot->CloneTree();
  newpottree->Write();
  */
  TCanvas c3;
  gStyle->SetOptStat(0);
  TH1F *h_base = new TH1F("h_base","Test BDT bkg;Score;Entries",40,-10,10);
  h_base->Draw();
  h_base->SetMaximum(std::max(sig_hist->GetMaximum(),bkg_hist->GetMaximum())*1.1);
  bkg_hist->SetLineColor(kRed);
  sig_hist->SetLineColor(kBlue);
  sig_hist->Draw("same");
  bkg_hist->Draw("same");
  sig_hist->SetMinimum(0);

  if ( b_data ) {
    //data_hist->Scale(0.12);
    data_hist->Scale(1.2);
    data_hist->SetMarkerStyle(kFullCircle);
    data_hist->SetMarkerSize(0.5);
    data_hist->Draw("same E0");
    for ( int i = 0; i < data_hist->GetNbinsX(); i++ ) {
      std::cout << i << " " << data_hist->GetBinContent(i) << std::endl;
    }
		

  }
  
  TString out_pdffile = out_filename;
  out_pdffile.ReplaceAll("samples/","img/");
  out_pdffile.ReplaceAll(".root",".pdf");
  c3.SaveAs(out_pdffile);

}
