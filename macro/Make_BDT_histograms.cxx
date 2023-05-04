void Make_BDT_histograms(TString filename,TString detvar="")
{
  if ( !filename.Contains("test") ) {
    std::cout << "Incorrect input file, needs to be test_ type" << std::endl;
    return;
  }
  TString filename_no_ext = filename;
  filename_no_ext.ReplaceAll(".root","");
  int mass;
  TObjArray *tok = filename.Tokenize("_");
  for ( int i = 0; i < tok->GetEntries(); i++ ) {
    TString t = ((TObjString*)(tok->At(i)))->String();
    if ( t.Contains("mev") ) {
      t.ReplaceAll("mev","");
      mass = t.Atoi();
    }
  }
  std::cout << "Mass: " << mass;
  bool fVerbose {false};

  
  TFile *f = new TFile(filename);

  TTree *pot;
  if ( detvar == "" ) {
    pot = (TTree*)f->Get("total_pot");
  }
  
  TTree* bkg_tree = (TTree*)f->Get("test_bkg");
  TTree* sig_tree = (TTree*)f->Get("test_sig");
  
  TH1F *bkg_hist = new TH1F("bkg_hist","Test BDT bkg;Score;Entries",40,-10,10);
  TH1F *sig_hist = new TH1F("sig_hist","Test BDT sig;Score;Entries",40,-10,10);
  double bin_width = 0.5;
  
  bkg_tree->Draw("bdt>>bkg_hist");
  sig_tree->Draw("bdt>>sig_hist");

  TString out_filename = filename;
  out_filename.ReplaceAll("test","hist");
  TFile outfile(out_filename,"recreate");
  // save the histograms, not the trees
  bkg_hist->Write();
  sig_hist->Write();
  TTree * newpottree;
  if ( detvar == "" ) {
    newpottree = pot->CloneTree();
    newpottree->Write();
  }
    
  TCanvas c3;
  gStyle->SetOptStat(0);
  TH1F *h_base = new TH1F("h_base","Test BDT bkg;Score;Entries",40,-10,10);
  h_base->Draw();
  h_base->SetTitle(Form("Mass %i MeV BDT score;BDT score;Entries",mass));
  if ( detvar != "" ) {
    h_base->SetTitle(Form("%iMeV - %s;BDT score;Entries",mass,detvar.Data()));
  }
  h_base->SetMaximum(std::max(sig_hist->GetMaximum(),bkg_hist->GetMaximum())*1.1);
  bkg_hist->SetLineColor(kRed);
  sig_hist->SetLineColor(kBlue);
  sig_hist->Draw("same");
  bkg_hist->Draw("same");
  sig_hist->SetMinimum(0);
  TString out_pdffile = out_filename;
  out_pdffile.ReplaceAll("root/","img/");
  out_pdffile.ReplaceAll("systematics/","img/");
  out_pdffile.ReplaceAll(".root",".pdf");
  tok = out_pdffile.Tokenize("/");
  gSystem->Exec("mkdir -p img/"+((TObjString*)(tok->At(1)))->String());
  c3.SaveAs(out_pdffile);

}
