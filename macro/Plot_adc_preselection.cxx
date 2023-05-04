void Plot_adc_preselection(TString filename)
{
  int threshold {1000};
  int nhits {2};
  int mass {100};
  TString det_thres = "";
  if ( filename.Contains("def-th") ) det_thres = "def-th";
  if ( filename.Contains("low-th") ) det_thres = "low-th";
  TString outdir;
  TString outfile;
  TObjArray *tok = filename.Tokenize("/");
  for ( int i = 0; i < tok->GetEntries(); i++ ) {
    TString t = ((TObjString*)(tok->At(i)))->String();
    if ( i < tok->GetEntries()-1 ) outdir += (t+"/");
    if ( i == tok->GetEntries()-1 ) outfile = t;
    if ( t.Contains("_") ) {
      TObjArray *tok2 = t.Tokenize("_");
      for ( int j = 0; j < tok2->GetEntries(); j++ ) {
	TString t2 = ((TObjString*)(tok2->At(j)))->String();
	TString t_temp = t2;
	if ( t_temp.Contains("kev") ) {
	  t_temp.ReplaceAll("kev","");
	  threshold = t_temp.Atoi();
	} else if ( t_temp.Contains("mev") ) {
	  t_temp.ReplaceAll("mev","");
	  mass = t_temp.Atoi();
	} else if ( t_temp.Contains("hits") ) {
	  t_temp.ReplaceAll("hits","");
	  nhits = t_temp.Atoi();
	}
      }
    }
  }
  outdir.ReplaceAll("root","img");
  outdir.ReplaceAll("systematics","img");
  outfile.ReplaceAll(".root",".pdf");
  outfile.ReplaceAll("spacepoints","spacepoints_adc");
  std::cout << "outdir, outfile " << outdir << " " << outfile << std::endl;
  TFile *f = new TFile(filename);

  // POT information
  double tot_pot = 0;
  double pot = 0;
  TTree *pot_tree = (TTree*)f->Get("ana/pottree");
  if ( pot_tree != nullptr ) {
    pot_tree->SetBranchAddress("totpot",&pot);
    for(int ievent = 0; ievent < pot_tree->GetEntries(); ++ievent) {
      pot_tree->GetEntry(ievent);
      tot_pot += pot;
    }
  }
  
  TTree *tree = (TTree*)f->Get("ana/t");
  if ( tree == nullptr ) {
    std::cout << "WARNING: using old gallery script instead of analyzer\n"
	      << "No POT information will be saved" << std::endl;
    tree = (TTree*)f->Get("t");
  }
  int tot_evt = (int)tree->GetEntries();
  auto c1 = new TCanvas();
  c1->SetLogy();
  gStyle->SetOptStat(0);

  TH1F * h0 = new TH1F("h0","",100,0,4);
  TH1F * bkg = new TH1F("bkg",";Log_{10}(plane 2 ADC)",100,0,4);
  bkg->SetLineColor(kRed);
  TH1F * sig = new TH1F("sig","",100,0,4);

  auto n_bkg = tree->Draw("log10(pl2_integs)>>bkg","run > 0");
  h0 = (TH1F*)bkg->Clone(); // get correct vertical scale for histo
  TString parameters = Form("%s %ikeV %ihits %iMeV",det_thres.Data(),threshold, nhits,mass);
  h0->SetTitle(Form("%s;Log_{10}(plane 2 ADC)",parameters.Data()));
  h0->Draw();
  tree->Draw("log10(pl2_integs)>>bkg","run > 0","same");
  auto n_sig = tree->Draw("log10(pl2_true_integs)>>sig","log10(pl2_true_integs)!=0","same");

  TLegend *myleg = new TLegend(0.12,0.65,0.35,0.9);
  //  bkg->SetTitle(Form("bkg (%lli pairs)",n_bkg));
  //  sig->SetTitle(Form("sig (%lli pairs)",n_sig));
  //myleg = c1->BuildLegend(0.1,0.7,0.3,0.9,"","");
  myleg->AddEntry("bkg",Form("bkg (%lli pairs)",n_bkg));
  myleg->AddEntry("sig",Form("sig (%lli pairs)",n_sig));
  myleg->AddEntry("",Form("MC POT: %.3e",tot_pot));
  myleg->AddEntry("",Form("MC Events: %i",tot_evt));
  myleg->Draw();

  gSystem->Exec("mkdir -p "+outdir);
  c1->SaveAs(outdir+outfile);
  outfile.ReplaceAll(".pdf",".png");
  c1->SaveAs(outdir+outfile);

  std::cout << "calculating best sensitivity" << std::endl;
  double range_bkg = 0;
  double range_sig = 0;
  double sensitivity = 0;
  int best_low = 0;
  int best_high = 0;
  for ( int low = 1; low < 101; low++ ) {
    for ( int high = low; high < 100; high++ ) {
      range_bkg = bkg->Integral(low,high);
      range_sig = sig->Integral(low,high);
      //std::cout << "low, high, range bkg range sig:\n";
      //std::cout << low << " " << high << " " <<  range_bkg << " " <<  range_sig << std::endl;
      if ( range_bkg == 0 ) continue;
      double sens = range_sig/sqrt(range_bkg);
      if ( sens > sensitivity ) {
	sensitivity = sens;
	best_low = low;
	best_high = high;
      }
      
    }
  }
  std::cout << Form("Sensitivity %.3f, for low %.3f, high %.3f",sensitivity, sig->GetBinLowEdge(best_low), sig->GetBinLowEdge(best_high+1)) << std::endl;
  std::cout << Form("sig %.3f, bkg %.3f\n", sig->Integral(best_low,best_high), bkg->Integral(best_low,best_high));
}
