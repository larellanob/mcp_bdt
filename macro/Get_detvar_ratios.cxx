void Draw_detvar_histos(std::vector<TH1F> ratios,TString label,TString dir);

  std::vector<TString> samples = {
    "CV",
    "wiremod_x",
    "wiremod_yz",
    "wiremod_anglexz",
    "wiremod_angleyz",
    "wiremod_dEdx",
    "ly_down",
    "ly_rayleigh",
    "ly_atten"
    //"sce",
    //"recomb"
  };
void Get_detvar_ratios(TString dir)
{

  TFile *f_cv = new TFile(dir+"/hist_BDT_scores_CV.root");
  TH1F *cv_sig = (TH1F*)f_cv->Get("sig_hist");
  TH1F *cv_bkg = (TH1F*)f_cv->Get("bkg_hist");

  // original histograms are 40 bins from -10 to 10, each bin 2.5 width
  // as I'm looking at it now, the histograms range from -5.0,2.5 (i.e. 
  TH1F *base_hist =  new TH1F("base_hist",";BDT Score;Entries",40,-10,10);
  
  std::vector<TH1F> sig_ratios;
  std::vector<TH1F> bkg_ratios;

  
  for ( auto s: samples ) {
    // open detvar file
    std::cout << "Opening " << s << " file" << std::endl;
    TFile f_dv = TFile(dir+"/hist_BDT_scores_"+s+".root");
    // get signal
    TH1F *h_dv_sig = (TH1F*)base_hist->Clone();
    h_dv_sig = (TH1F*)f_dv.Get("sig_hist");

    // get bkg
    TH1F *h_dv_bkg = (TH1F*)base_hist->Clone();
    h_dv_bkg = (TH1F*)f_dv.Get("bkg_hist");

    //
    TH1F cv_dv_sig = 100*((*cv_sig)-(*h_dv_sig))/(*cv_sig);
    TH1F cv_dv_bkg = 100*((*cv_bkg)-(*h_dv_bkg))/(*cv_bkg);

    sig_ratios.push_back(cv_dv_sig);
    bkg_ratios.push_back(cv_dv_bkg);
    
    
    
  }
  dir.ReplaceAll("systematics","img");
  Draw_detvar_histos(sig_ratios,"signal",dir);
  Draw_detvar_histos(bkg_ratios,"background",dir);
}
void Draw_detvar_histos(std::vector<TH1F> ratios,TString label,TString dir)
{
  auto c1 = new TCanvas();
  if ( ratios.size() < 1 ) return;
  gStyle->SetOptStat(0);
  TH1F *base_hist = (TH1F*)ratios[0].Clone();
  base_hist->Reset();
  base_hist->Draw();
  base_hist->SetTitle(Form("Detector variations (%s);BDT Score;100x(CV-detvar)/CV",label.Data()));
  base_hist->SetMinimum(-100);
  base_hist->SetMaximum(100);
  // horiz line 25%
  TLine *line25p = new TLine(-10,25,10,25);
  TLine *line25m = new TLine(-10,-25,10,-25);
  line25p->SetLineStyle(7);
  line25p->SetLineColor(kGreen+1);
  line25p->SetLineWidth(2);
  line25p->Draw();
  line25m->SetLineStyle(7);
  line25m->SetLineColor(kGreen+1);
  line25m->SetLineWidth(2);
  line25m->Draw();
  // horiz line 50%
  TLine *line50p = new TLine(-10,50,10,50);
  TLine *line50m = new TLine(-10,-50,10,-50);
  line50p->SetLineStyle(7);
  line50p->SetLineColor(kOrange+1);
  line50p->SetLineWidth(2);
  line50p->Draw();
  line50m->SetLineStyle(7);
  line50m->SetLineColor(kOrange+1);
  line50m->SetLineWidth(2);
  line50m->Draw();
  

  // quadrature
  TH1F * h_quad = (TH1F*)base_hist->Clone();
  TH1F * h_quad_m = (TH1F*)base_hist->Clone();
  for ( int i = 1; i <= base_hist->GetNbinsX(); i++ ) {
    double quad = 0;
    for ( int s = 0; s < ratios.size(); s++ ){
      quad += pow(ratios[0].GetBinContent(i)-ratios[s].GetBinContent(i),2.);
    }
    quad = sqrt(quad);
    h_quad->SetBinContent(i,quad);
    h_quad_m->SetBinContent(i,-quad);
  }


  TLegend leg(0.7,0.1,0.9,0.7);
  leg.AddEntry(line25p,"25%","L");
  leg.AddEntry(line50p,"50%","L");
  for ( int i = 0; i < ratios.size(); i++ ) {
    ratios[i].SetLineColor(1+i);
    ratios[i].Draw("HE0 same");
    leg.AddEntry(&ratios[i],samples[i]);
    
  }
  h_quad->SetLineWidth(2);
  h_quad->SetLineColor(kBlack);
  h_quad_m->SetLineWidth(2);
  h_quad_m->SetLineColor(kBlack);
  leg.AddEntry(h_quad,"Quadrature sum");
  h_quad->Draw("same");
  h_quad_m->Draw("same");
  leg.Draw("same");
  
  c1->SaveAs(Form("%s/diffs_%s.pdf",dir.Data(),label.Data()));
}
