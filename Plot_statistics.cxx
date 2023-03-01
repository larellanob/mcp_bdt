void Draw_Graphs(std::vector<TGraph *> gr,TLegend * tleg,TString xaxis, TString filename);
void Plot_statistics(TString model = "model_221207_combinedmass_pot")
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

  std::vector<TGraph *> graphs_eff;
  std::vector<TGraph *> graphs_pur;
  std::vector<TGraph *> graphs_s_b;

  bool b_reverse = true;
  
  for ( auto mass: masses ) {
    if ( fVerbose ) std::cout << "\n\n\n\nMass \n\n\n\n" << mass << std::endl;

    TString filename = Form("root/%s/test_1000kev_2hits_%imev_BDT_scores.root",model.Data(),mass);
    TFile *f = new TFile(filename);
    //TFile f(filename);
    
    TTree* bkg_tree = (TTree*)f->Get("test_bkg");
    TTree* sig_tree = (TTree*)f->Get("test_sig");
    
        
    TH1F *bkg_hist = new TH1F("bkg_hist","Test BDT bkg;Score;Entries",40,-10,10);
    TH1F *sig_hist = new TH1F("sig_hist","Test BDT sig;Score;Entries",40,-10,10);
    double bin_width = 0.5;
    
    bkg_tree->Draw("bdt>>bkg_hist");
    sig_tree->Draw("bdt>>sig_hist");
    
    int last_nonempty_bin_bkg {0};
    int last_nonempty_bin_sig {0};
    int first_nonempty_bin_bkg {0};
    int first_nonempty_bin_sig {0};
    
    
    bool leading_empty_bins_bkg {true};
    bool leading_empty_bins_sig {true};
    bool trailing_empty_bins_bkg {false};
    bool trailing_empty_bins_sig {false};
    for ( int i = 0; i < 40; i++ ) {
      if ( sig_hist->GetBinContent(i) == 0 && !trailing_empty_bins_sig ) first_nonempty_bin_sig++;
      if ( bkg_hist->GetBinContent(i) == 0 && !trailing_empty_bins_bkg ) first_nonempty_bin_bkg++;
      if ( bkg_hist->GetBinContent(i) != 0 ) {
	leading_empty_bins_bkg = false;
	trailing_empty_bins_bkg = true;
      }
      if ( sig_hist->GetBinContent(i) != 0 ) {
	leading_empty_bins_sig = false;
	trailing_empty_bins_sig = true;
      }
      if ( bkg_hist->GetBinContent(i) != 0 || leading_empty_bins_bkg ) last_nonempty_bin_bkg++;
      if ( sig_hist->GetBinContent(i) != 0 || leading_empty_bins_sig ) last_nonempty_bin_sig++;
    }
    double last_nonempty_bin = max(last_nonempty_bin_bkg,last_nonempty_bin_sig);
    double first_nonempty_bin = min(first_nonempty_bin_bkg,first_nonempty_bin_sig);
    double low_edge = sig_hist->GetBinLowEdge(first_nonempty_bin);
    double up_edge = bkg_hist->GetBinLowEdge(last_nonempty_bin);

    if ( fVerbose ) {
      std::cout << Form("Low edge %.1f, up edge: %.1f, 1st nonempty bin: %.1f, last non0 bin: %.1f ",
			low_edge, up_edge,first_nonempty_bin, last_nonempty_bin) << std::endl;

      std::cout << Form("last bkg %i, last sig %i, 1st bkg %i, 1st sig %i",
			last_nonempty_bin_bkg,last_nonempty_bin_sig,
			first_nonempty_bin_bkg,first_nonempty_bin_sig)
		<< std::endl;
      
    }
    
    TGraph * gr_eff = new TGraph();
    TGraph * gr_pur = new TGraph();
    TGraph * gr_s_b = new TGraph();
    Double_t cut_score {-6};

    /////////////////////////////////////////////////
    // find optimum BDT score with metric s/sqrt(s+b)
    /////////////////////////////////////////////////

    double max_s_b{0};
    while ( cut_score < up_edge ) {
      // bdt cut histograms
      int nbins = (up_edge-cut_score)/bin_width;
      TH1F *bkg_cut = new TH1F("bkg_cut","Test BDT bkg (cut);Score;Entries",nbins,cut_score,up_edge);
      TH1F *sig_cut = new TH1F("sig_cut","Test BDT sig (cut);Score;Entries",nbins,cut_score,up_edge);
      // fill bdt cut histograms
      bkg_tree->Draw("bdt>>bkg_cut",Form("bdt>%.2f",cut_score));
      sig_tree->Draw("bdt>>sig_cut",Form("bdt>%.2f",cut_score));

      if ( b_reverse ) {
	delete sig_cut;
	delete bkg_cut;
	nbins = (cut_score-low_edge)/bin_width;
	std::cout << "cut score, lowedge, nbins: " << cut_score << " " << low_edge << " " << nbins << std::endl;
	bkg_cut = new TH1F("bkg_cut","Test BDT bkg (cut);Score;Entries",nbins,low_edge,cut_score);
	sig_cut = new TH1F("sig_cut","Test BDT sig (cut);Score;Entries",nbins,low_edge,cut_score);
	sig_tree->Clear();
	bkg_tree->Clear();
	bkg_tree->Draw("bdt>>bkg_cut",Form("bdt<%.2f",cut_score));
	sig_tree->Draw("bdt>>sig_cut",Form("bdt<%.2f",cut_score));
      }
      
      int bkg_n = bkg_tree->GetEntries();
      int sig_n = sig_tree->GetEntries();
      int fp = bkg_cut->GetEntries();
      int tp = sig_cut->GetEntries();
      int tn = bkg_n - fp;
      int fn = sig_n - tp;

      if ( b_reverse && ((tp == 0 && fp == 0) || nbins <= 0)) {
	cut_score+= 0.5;
	delete bkg_cut;
	delete sig_cut;
	continue;
      }
      
      double eff = (float)tp/((int)(tp+fn));
      double pur = (float)tp/((int)(tp+fp));
      //double s_spb = (double)tp/(sqrt((double)(tp+fp))*(bkg_n+sig_n)); // 'normalized'
      double s_spb = (double)tp/sqrt((double)(tp+fp)); // 'unnormalized'

      if ( fVerbose ) {
	std::cout << Form("Cut %.1f, Efficiency %.3f, Purity %.3f, s/sqrt(s+b) %.3f",cut_score,eff,pur,s_spb) << std::endl;
	//std::cout << Form("\tTP %i, FP %i, TN %i, FN %i",tp,fp,tn,fn) << std::endl;
      }

      gr_eff->AddPoint(cut_score,eff);
      gr_pur->AddPoint(cut_score,pur);
      gr_s_b->AddPoint(cut_score,s_spb); // 'unnormalized'
      //std::cout << (double)tp/(sqrt((double)(tp+fp)*(bkg_n+sig_n))) << std::endl;


      // save histo to root file if it's maximum s/sqrt(s+b)
      //if ( max_s_b < TMath::MaxElement(gr_s_b->GetN(),gr_s_b->GetY()) && mass == 350 ) {
      if ( max_s_b < TMath::MaxElement(gr_s_b->GetN(),gr_s_b->GetY()) ) {
	max_s_b = TMath::MaxElement(gr_s_b->GetN(),gr_s_b->GetY());
	TFile outfile(Form("hist/%s/post_bdtcut_s_b_mass_%i.root",model.Data(),mass),"recreate");
	bkg_cut->Write();
	sig_cut->Write();

	TCanvas c3;
	sig_cut->SetTitle(Form("Mass %i MeV after BDT score selection (%.1f);BDT score;Entries",mass,cut_score));
	bkg_cut->SetLineColor(kRed);
	sig_cut->SetLineColor(kBlue);
	sig_cut->Draw();
	bkg_cut->Draw("same");
	sig_cut->SetMinimum(0);
	c3.SaveAs(Form("img/m%i_cut0.pdf",mass));
	f->cd();
      }

      cut_score+= 0.5;
      delete bkg_cut;
      delete sig_cut;
    }
    graphs_eff.push_back(gr_eff);
    graphs_pur.push_back(gr_pur);
    graphs_s_b.push_back(gr_s_b);
    //delete sig_tree;
    //delete bkg_tree;
    //delete sig_hist;
    //delete bkg_hist;
    std::cout << "deleted garbage" << std::endl;
    //f->Close();
    //delete f;  
  }



  // plot efficiency, purity, s/sqrt(s+b)
  TLegend *leg = new TLegend(0.35,0.1,0.65,0.4);
  //void Draw_Graphs(std::vector<TGraph *> gr,TLegend * tleg,TString xaxis, TString filename);
  Draw_Graphs(graphs_eff,leg,"Efficiency","graphs_eff");
  //leg->SetX1NDC(0.3);
  //leg->SetX2NDC(0.6);
  Draw_Graphs(graphs_pur,leg,"Purity","graphs_pur");
  Draw_Graphs(graphs_s_b,leg,"s/#sqrt{s+b}","graphs_s_b");
  delete leg;

  // mass 350 = masses[4]
  /*
  TFile *f2 = new TFile("hist/post_bdtcut_s_b.root");
  TH1F* bkg_cut = (TH1F*)f2->Get("bkg_cut");
  TH1F* sig_cut = (TH1F*)f2->Get("sig_cut");

  // bkg only
  double likelihood {0};
  double poisson {0};
  if ( !b_reverse ) {
    for ( int bin{1}; bin < bkg_cut->GetNbinsX()+1; bin++ ) {
      double expected = bkg_cut->GetBinContent(bin);
      double b = expected;
      poisson = 0;
      double sum {0};
      for ( int n = 0; n <= expected; n++ ) {
	sum += pow(b,n)/TMath::Factorial(n);
      }
      poisson = exp(-b)*sum;
      std::cout << expected << " " << sum << " " << poisson<< std::endl;
    }
  } else {
    std::cout << "finished" << std::endl;
  }
  */
  std::cout << "finished" << std::endl;
}

void Draw_Graphs(std::vector<TGraph *> gr,TLegend * leg,TString xaxis, TString filename)
{
  std::vector<int> masses {
    100,
    150,
    200,
    300,
    350,
    400
  };
  gStyle->SetOptStat(0);
  TString hist_title = Form(";BDT score Cut;%s",xaxis.Data());
  TH1F h1("h1",hist_title,100,-5.5,5);
  TCanvas c1;
  h1.Draw();

  leg->Clear();
  double max {0};
  for ( int i{0}; i < gr.size(); i++ ) {
    gr[i]->Draw("same");
    max = std::max(max,TMath::MaxElement(gr[i]->GetN(),gr[i]->GetY()));
    gr[i]->SetMarkerStyle(i+20);
    gr[i]->SetMarkerColor(i+2);
    gr[i]->SetLineColor(i+2);
    gr[i]->SetTitle(Form("M = %i MeV",masses[i]));
    leg->AddEntry(gr[i]);
  }
  h1.SetMaximum(max*1.1);
  leg->Draw();
  c1.SaveAs("img/"+filename+".pdf");

}
