void Plot_normalized_signal(TString mass = "200")
{
  gStyle->SetOptStat(0);

  // normalization type

  // efficiency normalization: randomly fill to histograms according
  // to number of missing reco blips
  Bool_t do_eff_norm = false;
  // pot normalization: normalize according to signal POT
  Bool_t do_pot_norm = false;
  // tail normalization: normalize according to integral of the tail
  // of the histogram, where low-threshold stuff shouldn't matter
  Bool_t do_tail_norm = true;
  float tail_value = 2.5;
  
  std::vector<TString> gen_thresholds =
    {
      "100",
      "500",
      "1000"
    };
  TH1F * h0 = new TH1F("h0","",100,0,4);
  TH1F * h1 = new TH1F("h1","",100,0,4);
  TH1F * h2 = new TH1F("h2","",100,0,4);
  TH1F * h3 = new TH1F("h3","",100,0,4);

  const char *blips[4] = {"0","1", "2", "3+"};
  TH1F * eff1 = new TH1F("eff1","",4,-0.5,3.5);
  TH1F * eff2 = new TH1F("eff2","",4,-0.5,3.5);
  TH1F * eff3 = new TH1F("eff3","",4,-0.5,3.5);

  std::vector<TH1F*> h_v { h1, h2, h3 };
  std::vector<TH1F*> eff_v { eff1, eff2, eff3 };
  std::vector<float> scale_v;

  // insert empty data to get bin names in correct order
  for ( int i = 0; i < 3; i++ ) {
    eff_v[i]->Fill(blips[0],0);
    eff_v[i]->Fill(blips[1],0);
    eff_v[i]->Fill(blips[2],0);
    eff_v[i]->Fill(blips[3],0);
    eff_v[i]->Reset();
  }
  

  h1->SetLineColor(2);
  h2->SetLineColor(3);
  h3->SetLineColor(4);
  auto c1 = new TCanvas();
  c1->SetLogy();
  h0->SetMaximum(0.03);
  if ( mass == "400" ) h0->SetMaximum(0.003);
  if ( !do_pot_norm ) h0->SetMaximum(5000);
  TString basetitle = Form("True blips #leq 2 - %s MeV",mass.Data());
  TString fulltitle = Form("%s;log_{10}(true integrals);Reco blips",basetitle.Data());
  h0->SetTitle(fulltitle);
  h0->Draw();

  TLegend *myleg = new TLegend(0.12,0.65,0.35,0.9);
  
  int counter_signal=-1;
  for ( TString gth: gen_thresholds ) {
    counter_signal++;

    TFile *f;
    if ( gth == "1000" ) {
      f = new TFile(Form("root/def-th_%skev_2hits/spacepoints_1000kev_2hits_%smev.root",gth.Data(),mass.Data()));
    } else {
      f = new TFile(Form("root/def-th_%skev_2hits/spacepoints_%smev.root",gth.Data(),mass.Data()));
    }

    // get pot information
    double tot_pot = 0;
    double pot = 0;
    TTree *pot_tree = (TTree*)f->Get("ana/pottree");
    if ( pot_tree != nullptr ) {
      pot_tree->SetBranchAddress("totpot",&pot);
      for(int ievent = 0; ievent < pot_tree->GetEntries(); ++ievent) {
	pot_tree->GetEntry(ievent);
	tot_pot += pot;
      }
    } else { // 1000kev file doesn't have pot info
      std::cout << "WARNING: using old gallery script instead of analyzer\n"
		<< "Only hardcoded POT information will be available" << std::endl;
      if ( mass == "100" ) tot_pot = 1.47608e27;
      else if ( mass == "150" ) tot_pot = 2.43936e27;
      else if ( mass == "200" ) tot_pot = 2.55958e27;
      else if ( mass == "300" ) tot_pot = 4.76117e27;
      else if ( mass == "350" ) tot_pot = 9.75715e27;
      else if ( mass == "400" ) tot_pot = 2.28301e29;
    } // end getting pot
    std::cout << tot_pot << std::endl;
    double target_pot = 1.5e21;

    // get signal
    TTree *tree = (TTree*)f->Get("ana/t");
    if ( tree == nullptr ) {
      std::cout << "WARNING: using old gallery script instead of analyzer\n"
		<< "No POT information will be available" << std::endl;
      tree = (TTree*)f->Get("t");
    }
    TTreeReader tr(tree);
    TTreeReaderArray<double> true_integs(tr,"pl2_true_integs");
    TTreeReaderValue<int> run(tr,"run");
    TTreeReaderValue<int> evt(tr,"evt");
    
    int counter_events = 0;
    int n_zero_blips = 0;
    int n_one_blip = 0;
    int n_two_blips  = 0;
    int n_twoplus_blips = 0;
    
    while (tr.Next() ) {
      std::vector<float> adc;
      counter_events++;
      int current_nblips = 0;
      for ( int i = 0; i < true_integs.GetSize(); i++ ) {
	if ( true_integs[i] > 0.0 ) {
	  current_nblips++;
	  adc.push_back(log10(true_integs[i]));
	}
      }

      // signal distribution histogram
      for ( int i = 0; i < adc.size(); i++ ) {
	// ignore high blip events
	if ( adc.size() > 2 ) continue;
	h_v[counter_signal]->Fill(adc[i]);
      }

      // efficiency (bin migration histogram)
      if ( current_nblips == 0 ) {
	n_zero_blips++;
	eff_v[counter_signal]->Fill(blips[0],1);
      } else if ( current_nblips == 1 ) {
	n_one_blip++;
	eff_v[counter_signal]->Fill(blips[1],1);
      } else if ( current_nblips == 2 ) {
	n_two_blips++;
	eff_v[counter_signal]->Fill(blips[2],1);
      } else if ( current_nblips >2 ) {
	n_twoplus_blips++;
	eff_v[counter_signal]->Fill(blips[3],1);
      }
    } // end while ttreereader
    
    //tree->Draw("log10(pl2_true_integs)>>h","log10(pl2_true_integs)!=0","same");
    //TH1F *h = (TH1F*)gDirectory->Get("h");
    std::cout << "threshold: " << gth << std::endl;
    std::cout << counter_signal << " "  << " " << counter_events  << std::endl;
    std::cout << "number of 0,1,2,2+ blips: " << std::endl;
    std::cout << n_zero_blips << std::endl;
    std::cout << n_one_blip << std::endl;
    std::cout << n_two_blips << std::endl;
    std::cout << n_twoplus_blips << std::endl;
    std::cout << "--------------------------" << std::endl;
    h_v[counter_signal]->Draw("same");
    if ( do_pot_norm ) {
      h_v[counter_signal]->Scale(target_pot/tot_pot);
      h0->SetTitle(basetitle+"- POT norm.");
    }
    //if ( do_pot_norm ) h_v[counter_signal]->Scale(float(counter_events)/float(n_one_blip));
    if ( do_eff_norm ) {
      TH1F *hclone = (TH1F*)h_v[counter_signal]->Clone();
      h_v[counter_signal]->FillRandom(hclone,2*n_zero_blips+1*n_one_blip);
    }
    std::cout << "events: " << counter_events << " " << 10850./(float)counter_events << std::endl;
    int n_123 = n_zero_blips+n_one_blip+n_two_blips;
    std::cout << "0+1+2 : " << n_123 << " " << 10850./(float)n_123 << std::endl;
    float scaling_factor = 10850./(float)n_123;
    scale_v.push_back(scaling_factor);
    if ( do_eff_norm ) {
      h_v[counter_signal]->Scale(scaling_factor);
      std::cout << "scaling factor: " << scaling_factor << std::endl;
      myleg->AddEntry(h_v[counter_signal],Form("%s keV signal (%.0f)#times%.2f",
					       gen_thresholds[counter_signal].Data(),
					       h_v[counter_signal]->GetEntries(),
					       scaling_factor
					       ));
      

    }
  } // end looping over files/gen-thresholds


  // scale by the "tail" after 2
  // histograms contained within h_v[]
  if ( do_eff_norm ) {
    h0->SetTitle(basetitle+"- Eff. norm.");
  }
  if ( do_tail_norm ) {
    h0->SetTitle(basetitle+Form("- Normalized to integral of bins > %.1f",tail_value));
  }
  int starting_bin = -1;
  for ( int i = 0; i < h_v[0]->GetNbinsX()+1; i++ ) {
    std::cout << "hey " << h_v[0]->GetBinLowEdge(i) << std::endl;
    if ( h_v[0]->GetBinLowEdge(i) > tail_value ) {
      starting_bin = i;
      break;
    }
  }
  std::cout << "starting bin: " << starting_bin << std::endl;
  float tail_2 = h_v[2]->Integral(starting_bin,h_v[2]->GetNbinsX()+1);
  float tail_1 = h_v[1]->Integral(starting_bin,h_v[1]->GetNbinsX()+1);
  float tail_0 = h_v[0]->Integral(starting_bin,h_v[0]->GetNbinsX()+1);
  std::cout << Form("h_v[2] after %.1f: ",tail_value) << tail_2 << std::endl;
  std::cout << Form("h_v[1] after %.1f: ",tail_value) << tail_1 << std::endl;
  std::cout << Form("h_v[0] after %.1f: ",tail_value) << tail_0 << std::endl;

  std::vector<float> scaling_v { tail_2/tail_0, tail_2/tail_1, tail_2/tail_2};

  if ( do_tail_norm ) {
    for ( int i = 0; i < 3; i++ ) {
      h_v[i]->Scale(scaling_v[i]);
      myleg->AddEntry(h_v[i],Form("%s keV signal (%.0f)#times%.2f",
				  gen_thresholds[i].Data(),
				  h_v[i]->GetEntries(),
				  scaling_v[i]
				  ));
      
    }
    
    std::cout << "scale factors 0, 1: " << tail_2/tail_0 << " " << tail_2/tail_1 << std::endl;
    std::cout << "new tail 0, 1: " << h_v[0]->Integral(starting_bin,h_v[0]->GetNbinsX()+1)
	      << " " << h_v[1]->Integral(starting_bin,h_v[1]->GetNbinsX()+1) << std::endl;
  }
  
  myleg->Draw();
  if ( do_pot_norm ) c1->SaveAs(Form("img/norm-pot_signal_comparison_ttreereader_%smev.pdf",mass.Data()));
  if ( do_eff_norm ) c1->SaveAs(Form("img/norm-eff_signal_comparison_ttreereader_%smev.pdf",mass.Data()));
  if ( do_tail_norm ) c1->SaveAs(Form("img/norm-tail_signal_comparison_ttreereader_%smev.pdf",mass.Data()));
  

  gStyle->SetPadTopMargin(0.01);
  gStyle->SetPadBottomMargin(0.12);
  eff_v[0]->SetTitle(";Number of reconstructed blips; Events");
  eff_v[0]->SetLabelSize(0.07);
  eff_v[0]->SetLabelOffset(0.01);
  eff_v[0]->GetXaxis()->SetTitleOffset(1.8);
  //eff_v[0]->GetXaxis()->CenterTitle(true);
  //eff_v[0]->GetYaxis()->CenterTitle(true);
  TLegend *leg = new TLegend(0.4,0.85,0.9,1.0);
  auto *c2 = new TCanvas();
  //c2->SetLogy();
  for ( int i = 0; i < 3; i++ ) {
    std::cout << eff_v[i]->GetEntries() << std::endl;
    for ( int b = 1; b <= eff_v[i]->GetNbinsX(); b++ ) {
      std::cout << eff_v[i]->GetBinContent(b) << std::endl;
    }
    std::cout << "----------" << std::endl;
    eff_v[i]->SetMarkerStyle(20+i);
    eff_v[i]->SetMarkerColor(2+i);
    eff_v[i]->SetMarkerSize(2);
    eff_v[i]->SetLineColor(2+i);
    eff_v[i]->SetLineWidth(2);
    eff_v[i]->Draw("same PH");
    leg->AddEntry(eff_v[i],Form("Gen. threshold %skeV (%.0f events)",
				gen_thresholds[i].Data(),eff_v[i]->GetEntries()));
  }
  leg->Draw();
  c2->SaveAs("img/eff_combined.pdf");
  c2->SaveAs("img/eff_combined.png");
  
}
