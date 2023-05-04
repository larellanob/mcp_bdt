void Plot_training_variables(TString tag = "def-th_1000kev_2hits", int mass = 100)
{
  TString sps_or_pre = "preselection";
  TString filename = Form("root/%s/%s_%imev.root",tag.Data(),sps_or_pre.Data(),mass);

  TFile *f = new TFile(filename);
  auto c1 = new TCanvas();
  //c1->SetLogy();

  std::vector<TString> trees_v;
  
  if ( sps_or_pre == "spacepoints" ){
    trees_v.push_back("ana/t");

  }

  if ( sps_or_pre == "preselection" ) {
    //trees_v.push_back("tbg_train");
    trees_v.push_back("tbg_test");
    //trees_v.push_back("tsig_train");
    trees_v.push_back("tsig_test");
  }

  gSystem->Exec(Form("mkdir -p img/%s/%s_%imev",
		     tag.Data(),
		     sps_or_pre.Data(),
		     mass)
		);
  //TTree *tree = (TTree*)f->Get("ana/t");
  TTree *tree_sig = (TTree*)f->Get("tsig_test");
  TTree *tree_bkg = (TTree*)f->Get("tbg_test");
  int n_branches = tree_bkg->GetListOfBranches()->GetEntries();
  for ( int i = 0; i < n_branches; i++ ) {
    TString b_name = tree_bkg->GetListOfBranches()->At(i)->GetName();
    if ( b_name.Contains("nhit") ) {
      c1->SetLogx();
      c1->SetLogy();
    }
    else  {
      c1->SetLogx(0);
      c1->SetLogy(0);
    }
    //TString t_name = tree.Data();
    //t_name.ReplaceAll("/","-");
    std::cout << i << " " << b_name << std::endl;
    //if ( tree != nullptr ) tree->Draw(b_name);
    if ( tree_bkg != nullptr ) tree_bkg->Draw(b_name);
    if ( tree_sig != nullptr ) tree_sig->Draw(b_name,"","hist same");
    TH1F * bkg = ((TH1F*)(gPad->GetListOfPrimitives()->At(0)));
    TH1F * sig = ((TH1F*)(gPad->GetListOfPrimitives()->At(1)));



    bkg->SetFillColor(kBlue+3);
    sig->SetLineColor(kRed);
    sig->SetLineWidth(3);
    //float sig_scaling = bkg->GetMaximum()/(float)sig->GetMaximum();
    float sig_scaling = bkg->Integral()/(float)sig->Integral();
    TString h_title = Form("Mass %i - %s",mass,bkg->GetTitle());
    std::cout << h_title << " " << sig->GetMaximum() << std::endl;
    bkg->SetTitle(Form("Background - %lli entries",tree_bkg->GetEntries()));
    sig->SetTitle(Form("Signal - %lli entries (x%.2f)",tree_sig->GetEntries(),sig_scaling));
    sig->Scale(sig_scaling);
    bkg->SetMaximum(std::max(bkg->GetMaximum(),sig->GetMaximum())*1.2);
    if ( b_name.Contains("nhit") ) {
      bkg->SetMaximum(std::max(bkg->GetMaximum(),sig->GetMaximum())*10.5);
    }
    TLegend *myleg = new TLegend();
    myleg = c1->BuildLegend(0.5,0.8,0.9,0.9,"","");

    bkg->SetTitle(h_title.Data());
    
    std::cout << b_name.Data() << std::endl;
    c1->SaveAs(Form("img/%s/%s_%imev/%i_%s.png",
		    tag.Data(),
		    sps_or_pre.Data(),
		    mass,
		    //t_name.Data(),
		    i,
		    b_name.Data()
		    )
	       );
  }
  /*
  TH1F * bkg = new TH1F("bkg",";Log_{10}(plane 2 ADC)",100,0,4);
  bkg->SetLineColor(kRed);

  TH1F * sig = new TH1F("sig","",100,0,4);
  

  tree->Draw("log10(pl2_integs)>>bkg");
  tree->Draw("log10(pl2_true_integs)>>sig","log10(pl2_true_integs)!=0","same");

  TLegend *myleg = new TLegend();
  myleg = c1->BuildLegend(0.1,0.7,0.3,0.9,"","");

  */

  
  //c1->SaveAs(Form("img/%s/Plane_integrals_%imev.pdf",tag.Data(),mass));


}

void Plot_all_training_variables(TString tag = "def-th_1000kev_2hits")
{
  std::vector<int> masses {100,150,200,300,350,400};
  for ( auto m: masses ){
    Plot_training_variables(tag,m);
  }
}
