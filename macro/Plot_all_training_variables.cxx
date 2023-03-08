void Plot_all_training_variables(int mass = 100)
{
  int threshold {1000};
  int nhits {2};
  TString sps_or_pre = "preselection";
  
  TString filename = Form("ntuple/%s_%ikev_%ihits_%imev.root",sps_or_pre.Data(),threshold,nhits,mass);
  if ( mass == 0 ) {
    std::cout << "using combined mass " << std::endl;
    filename = Form("ntuple/spacepoints_%ikev_%ihits_combinedmass.root",threshold,nhits);
  }
  TFile *f = new TFile(filename);
  //TH2F * full_angles = new TH2F("full_angles",
  TTree *tree;
  auto c1 = new TCanvas();
  //c1->SetLogy();

  std::vector<TString> trees_v;
  
  if ( sps_or_pre == "spacepoints" ){
    trees_v.push_back("t");

  }

  if ( sps_or_pre == "preselection" ) {
    trees_v.push_back("tbg_train");
    trees_v.push_back("tbg_test");
    trees_v.push_back("tsig_train");
    trees_v.push_back("tsig_test");
  }
  
  gSystem->Exec(Form("mkdir -p img/%s_%ikev_%ihits_%imev",sps_or_pre.Data(),threshold,nhits,mass));
  for ( auto t: trees_v ) {
    tree = (TTree*)f->Get(t);
    int n_branches = tree->GetListOfBranches()->GetEntries();
    for ( int i = 0; i < n_branches; i++ ) {
      TString b_name = tree->GetListOfBranches()->At(i)->GetName();
      //std::cout << i << " " << tree->GetListOfBranches()->At(i)->GetName() << std::endl;
      std::cout << i << " " << b_name << std::endl;
      tree->Draw(b_name);
      c1->SaveAs(Form("img/%s_%ikev_%ihits_%imev/%s_%i_%s.png",
		      sps_or_pre.Data(),threshold,nhits,mass,t.Data(),i,b_name.Data()));
    }
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

  
  if ( mass != 100 ) {
    c1->SaveAs(Form("img/Plane_integrals_%ikeV_%ihits_%iMeV.pdf",threshold,nhits,mass));
  } else if ( mass == 0 ) {
    c1->SaveAs(Form("img/Plane_integrals_%ikeV_%ihits_combinedmass.pdf",threshold,nhits));
  }


}
