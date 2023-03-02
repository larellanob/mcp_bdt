void Plot_Raw(int mass = 0)
{
  int threshold {1000};
  int nhits {2};
  TString model = "230123_normal_thresholds";
  TString filename = Form("root/model_%s/spacepoints_%ikev_%ihits_%imev.root",model.Data(),threshold,nhits,mass);
  if ( mass == 0 ) {
    std::cout << "using combined mass " << std::endl;
    filename = Form("ntuple/spacepoints_%ikev_%ihits_combinedmass.root",threshold,nhits);
  }
  TFile *f = new TFile(filename);
  //TH2F * full_angles = new TH2F("full_angles",
  TTree *tree = (TTree*)f->Get("t");
  auto c1 = new TCanvas();
  c1->SetLogy();

  TH1F * bkg = new TH1F("bkg",";Log_{10}(plane 2 ADC)",100,0,4);
  bkg->SetLineColor(kRed);

  TH1F * sig = new TH1F("sig","",100,0,4);
  

  tree->Draw("log10(pl2_integs)>>bkg");
  tree->Draw("log10(pl2_true_integs)>>sig","log10(pl2_true_integs)!=0","same");

  TLegend *myleg = new TLegend();
  myleg = c1->BuildLegend(0.1,0.7,0.3,0.9,"","");

  if ( mass != 0 ) {
    c1->SaveAs(Form("img/Plane_integrals_%ikeV_%ihits_%iMeV.pdf",threshold,nhits,mass));
  } else if ( mass == 0 ) {
    c1->SaveAs(Form("img/Plane_integrals_%ikeV_%ihits_combinedmass.pdf",threshold,nhits));
  }

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
