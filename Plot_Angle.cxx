void Plot_Angle()
{
  TFile *f = new TFile("spacepoints.root");
  //TH2F * full_angles = new TH2F("full_angles",
  TTree *tree = (TTree*)f->Get("t");

  TH2F * bkg = new TH2F("bkg","angles background;#Delta#theta;#Delta#phi",50,0,90,50,-180,180);
  //  bkg->SetLineColor(kRed);

  TH2F * sig = new TH2F("sig","angles signal;theta;phi",50,0,90,50,-180,180);
  TH2F * bkg_and_sig = new TH2F("bkg_and_sig","angles bkg+sig;#Delta#theta;#Delta#phi",50,0,90,50,-180,180);

  TH1F * dist = new TH1F("dist","dist to true hits;distance",500,0,1200);

  vector<double> *pxs = 0, *pys = 0, *pzs = 0, *adcs = 0;
  vector<double> *true_xs = 0, *true_ys = 0, *true_zs = 0, *true_adcs = 0, *true_enes = 0;

  // branches in quotes exist within the input root file spacepoints.root/t
  tree->SetBranchAddress("sps_x",&pxs);
  tree->SetBranchAddress("sps_y",&pys);
  tree->SetBranchAddress("sps_z",&pzs);
  tree->SetBranchAddress("elec_x",&true_xs);
  tree->SetBranchAddress("elec_y",&true_ys);
  tree->SetBranchAddress("elec_z",&true_zs);
  tree->SetBranchAddress("elec_E",&true_enes);
  tree->SetBranchAddress("pl2_integs",&adcs);
  tree->SetBranchAddress("pl2_true_integs",&true_adcs);

  // true and direction (?) angles
  double ot_truetheta;
  double ot_truephi;
  double ot_dirtheta;
  double ot_dirphi;


  for(int ievent = 0; ievent < tree->GetEntries(); ++ievent) {
    tree->GetEntry(ievent);

    const double tx1 = true_xs->size() == 2 ? true_xs->at(0) : 0.;
    const double ty1 = true_ys->size() == 2 ? true_ys->at(0) : 0.;
    const double tz1 = true_zs->size() == 2 ? true_zs->at(0) : 0.;
    const double te1 = true_zs->size() == 2 ? true_enes->at(0) : 0.;
    const double tx2 = true_xs->size() == 2 ? true_xs->at(1) : 0.;
    const double ty2 = true_ys->size() == 2 ? true_ys->at(1) : 0.;
    const double tz2 = true_zs->size() == 2 ? true_zs->at(1) : 0.;
    const double te2 = true_zs->size() == 2 ? true_enes->at(1) : 0.;
    const TVector3 tp1{tx1,ty1,tz1};
    const TVector3 tp2{tx2,ty2,tz2};
    if ( true_xs->size() > 2 ) {
      std::cout << true_xs->size() << "true hits!" << std::endl;
    }
    // calculate the angles (sort them according to z position)
    ot_truetheta = true_zs->size() == 2 ? (tz1 < tz2 ? (tp2-tp1).Unit().Theta() : (tp1-tp2).Unit().Theta()) : 0.;
    ot_truephi = true_zs->size() == 2 ? (tz1 < tz2 ? (tp2-tp1).Unit().Phi() : (tp1-tp2).Unit().Phi()) : 0.;
    //bkg->Fill(
    sig->Fill(ot_truetheta*TMath::RadToDeg(),ot_truephi*TMath::RadToDeg());

    // loop over all reconstructed spacepoints
    for(size_t isps = 0; isps < pxs->size(); ++isps) {
      const double adc1 = adcs->at(isps);
      if(adc1 < 100 || adc1 > 1000) continue;

      const double px1 = pxs->at(isps);
      const double py1 = pys->at(isps);
      const double pz1 = pzs->at(isps);
      const TVector3 p1{px1,py1,pz1};

      dist->Fill((p1-tp1).Mag());
      dist->Fill((p1-tp2).Mag());
      
      for(size_t jsps = isps+1; jsps < pxs->size(); ++jsps) {
	const double adc2 = adcs->at(jsps);
        if(adc2 < 100 || adc2 > 1000) continue;
	const double px2 = pxs->at(jsps);
        const double py2 = pys->at(jsps);
        const double pz2 = pzs->at(jsps);
	const TVector3 p2{px2,py2,pz2};
	if ( (p2-tp1).Mag() > 10 &&
	     (p2-tp2).Mag() > 10 &&
	     (p1-tp1).Mag() > 10
	     && (p1-tp2).Mag() > 10 )
	  bkg->Fill(ot_dirtheta*TMath::RadToDeg(),ot_dirphi*TMath::RadToDeg());
	const TVector3& diff = p2.Z() > p1.Z() ? p2 - p1 : p1 - p2;
        const TVector3& dir = diff.Unit();

	bkg_and_sig->Fill(ot_dirtheta*TMath::RadToDeg(),ot_dirphi*TMath::RadToDeg());
	
	ot_dirtheta = dir.Theta();
        ot_dirphi = dir.Phi();


	
	
      }
    }
    
  }

  //tree->Draw("log10(pl2_integs)>>bkg");
  //tree->Draw("log10(pl2_true_integs)>>sig","log10(pl2_true_integs)!=0","same");

  auto c1 = new TCanvas();
  //c1->SetLogy();

  sig->Draw("colz");
  TLegend *myleg = new TLegend();
  myleg = c1->BuildLegend(0.1,0.7,0.3,0.9,"","");
  
  c1->SaveAs("Angles_signal.pdf");

  bkg->Draw("colz");
  c1->SaveAs("Angles_background.pdf");

  bkg_and_sig->Draw("colz");
  c1->SaveAs("Angles_background_and_signal.pdf");

  dist->Draw("");
  c1->SaveAs("Angles_distance_to_true.pdf");
  
}
