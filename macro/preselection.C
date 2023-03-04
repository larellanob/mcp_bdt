void preselection(const char* fn = "spacepoints_mc.root") {
  TFile *f = new TFile(fn);
  TTree *t = (TTree*)f->Get("ana/t");
  int run,  evt;
  t->SetBranchAddress("run",&run);
  t->SetBranchAddress("evt",&evt);
  vector<double> *pxs = 0, *pys = 0, *pzs = 0, *adcs = 0;
  vector<double> *true_xs = 0, *true_ys = 0, *true_zs = 0, *true_adcs = 0, *true_enes = 0;
  vector<int> *nh2 = 0, *nhO = 0, *nhM = 0;
  t->SetBranchAddress("sps_x",&pxs);
  t->SetBranchAddress("sps_y",&pys);
  t->SetBranchAddress("sps_z",&pzs);
  t->SetBranchAddress("pl2_integs",&adcs);
  t->SetBranchAddress("pl2_nhits",&nh2);
  t->SetBranchAddress("pl_othmax_nhits",&nhO);
  t->SetBranchAddress("pl2_othmin_nhits",&nhM);
  t->SetBranchAddress("elec_x",&true_xs);
  t->SetBranchAddress("elec_y",&true_ys);
  t->SetBranchAddress("elec_z",&true_zs);
  t->SetBranchAddress("elec_E",&true_enes);
  t->SetBranchAddress("pl2_true_integs",&true_adcs);

  TString output_filename {Form("%s_output.root",fn)};
  TString infile_compare = fn;
  if ( infile_compare != "spacepoints.root" ) {
    output_filename = infile_compare;
    output_filename.ReplaceAll("spacepoints","preselection");
  }
  std::cout << output_filename << std::endl;
  
  TFile *of = new TFile(output_filename,"recreate");

  // POT and normalization from spacepoints file
  TTree *pot_tree = (TTree*)f->Get("ana/pottree");
  TTree *o_pot_tree = new TTree("total_pot","total pot and events from sample");
  double tot_pot = 0;
  double pot = 0;
  double totevents = t->GetEntries();
  pot_tree->SetBranchAddress("totpot",&pot);
  for(int ievent = 0; ievent < pot_tree->GetEntries(); ++ievent) {
    pot_tree->GetEntry(ievent);
    tot_pot += pot;
  }
  std::cout << "tot_pot " << tot_pot << std::endl;
  auto b_pot = o_pot_tree->Branch("tot_pot",&tot_pot);
  auto b_evt = o_pot_tree->Branch("tot_evt",&totevents);
  o_pot_tree->Fill();
  o_pot_tree->Write();

  // events from spacepoints file
  TTree *ot1a = new TTree("tbg_train","background training");
  TTree *ot2a = new TTree("tsig_train","signal training");
  TTree *ot1b = new TTree("tbg_test","background");
  TTree *ot2b = new TTree("tsig_test","signal");

  double ot_px1;
  double ot_py1;
  double ot_pz1;
  double ot_adc1;
  double ot_dtp1;
  double ot_dtx1;
  double ot_dty1;
  double ot_dtz1;
  double ot_te1;
  size_t ot_it1;
  int ot_nhits1;
  int ot_nhitsO1;
  int ot_nhitsM1;

  double ot_px2;
  double ot_py2;
  double ot_pz2;
  double ot_adc2;
  double ot_dtp2;
  double ot_dtx2;
  double ot_dty2;
  double ot_dtz2;
  double ot_te2;
  size_t ot_it2;
  int ot_nhits2;
  int ot_nhitsO2;
  int ot_nhitsM2;

  double ot_dist;
  double ot_dirtheta;
  double ot_dirphi;

  double ot_truetheta;
  double ot_truephi;

  double ot_dist_closest_behind;
  double ot_dist_closest_mid;
  double ot_dist_closest_ahead;
  double ot_adc_closest_behind;
  double ot_adc_closest_mid;
  double ot_adc_closest_ahead;
  double ot_mindist_closest;
  double ot_mindist_adc_closest;
  int ot_nh2_closest_behind;
  int ot_nhO_closest_behind;
  int ot_nhM_closest_behind;
  int ot_nh2_closest_mid;
  int ot_nhO_closest_mid;
  int ot_nhM_closest_mid;
  int ot_nh2_closest_ahead;
  int ot_nhO_closest_ahead;
  int ot_nhM_closest_ahead;

  for(auto ot : { ot1a, ot1b, ot2a, ot2b }) {
  ot->Branch("run",&run);
  ot->Branch("evt",&evt);

  ot->Branch("px1",&ot_px1);
  ot->Branch("py1",&ot_py1);
  ot->Branch("pz1",&ot_pz1);
  ot->Branch("log10_adc1",&ot_adc1);
  ot->Branch("nhits1",&ot_nhits1);
  ot->Branch("nhitsO1",&ot_nhitsO1);
  ot->Branch("nhitsM1",&ot_nhitsM1);
  ot->Branch("dtp1",&ot_dtp1);
  ot->Branch("dtz1",&ot_dtx1);
  ot->Branch("dty1",&ot_dty1);
  ot->Branch("dtx1",&ot_dtz1);
  ot->Branch("ene1",&ot_te1);
  //ot->Branch("it1",&ot_it1);
  ot->Branch("px2",&ot_px2);
  ot->Branch("py2",&ot_py2);
  ot->Branch("pz2",&ot_pz2);
  ot->Branch("log10_adc2",&ot_adc2);
  ot->Branch("nhits2",&ot_nhits2);
  ot->Branch("nhitsO2",&ot_nhitsO2);
  ot->Branch("nhitsM2",&ot_nhitsM2);
  ot->Branch("dtp2",&ot_dtp2);
  ot->Branch("dtz2",&ot_dtx2);
  ot->Branch("dty2",&ot_dty2);
  ot->Branch("dtx2",&ot_dtz2);
  ot->Branch("ene2",&ot_te2);
  //ot->Branch("it2",&ot_it2);
  ot->Branch("dist",&ot_dist);
  ot->Branch("theta",&ot_dirtheta);
  ot->Branch("phi",&ot_dirphi);
  ot->Branch("true_theta",&ot_truetheta);
  ot->Branch("true_phi",&ot_truephi);

  ot->Branch("sqrt_dca_behind",&ot_dist_closest_behind);
  ot->Branch("sqrt_dca_mid",&ot_dist_closest_mid);
  ot->Branch("sqrt_dca_ahead",&ot_dist_closest_ahead);
  ot->Branch("dca_log10_adc_behind",&ot_adc_closest_behind);
  ot->Branch("dca_log10_adc_mid",&ot_adc_closest_mid);
  ot->Branch("dca_log10_adc_ahead",&ot_adc_closest_ahead);
  ot->Branch("dca_nh2_behind",&ot_nh2_closest_behind);
  ot->Branch("dca_nh2_mid",&ot_nh2_closest_mid);
  ot->Branch("dca_nh2_ahead",&ot_nh2_closest_ahead);
  ot->Branch("dca_nhO_behind",&ot_nhO_closest_behind);
  ot->Branch("dca_nhO_mid",&ot_nhO_closest_mid);
  ot->Branch("dca_nhO_ahead",&ot_nhO_closest_ahead);
  ot->Branch("dca_nhM_behind",&ot_nhM_closest_behind);
  ot->Branch("dca_nhM_mid",&ot_nhM_closest_mid);
  ot->Branch("dca_nhM_ahead",&ot_nhM_closest_ahead);
  ot->Branch("sqrt_mindist_closest",&ot_mindist_closest);
  ot->Branch("mindist_log10_adc_closest",&ot_mindist_adc_closest);
  }

  map<int,map<int,bool>> istrain_map_sig;
  map<int,map<int,bool>> istrain_map_bg;
  TRandom3 rand;
  auto is_train = [&rand](auto& istrain_map, const int run, const int evt, const double train_frac = 0.7) -> bool {
    auto r = istrain_map.find(run);
    if(r != istrain_map.end()) {
      auto const& rr = r->second;
      auto e = rr.find(evt);
      if(e != rr.end()) {
        return e->second;
      }
    }
    const bool train = (rand.Uniform() < train_frac);
    istrain_map[run][evt] = train;
    return train;
  };

  for(int ievent = 0; ievent < t->GetEntries(); ++ievent) {
    t->GetEntry(ievent);
    const bool has_truth = !true_xs->empty();
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
    ot_truetheta = true_zs->size() == 2 ? (tz1 < tz2 ? (tp2-tp1).Unit().Theta() : (tp1-tp2).Unit().Theta()) : 0.;
    ot_truephi = true_zs->size() == 2 ? (tz1 < tz2 ? (tp2-tp1).Unit().Phi() : (tp1-tp2).Unit().Phi()) : 0.;
    for(size_t isps = 0; isps < pxs->size(); ++isps) {
      const double adc1 = adcs->at(isps);
      if(adc1 < 100 || adc1 > 3630) continue;
      const double px1 = pxs->at(isps);
      const double py1 = pys->at(isps);
      const double pz1 = pzs->at(isps);
      ot_nhits1 = nh2->at(isps);
      ot_nhitsO1 = nhO->at(isps);
      ot_nhitsM1 = nhM->at(isps);
      const TVector3 p1{px1,py1,pz1};
      const bool is_true_hit1 = has_truth && true_adcs->at(isps) > 0.;
      const double dt1 = is_true_hit1 ? (p1-tp1).Mag() : -1.;
      const double dt2 = is_true_hit1 ? (p1-tp2).Mag() : -1.;
      const size_t it1 = is_true_hit1 ? (dt1 < dt2 ? 0 : 1) : -1;
      const double dtp1 = is_true_hit1 ? std::min(dt1,dt2) : -1.;
      const double dtx1 = is_true_hit1 ? (dt1 < dt2 ? px1 - tx1 : px1 - tx2) : -1.;
      const double dty1 = is_true_hit1 ? (dt1 < dt2 ? py1 - ty1 : py1 - ty2) : -1.;
      const double dtz1 = is_true_hit1 ? (dt1 < dt2 ? pz1 - tz1 : pz1 - tz2) : -1.;
      const double tene1 = is_true_hit1 ? (dt1 < dt2 ? te1 : te2) : -1.;

      ot_px1 = px1;
      ot_py1 = py1;
      ot_pz1 = pz1;
      ot_adc1 = log10(adc1);
      ot_dtx1 = dtx1;
      ot_dty1 = dty1;
      ot_dtz1 = dtz1;
      ot_dtp1 = dtp1;
      ot_te1 = tene1;
      ot_it1 = it1;

      for(size_t jsps = isps+1; jsps < pxs->size(); ++jsps) {
        const double adc2 = adcs->at(jsps);
        if(adc2 < 100 || adc2 > 3630) continue;
        const bool train_bg = rand.Uniform() < 0.7; //is_train(istrain_map_bg,run,evt,0.002);
        const bool train_sig = rand.Uniform() < 0.7; //is_train(istrain_map_sig,run,evt,0.7);
        const double px2 = pxs->at(jsps);
        const double py2 = pys->at(jsps);
        const double pz2 = pzs->at(jsps);
        ot_nhits2 = nh2->at(jsps);
        ot_nhitsO2 = nhO->at(jsps);
        ot_nhitsM2 = nhM->at(jsps);
        const TVector3 p2{px2,py2,pz2};
        const bool is_true_hit2 = has_truth && true_adcs->at(jsps) > 0.;
        const double dt1 = is_true_hit2 ? (p2-tp1).Mag() : -1.;
        const double dt2 = is_true_hit2 ? (p2-tp2).Mag() : -1.;
        const size_t it2 = is_true_hit2 ? (dt1 < dt2 ? 0 : 1) : -1;
        const double dtp2 = is_true_hit2 ? std::min(dt1,dt2) : -1.;
        const double dtx2 = is_true_hit2 ? (dt1 < dt2 ? px2 - tx1 : px2 - tx2) : -1.;
        const double dty2 = is_true_hit2 ? (dt1 < dt2 ? py2 - ty1 : py2 - ty2) : -1.;
        const double dtz2 = is_true_hit2 ? (dt1 < dt2 ? pz2 - tz1 : pz2 - tz2) : -1.;
        const double tene2 = is_true_hit2 ? (dt1 < dt2 ? te1 : te2) : -1.;

        ot_px2 = px2;
        ot_py2 = py2;
        ot_pz2 = pz2;
        ot_adc2 = log10(adc2);
        ot_dtx2 = dtx2;
        ot_dty2 = dty2;
        ot_dtz2 = dtz2;
        ot_dtp2 = dtp2;
        ot_te2 = tene2;
        ot_it2 = it2;

        const TVector3& diff = p2.Z() > p1.Z() ? p2 - p1 : p1 - p2;
        const TVector3& dir = diff.Unit();

        ot_dist = diff.Mag();
        ot_dirtheta = dir.Theta();
        ot_dirphi = dir.Phi();

        const double thetadeg = ot_dirtheta * TMath::RadToDeg();
        const double phideg = ot_dirphi * TMath::RadToDeg();

        if(thetadeg < 26 || thetadeg > 30 || phideg < 0 || phideg > 10) continue;


        ot_dist_closest_behind = 1e100, ot_dist_closest_mid = 1e100, ot_dist_closest_ahead = 1e100;
        ot_adc_closest_behind = 0., ot_adc_closest_mid = 0., ot_adc_closest_ahead = 0.;
        ot_nh2_closest_behind = 0, ot_nh2_closest_mid = 0, ot_nh2_closest_ahead = 0;
        ot_nhO_closest_behind = 0, ot_nhO_closest_mid = 0, ot_nhO_closest_ahead = 0;
        ot_nhM_closest_behind = 0, ot_nhM_closest_mid = 0, ot_nhM_closest_ahead = 0;
        ot_mindist_closest = 1e100;
        ot_mindist_adc_closest = 0.;
        for(size_t ksps = 0; ksps< pxs->size(); ++ksps) {
          if(ksps == isps || ksps == jsps) continue;
          const double adc3 = adcs->at(ksps);
          //if(adc3 < 100 || adc3 > 1000) continue;
          const double px3 = pxs->at(ksps);
          const double py3 = pys->at(ksps);
          const double pz3 = pzs->at(ksps);
          const int nhits2 = nh2->at(ksps);
          const int nhitsO = nhO->at(ksps);
          const int nhitsM = nhM->at(ksps);
          const TVector3 p3{px3,py3,pz3};
          const double dotp1 = (p3 - p1).Unit().Dot(dir);
          const double dotp2 = (p3 - p2).Unit().Dot(dir);
          const double dist = (p3-p1).Perp(dir);
	  const double sqrtdist = sqrt(dist);
          if(dotp1 < 0 && dotp2 < 0) {
            if(sqrtdist < ot_dist_closest_behind) {
              ot_dist_closest_behind = sqrtdist;
              ot_adc_closest_behind = log10(adc3);
	      //std::cout << adc3 << " " << log10(adc3) << std::endl;
              ot_nh2_closest_behind = nhits2;
              ot_nhO_closest_behind = nhitsO;
              ot_nhM_closest_behind = nhitsM;
            }
          }
          else if(dotp1 > 0 && dotp2 > 0) {
            if(sqrtdist < ot_dist_closest_ahead) {
              ot_dist_closest_ahead = sqrtdist;
              ot_adc_closest_ahead = log10(adc3);
              ot_nh2_closest_ahead = nhits2;
              ot_nhO_closest_ahead = nhitsO;
              ot_nhM_closest_ahead = nhitsM;
            }
          }
          else {
            if(sqrtdist < ot_dist_closest_mid) {
              ot_dist_closest_mid = sqrtdist;
              ot_adc_closest_mid = log10(adc3);
              ot_nh2_closest_mid = nhits2;
              ot_nhO_closest_mid = nhitsO;
              ot_nhM_closest_mid = nhitsM;
            }
          }
          if(sqrtdist < ot_mindist_closest) {
            ot_mindist_closest = sqrtdist;
            ot_mindist_adc_closest = log10(adc3);
          }

        }

        auto const& ot = (it1 != it2 && it1 < 2 && it2 < 2) ? (train_sig ? ot2a : ot2b) : (train_bg ? ot1a : ot1b);
        ot->Fill();
      }
    }
  }
  of->cd();
  for(auto ot : { ot1a, ot1b, ot2a, ot2b }) {
    ot->SetAlias("distinct_truth","it1+it2==1");
    ot->SetAlias("thetadeg","theta*180/pi");
    ot->SetAlias("phideg","phi*180/pi");
    ot->SetAlias("true_mean_thetadeg","2.76132e+01");
    ot->SetAlias("true_mean_sintheta_phideg","2.79794e+00");
    ot->SetAlias("sigma_thetadeg","4.28918e-01");
    ot->SetAlias("mean_thetadeg","2.17494e-01");
    ot->SetAlias("sigma_sintheta_phideg","8.89466e-01");
    ot->SetAlias("mean_sintheta_phideg","1.24913e-01");
    ot->SetAlias("thetapull","(thetadeg-true_mean_thetadeg-mean_thetadeg)/sigma_thetadeg");
    ot->SetAlias("phipull","(sin(theta)*phideg - true_mean_sintheta_phideg - mean_sintheta_phideg)/sigma_sintheta_phideg");
    ot->SetAlias("anglepull","sqrt(pow(thetapull,2)+pow(phipull,2))");
    //ot->SetAlias("mindist_closest","min(min(dca_behind,dca_mid),dca_ahead)");
    //ot->SetAlias("mindist_adc_closest","(mindist_closest==dca_behind)*dca_adc_behind+(mindist_closest==dca_mid)*dca_adc_mid+(mindist_closest==dca_ahead)*dca_adc_ahead");
    ot->Write();
  }
  delete of;
  delete f;
  std::cout << "Saved file " << output_filename << std::endl;
}
