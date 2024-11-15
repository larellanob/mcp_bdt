void blips_make_pairs(TString g_directory, int g_run, TString g_sample_type, bool g_validation = false, int g_sig_mass = 100)
{
  /** makes blip pair ntuples from BlipAna files

      g_directory and g_run are the working directory and run (1, 2,
      3, 4, 5).

      a working directory structure will be used for the analysis
      framework, so BlipAna files must be in ${g_directory}/BlipAna


      and example fo the structure is (separated by run ${g_run})

      ${g_directory}/${g_run}/
      - pot/
        - tpot_hist_BlipAna_numi_${g_run}_data.txt
	- tpot_hist_BlipAna_numi_${g_run}_ext.txt
	- tpot_hist_BlipAna_numi_${g_run}_fhc_nu_overlay.txt
      - BlipAna/
        - hist_BlipAna_numi_${g_run}_data.root
        - hist_BlipAna_numi_${g_run}_ext.root
        - hist_BlipAna_numi_${g_run}_fhc_signal_30mev.root
        - hist_BlipAna_numi_${g_run}_fhc_nu_overlay.root
      - blippairs/
        - ntuple_blippairs_numi_${g_run}_data.root
        - ntuple_blippairs_numi_${g_run}_ext.root
        - ntuple_blippairs_numi_${g_run}_fhc_signal_30mev.root
        - ntuple_blippairs_numi_${g_run}_fhc_nu_overlay.root
	
      g_sample_type must be (data, ext, fhc_nu_overlay, sig)

      blip pair files can be used either for training, or for
      validation (g_validation flag)

      if you use g_sample_type = sig, you can also specify the mass
      with g_sig_mass. default is 100 (mev)

      if g_sig_mass = 0 it will use "combinedmasses", which is the
      combined masses root file used for training
   */

  
  TString s_validation = "training";
  if ( g_validation ) {
    s_validation = "validation";
  }
  const char * ch_validation = s_validation.Data();
  TString work_dir;
  work_dir = Form("%s/run%i/",
		  g_directory.Data(),
		  g_run
		  );

  TString filename;

  if ( g_sample_type == "sig" ) {
    g_sample_type = Form("%s_%imev",g_sample_type.Data(),g_sig_mass);
  }

  if ( !g_sample_type.Contains("sig") ) {
    filename = Form("%s/BlipAna/hist_BlipAna_numi_run%i_%s.root",
		    work_dir.Data(),
		    g_run,
		    g_sample_type.Data()
		    );
  } else { // if signal
    std::cout << "signal" << std::endl;
    if ( g_sig_mass != 0 ) { // if not combined masses
      std::cout << "signal not combined" << std::endl;

      filename = Form("%s/BlipAna/sig_%s/hist_BlipAna_numi_run%i_%s.root",
		      work_dir.Data(),
		      ch_validation,
		      g_run,
		      g_sample_type.Data());
    } else { // if combined masses
      std::cout << "signal combined" << std::endl;
      filename = Form("%s/BlipAna/sig_%s/hist_BlipAna_numi_run%i_sig_combinedmasses.root",
		      work_dir.Data(),
		      ch_validation,
		      g_run);
    }
  }
  std::cout << Form("Working on %s sample: %s\n",s_validation.Data(),filename.Data());
  TFile *f = new TFile(filename.Data());
  TTree *anatree;
  TTree *pottree;
  if ( g_validation) {
    anatree = (TTree *)f->Get("blipanapot/anatree");
    pottree = (TTree *)f->Get("blipanapot/pottree");
  } else {
    anatree = (TTree *)f->Get("blipana/anatree");
    pottree = (TTree *)f->Get("blipana/pottree");
  }
  TTreeReader reader(anatree);

  TTreeReaderValue<int> evt(reader,"event");
  TTreeReaderValue<int> run(reader,"run");
  TTreeReaderValue<int> subrun(reader,"subrun");
  
  TTreeReaderArray<float> x(reader,"blip_x");
  TTreeReaderArray<float> y(reader,"blip_y");
  TTreeReaderArray<float> z(reader,"blip_z");
  TTreeReaderArray<float> dx(reader,"blip_dx");
  TTreeReaderArray<float> dyz(reader,"blip_dyz");
  TTreeReaderArray<float> size(reader,"blip_size");
  TTreeReaderArray<int> charge(reader,"blip_charge");
  TTreeReaderArray<float> energy(reader,"blip_energy");
  TTreeReaderArray<float> yzcorr(reader,"blip_yzcorr");
  TTreeReaderArray<bool> incylinder(reader,"blip_incylinder");
  TTreeReaderArray<float> proxtrkdist(reader,"blip_proxtrkdist");
  TTreeReaderArray<int> edepid(reader,"blip_edepid");
  // cluster vars
  TTreeReaderArray<int> clust0(reader,"blip_pl0_clustid");
  TTreeReaderArray<int> clust1(reader,"blip_pl1_clustid");
  TTreeReaderArray<int> clust2(reader,"blip_pl2_clustid");
  TTreeReaderArray<int> clust_plane(reader,"clust_plane");
  TTreeReaderArray<int> clust_nhits(reader,"clust_nhits");
  TTreeReaderArray<float> clust_amp(reader,"clust_amp");

  
  
  int nevents = anatree->GetEntries();
  int event = -1;

  // for the zoomed in region
  double th_low = 27.0;
  double th_hig = 28.5;
  double ph_low = 3.0;
  double ph_hig = 8.0;
  int th_bins = 30;
  int ph_bins = 20;
  
  
  

  TH2F * h2_angle_sig = new TH2F("h2_angle_sig","Signal blip pairs;#theta;#phi",100,0,90,100,-180,180);
  TH2F * h2_angle_bkg = new TH2F("h2_angle_bkg","Background blip-pairs;#theta;#phi",100,0,90,100,-180,180);
  TH2F * h2_angle_mix = new TH2F("h2_angle_mix","Mixed blip-pairs;#theta;#phi",100,0,90,100,-180,180);
  TH2F * h2_angle_dat = new TH2F("h2_angle_dat","Data blip pairs;#theta;#phi",100,0,90,100,-180,180);
  TH2F * h2_angle_sig_zoom = new TH2F("h2_angle_sig_zoom",
				      "Signal blip pairs;#theta;#phi",
				      th_bins,th_low,th_hig,
				      ph_bins,ph_low,ph_hig);
  TH2F * h2_angle_bkg_zoom = new TH2F("h2_angle_bkg_zoom",
				      "Background blip-pairs;#theta;#phi",
				      th_bins,th_low,th_hig,
				      ph_bins,ph_low,ph_hig);
  TH2F * h2_angle_mix_zoom = new TH2F("h2_angle_mix_zoom",
				      "Mixed blip-pairs;#theta;#phi",
				      th_bins,th_low,th_hig,
				      ph_bins,ph_low,ph_hig);
  TH2F * h2_angle_dat_zoom = new TH2F("h2_angle_dat_zoom",
				      "Data blip-pairs;#theta;#phi",
				      th_bins,th_low,th_hig,
				      ph_bins,ph_low,ph_hig);

  // blip-pair variables that will go into the bdt training
  int bp_run;
  int bp_subrun;
  int bp_event;
  int bp_blipid1;
  int bp_blipid2;
  // // positions of blip1 and 2
  double bp_px1;
  double bp_py1;
  double bp_pz1;
  double bp_px2;
  double bp_py2;
  double bp_pz2;
  // // extent of blip1 and 2
  double bp_dx1;
  double bp_dyz1;
  double bp_dx2;
  double bp_dyz2;
  // // size and e of min and max
  double bp_emin;
  double bp_emax;
  double bp_sizemin;
  double bp_sizemax;
  // // min max hits in planes
  int bp_nhits_p2_max;
  int bp_nhits_p2_min;
  int bp_nhits_p01_max;
  int bp_nhits_p01_min;
  // // pair geometry
  double bp_dist;
  double bp_th;
  double bp_ph;
  // min and max dist to tracks
  double bp_proxtrkdist_min;
  double bp_proxtrkdist_max;
  // transverse distance to 3rd point
  double bp_perp_dist_closest_behind;
  double bp_perp_dist_closest_between;
  double bp_perp_dist_closest_ahead;
  // comparisons to third blips
  double bp_perp_dist_closest_overall;
  double bp_dist_closest_overall;
  
  // output file
  TString filename_out = "";
  gSystem->Exec(Form("mkdir -p %s/blippairs/sig_validation",
		     work_dir.Data()));
  gSystem->Exec(Form("mkdir -p %s/blippairs/sig_training",
		     work_dir.Data()));
  if ( g_sample_type != "sig" ) {
    filename_out = Form("%s/blippairs/hist_blippairs_numi_run%i_%s.root",
			work_dir.Data(),
			g_run,
			g_sample_type.Data()
			);
  } else { // if signal
    if ( g_sig_mass != 0 ) { // if not combined masses
      filename_out = Form("%s/blippairs/sig_%s/hist_blippairs_numi_run%i_%s.root",
			  work_dir.Data(),
			  ch_validation,
			  g_run,
			  g_sample_type.Data());
    } else { // if combined masses
      filename_out = Form("%s/blippairs/sig_%s/hist_blippairs_numi_run%i_sig_combinedmasses.root",
			  work_dir.Data(),
			  ch_validation,
			  g_run);
    }
  }
  TFile *fout = new TFile(filename_out.Data(),"recreate");

  // declare, fill, and write pottree right here, if you have it
  TTree *out_pottree = new TTree("total_pot","Total pot from sample");
  double tot_pot = 0;
  double pot = 0;
  // this is a short tree, so it's ok to use this inefficient looping
  if ( pottree != nullptr ) {
    pottree->SetBranchAddress("totpot",&pot);
    for(int ievent = 0; ievent < pottree->GetEntries(); ++ievent) {
      pottree->GetEntry(ievent);
      tot_pot += pot;
    }
  } else { 
    std::cout << "No POT Tree found\n";
  }
  std::cout << "tot_pot " << tot_pot << std::endl;
  auto b_pot = out_pottree->Branch("tot_pot",&tot_pot);
  out_pottree->Fill();
  out_pottree->Write();

  // now sig/bkg/mix trees  
  TTree *tsig = new TTree("sig_blip_pairs","Both blips signal");
  TTree *tbkg = new TTree("bkg_blip_pairs","Both blips background");
  TTree *tmix = new TTree("mix_blip_pairs","One blip signal, one blip background");
  TTree *tdat = new TTree("dat_blip_pairs","Data blip pairs");

  std::vector<TTree*> v_trees;
  bool is_data = false;
  if ( g_sample_type == "data" || g_sample_type == "ext" ) {
    is_data = true;
    v_trees.push_back(tdat);
  }  else {
    v_trees.push_back(tsig);
    v_trees.push_back(tbkg);
    v_trees.push_back(tmix);
  }


  for ( TTree * tout :  v_trees ) {
    // // non-training info
    tout->Branch("run",&bp_run);
    tout->Branch("subrun",&bp_subrun);
    tout->Branch("event",&bp_event);
    tout->Branch("bp_blipid1",&bp_blipid1);
    tout->Branch("bp_blipid2",&bp_blipid2);
    // // positions of blip1 and 2
    tout->Branch("bp_px1",&bp_px1);
    tout->Branch("bp_py1",&bp_py1);
    tout->Branch("bp_pz1",&bp_pz1);
    tout->Branch("bp_px2",&bp_px2);
    tout->Branch("bp_py2",&bp_py2);
    tout->Branch("bp_pz2",&bp_pz2);
    // // extent of blip1 and 2
    tout->Branch("bp_dx1",&bp_dx1);
    tout->Branch("bp_dyz1",&bp_dyz1);
    tout->Branch("bp_dx2",&bp_dx2);
    tout->Branch("bp_dyz2",&bp_dyz2);
    // // size and e of min and max
    tout->Branch("bp_emin",&bp_emin);
    tout->Branch("bp_emax",&bp_emax);
    tout->Branch("bp_sizemin",&bp_sizemin);
    tout->Branch("bp_sizemax",&bp_sizemax);
    // // min max hits in planes
    tout->Branch("bp_nhits_p2_max",&bp_nhits_p2_max);
    tout->Branch("bp_nhits_p2_min",&bp_nhits_p2_min);
    tout->Branch("bp_nhits_p01_max",&bp_nhits_p01_max);
    tout->Branch("bp_nhits_p01_min",&bp_nhits_p01_min);
    // // pair geometry
    tout->Branch("bp_dist",&bp_dist);
    tout->Branch("bp_th",&bp_th);
    tout->Branch("bp_ph",&bp_ph);
    // min and max dist to tracks
    tout->Branch("bp_proxtrkdist_min",&bp_proxtrkdist_min);
    tout->Branch("bp_proxtrkdist_max",&bp_proxtrkdist_max);
    // transverse distance to 3rd point
    tout->Branch("bp_perp_dist_closest_behind",&bp_perp_dist_closest_behind);
    tout->Branch("bp_perp_dist_closest_between",&bp_perp_dist_closest_between);
    tout->Branch("bp_perp_dist_closest_ahead",&bp_perp_dist_closest_ahead);
    // comparisons to third blips
    tout->Branch("bp_perp_dist_closest_overall",&bp_perp_dist_closest_overall);
    tout->Branch("bp_dist_closest_overall",&bp_dist_closest_overall);
  }

  std::cout << "Will loop over events\n";  
  
  //gStyle->SetOptStat(0);
  while ( reader.Next() ) {
    event++;
    bp_event = *evt;
    bp_run = *run;
    bp_subrun = *subrun;
    //if  ( event > 200 ) break;

    // fiducial cut lambda
    auto cut_fiducial = [&x,&y,&z](int i) {
      return
	(x[i] > 0.0 && x[i] < 260.0)
	&& (y[i] > -90.0 && y[i] < 90.0)
	&& (z[i] > 40.0 && z[i] < 1000.0);
    };

    for ( int i = 0; i < x.GetSize(); i++ ) {
      // fiducial cut
      if ( !cut_fiducial(i) ) {
	continue;
      }

      // blip pair classification
      for ( int j = i+1; j < x.GetSize(); j++ ) {

	// fiducial cut	
	if ( !cut_fiducial(j) ) {
	  continue;
	}
	bool is_sig = false;
	bool is_bkg = false;
	bool is_mix = false;

	// reset blip pair variables
	bp_blipid1 = -99;
	bp_blipid2 = -99;
	bp_px1 = -99;
	bp_py1 = -99;
	bp_pz1 = -99;
	bp_px2 = -99;
	bp_py2 = -99;
	bp_pz2 = -99;
	bp_dx1 = -99;
	bp_dyz1 = -99;
	bp_dx2 = -99;
	bp_dyz2 = -99;
	bp_emin = -99;
	bp_emax = -99;
	bp_sizemin = -99;
	bp_sizemax = -99;
	bp_nhits_p2_max = -99;
	bp_nhits_p2_min = -99;
	bp_nhits_p01_max = -99;
	bp_nhits_p01_min = -99;
	bp_dist = -99;
	bp_th = -99;
	bp_ph = -99;
	bp_proxtrkdist_min = -99;
	bp_proxtrkdist_max = -99;
	bp_perp_dist_closest_behind = -99;
	bp_perp_dist_closest_between = -99;
	bp_perp_dist_closest_ahead = -99;
	bp_perp_dist_closest_overall = -99;
	bp_dist_closest_overall = -99;
	
	TVector3 v_b1 (x[i],y[i],z[i]);
	TVector3 v_b2 (x[j],y[j],z[j]);
	// force v_b2.Z > v_b1.Z by flipping 1 and 2 if necessary
	// need to have new indices
	int ip = i;
	int jp = j;
	if ( z[j] < z[i] ) {
	  // flip the vectors
	  TVector3 v_temp;
	  v_temp = v_b1;
	  v_b1 = v_b2;
	  v_b2 = v_temp;
	  // flipped indices
	  ip = j;
	  jp = i;
	} 
	TVector3 v_bp12 = v_b2-v_b1; // always positive z
	TVector3 v_bp_midpoint = 0.5*(v_b1+v_b2);

	double th = v_bp12.Theta()*TMath::RadToDeg();
	double ph = v_bp12.Phi()*TMath::RadToDeg();

	bool cut_angle = th >= th_low && th <= th_hig && ph >= ph_low && ph <= ph_hig;
	// if cut at this point or earlier, is_sig, is_mix and is_bkg
	// will all be false, and we don't fill any ttrees

	// actually the cut happens after these plots

	// blip-pair sig bkg mix
	if ( !is_data ) {
	  if ( edepid[i] != -9 && edepid[j] != -9 ) {
	    is_sig = true;
	  } else if ( edepid[i] != -9 || edepid[j] != -9 ) {
	    is_mix = true;
	  } else {
	    is_bkg = true;
	  }
	}
	
	// angle cut plots
	if ( is_data ) {
	  h2_angle_dat->Fill(th,ph);
	  if ( cut_angle ) {
	    h2_angle_dat_zoom->Fill(th,ph);
	  }
	} else {
	  if ( is_sig ) {
	    h2_angle_sig->Fill(th,ph);
	    if ( cut_angle ) {
	      h2_angle_sig_zoom->Fill(th,ph);
	    }
	  } else if ( is_mix ) {
	    h2_angle_mix->Fill(th,ph);
	    if ( cut_angle ) {
	      h2_angle_mix_zoom->Fill(th,ph);
	    }
	  } else if ( is_bkg ) {
	    h2_angle_bkg->Fill(th,ph);
	    if ( cut_angle ) {
	      h2_angle_bkg_zoom->Fill(th,ph);
	    }
	  } else {
	    std::cout << Form("WARNING: Event %i pair (%i,%i) is neither sig, bkg, mix, nor data\n",event,i,j);
	  }
	}

	// here is the cut, we set flags to false again
	if ( !cut_angle ) {
	  is_sig = false;
	  is_mix = false;
	  is_bkg = false;
	  continue;
	}

	// blip-pair tree
	// v_b1         : pos of blip1
	// v_b2         : pos of blip2
	// v_bp12       : vector joining blip-pair (positive z)
	// v_bp_midpoint: mid point position
	double dist = v_bp12.Mag(); // distance between blip-pair

	bp_blipid1 = ip; // ip: index of b1 (after sorting)
	bp_blipid2 = jp; // jp: index of b2 (after sorting)
		
	// blip pos
	bp_px1 = v_b1.X();
	bp_py1 = v_b1.Y();
	bp_pz1 = v_b1.Z();
	bp_px2 = v_b2.X();
	bp_py2 = v_b2.Y();
	bp_pz2 = v_b2.Z();
	bp_dx1 = dx[ip];
	bp_dyz1 = dyz[ip];
	bp_dx2 = dx[jp];
	bp_dyz2 = dyz[jp];
	// sort by energy
	if ( energy[ip] < energy[jp] ) {
	  bp_emin = energy[ip];
	  bp_emax = energy[jp];
	} else {
	  bp_emin = energy[jp];
	  bp_emax = energy[ip];
	}
	// sort by size
	if ( size[ip] < size[jp] ) {
	  bp_sizemin = size[ip];
	  bp_sizemax = size[jp];
	} else {
	  bp_sizemin = size[jp];
	  bp_sizemax = size[ip];
	}
	// clusters of blip1 and blip2
	// use these as indices when accessing cluster variables
	double c0_b1 = clust0[ip];
	double c1_b1 = clust1[ip];
	double c2_b1 = clust2[ip];
	double c0_b2 = clust0[jp];
	double c1_b2 = clust1[jp];
	double c2_b2 = clust2[jp];
	// min max hits in planes p2 (collection plane)
	if ( clust_nhits[c2_b1] < clust_nhits[c2_b2] ) {
	  bp_nhits_p2_min = clust_nhits[c2_b1];
	  bp_nhits_p2_max = clust_nhits[c2_b2];
	} else {
	  bp_nhits_p2_min = clust_nhits[c2_b2];
	  bp_nhits_p2_max = clust_nhits[c2_b1];
	}
	// min max hits in planes p0 p1 (induction planes)
	/**
	   Blips are constructed by matching clusters in the
	   collection plane to clusters on one, or both induction
	   planes (so called 2D or 3D blips, respectively). So one of
	   the indices c0_b1 c1_b1 could be the default value of -9,
	   in which case there is no cluster match in that induction
	   plane.

	   We will use the sum of the hits in both induction planes
	   (p01), so first we check if there's a cluster, then we sum,
	   otherwise ignore the plane that doesn't have a cluster. We
	   assume one of them will be matched, otherwise it wouldn't
	   be a blip.
	*/
	
	double sum_p01_nhits_b1 =
	  ( c0_b1 != -9 && c1_b1 != -9)           // if both planes match
	  ? clust_nhits[c0_b1]+clust_nhits[c1_b1] // then sum
	  : (( c0_b1 != -9 ) ? clust_nhits[c0_b1] // else if p0 matches use that
	     : clust_nhits[c1_b1]);               // else use p1
	double sum_p01_nhits_b2 =
	  ( c0_b2 != -9 && c1_b2 != -9)           // if both planes match
	  ? clust_nhits[c0_b2]+clust_nhits[c1_b2] // then sum
	  : (( c0_b2 != -9 ) ? clust_nhits[c0_b2] // else if p0 matches use that
	     : clust_nhits[c1_b2]);               // else use p1

	if ( sum_p01_nhits_b1 < sum_p01_nhits_b2 ) {
	  bp_nhits_p01_min = sum_p01_nhits_b1;
	  bp_nhits_p01_max = sum_p01_nhits_b2;
	} else {
	  bp_nhits_p01_min = sum_p01_nhits_b2;
	  bp_nhits_p01_max = sum_p01_nhits_b1;
	}

	// pair geometry
	bp_dist = dist;
	bp_th = th; // already in degrees
	bp_ph = ph; // already in degrees

	// min and max dist to track
	if ( proxtrkdist[ip] <  proxtrkdist[jp] ) {
	  bp_proxtrkdist_min = proxtrkdist[ip];
	  bp_proxtrkdist_max = proxtrkdist[jp];
	} else {
	  bp_proxtrkdist_min = proxtrkdist[jp];
	  bp_proxtrkdist_max = proxtrkdist[ip];
	}

	// transverse distance to 3rd point
	double perp_dist_closest_behind  = 1e6;
	double perp_dist_closest_between     = 1e6;
	double perp_dist_closest_ahead   = 1e6;
	double perp_dist_closest_overall = 1e6;
	double dist_closest_overall = 1e6;
	// we loop over 3rd points
	for ( int k = 0; k < x.GetSize(); k++ ) {
	  if ( k == i || k == j ) continue;
	  /**
	     We take the direction of our blip pair as axis and
	     separate those blips "behind" the first blip b1, those
	     "ahead" of blip b2, and those in between b1 and b2
	     ("between").

	     When considering a blip b3, its longitudinal component
	     w.r.t. the blip-pair direction is given by the dot
	     product, and transverse by Perp().
	  */
	  TVector3 v_b3 (x[k],y[k],z[k]);
	  TVector3 v_mid_to_b3 = v_b3-v_bp_midpoint;
	  double dist_mid_to_b3 = v_mid_to_b3.Mag();
	  // reminder: v_bp12 is the direction from b1 to b2
	  // we use unit vectors for it to not change magnitdues
	  double long_comp = v_mid_to_b3.Dot(v_bp12.Unit());
	  double perp_comp = v_mid_to_b3.Perp(v_bp12.Unit());
	  // compare longitudinal distance of b3 to the half distance
	  // connecting b1 and b2 "dist"/2
	  if ( std::abs(long_comp) < dist/2. ) {
	    // between
	    if ( perp_comp < perp_dist_closest_between ) {
	      perp_dist_closest_between = perp_comp;
	    }
	  } else { // if not in between
	    if ( long_comp < 0 ) {
	      // behind
	      if ( perp_comp < perp_dist_closest_behind ) {
		perp_dist_closest_behind = perp_comp;
	      }
	    } else {
	      // ahead
	      if ( perp_comp < perp_dist_closest_ahead ) {
		perp_dist_closest_ahead = perp_comp;
	      }
	    }
	  }
	  if ( perp_comp < perp_dist_closest_overall ) {
	    perp_dist_closest_overall = perp_comp;
	  }
	  if ( dist_mid_to_b3 < dist_closest_overall ) {
	    dist_closest_overall = dist_mid_to_b3;
	  }
	  // don't push_back yet, we need to loop over all b3s.
	} // end looping over blip k

	// these values can still be unchanged from default if nothing
	// behind, between or ahead, which ruins scale of my plots. we
	// set it to a closer value that is still clearly anomalous
	if ( perp_dist_closest_behind == 1e6 ) {
	  perp_dist_closest_behind = -4.0;
	}
	if ( perp_dist_closest_between == 1e6 ) {
	  perp_dist_closest_between = -4.0;
	}
	if ( perp_dist_closest_ahead == 1e6 ) {
	  perp_dist_closest_ahead = -4.0;
	}
	// these ones are extremely unlikely but not impossible, there
	// would need to be literally two blips in the entire event
	if ( perp_dist_closest_overall == 1e6 ) {
	  perp_dist_closest_overall = -4.0;
	}
	if ( dist_closest_overall == 1e6 ) {
	  dist_closest_overall = -4.0;
	}
	
	bp_perp_dist_closest_behind = perp_dist_closest_behind;
	bp_perp_dist_closest_between = perp_dist_closest_between;
	bp_perp_dist_closest_ahead = perp_dist_closest_ahead;
	bp_perp_dist_closest_overall = perp_dist_closest_overall;
	bp_dist_closest_overall = dist_closest_overall;
	
	// fill the trees (if not empty?)
	// if it's "empty" it will be neither data, sig, mix or bkg
	if ( is_data ) {
	  tdat->Fill();
	} else {
	  if ( is_sig ) {
	    tsig->Fill();
	  }
	  if ( is_mix ) {
	    tmix->Fill();
	  }
	  if ( is_bkg ) {
	    tbkg->Fill();
	  }
	}
	
      } // end looping over blip j
      
    } // end looping over blip i
  } // end looping over file

  if ( is_data ) {
    tdat->Write();
  } else {
    tsig->Write();
    tmix->Write();
    tbkg->Write();
  }
  
  TCanvas *c1 = new TCanvas();
  TString img_dir = Form("%s/run%i/blippairs/img/%s",
			 g_directory.Data(),
			 g_run,
			 g_sample_type.Data());
  if ( g_sample_type.Contains("sig") ) {
    img_dir = Form("%s/run%i/blippairs/img/sig_%s",
		   g_directory.Data(),
		   g_run,
		   ch_validation);
  }
  gSystem->Exec(Form("mkdir -p %s/",img_dir.Data()));

  c1->Print(Form("%s/pairs_ang_%s.pdf[",img_dir.Data(),g_sample_type.Data()));
  h2_angle_sig->Draw("colz");
  c1->Print(Form("%s/pairs_ang_%s.pdf",img_dir.Data(),g_sample_type.Data()));
  h2_angle_bkg->Draw("colz");
  c1->Print(Form("%s/pairs_ang_%s.pdf",img_dir.Data(),g_sample_type.Data()));
  h2_angle_mix->Draw("colz");
  c1->Print(Form("%s/pairs_ang_%s.pdf",img_dir.Data(),g_sample_type.Data()));
  h2_angle_sig_zoom->Draw("colz");
  c1->Print(Form("%s/pairs_ang_%s.pdf",img_dir.Data(),g_sample_type.Data()));
  h2_angle_bkg_zoom->Draw("colz");
  c1->Print(Form("%s/pairs_ang_%s.pdf",img_dir.Data(),g_sample_type.Data()));
  h2_angle_mix_zoom->Draw("colz");
  c1->Print(Form("%s/pairs_ang_%s.pdf",img_dir.Data(),g_sample_type.Data()));
  c1->Print(Form("%s/pairs_ang_%s.pdf]",img_dir.Data(),g_sample_type.Data()));


  if ( is_data ) {
    std::cout << "| data |" << std::endl;
    std::cout << Form("| %lli  |\n",
		      tdat->GetEntries());
  } else {
    std::cout << "| mass |  sig |  bkg | mix |" << std::endl;
    std::cout << Form("|  %i | %lli | %lli | %lli |\n",
		      g_sig_mass,
		      tsig->GetEntries(), 
		      tbkg->GetEntries(), 
		      tmix->GetEntries());
  }
  
}
