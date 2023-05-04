bool compare_pair(std::pair<int,int> p1, std::pair<int,int> p2)
{
  if ( p1.first == p2.first ) {
    return ( p1.second < p2.second );
  } else {
    return ( p1.first < p2.first );
  }
}

void match_detvar_events(int mass = 100, TString tag = "def-th_1000kev_2hits")
{
  std::vector<TString> detvar_v
    {
      "wiremod_x",
      "wiremod_yz",
      "wiremod_anglexz",
      "wiremod_angleyz",
      "wiremod_dEdx",
      "ly_down",
      "ly_rayleigh",
      "ly_atten",
      "sce",
      "recomb"
    };
  
  
  // open cv file
  TString cv_filename = Form("systematics/detvar_%s_%imev_v08_00_00_61/spacepoints_CV_raw.root",
			     tag.Data(),mass);
  TFile *f_cv = new TFile(cv_filename);

  TTreeReader reader("ana/t", f_cv);
  TTreeReaderValue<Int_t> run(reader, "run");
  TTreeReaderValue<Int_t> eve(reader, "evt");

  // gather list of events
  std::vector<std::pair<int,int>> cv_evt;
  while ( reader.Next() ) {
    std::pair<int,int> evt_run { *run,*eve };
    cv_evt.push_back(evt_run);
  }
  // remove CV duplicates, if any  
  sort(cv_evt.begin(), cv_evt.end(), compare_pair);
  auto duplicate = std::adjacent_find(cv_evt.begin(), cv_evt.end());
  while ( duplicate != cv_evt.end() ) {
    cv_evt.erase(duplicate);
    duplicate = std::adjacent_find(cv_evt.begin(), cv_evt.end());
  }
  if ( duplicate != cv_evt.end() ) {
    // should have been removed by this point
    std::cout << "duplicate found at: " << duplicate->first << " " << duplicate->second << std::endl;
  } else {
    std::cout << "CV: no duplicates" << std::endl;
  }
  
  // go through detvars and alter list
  for ( auto detvar: detvar_v ) {
    std::cout << "-----------------" << std::endl;
    std::cout << "doing detvar: " << detvar << std::endl;
    TString detvar_filename = Form("systematics/detvar_%s_%imev_v08_00_00_61/spacepoints_%s_raw.root",
				   tag.Data(),mass, detvar.Data());
    TFile *f_detvar = new TFile(detvar_filename);
    TTreeReader reader_dv("ana/t", f_detvar);
    TTreeReaderValue<Int_t> run_dv(reader_dv, "run");
    TTreeReaderValue<Int_t> eve_dv(reader_dv, "evt");

    std::vector<std::pair<int,int>> dv_evt;

    // store detvar vector
    while ( reader_dv.Next() ) {
      std::pair<int,int> evt_run_dv {*run_dv,*eve_dv};
      dv_evt.push_back(evt_run_dv);
      if ( evt_run_dv.first == 6502 && evt_run_dv.second == 2702 && detvar == "wiremod_x" ) {
	std::cout << "HERE" << std::endl;
      }

    }
    // sort detvar vector
    sort(dv_evt.begin(), dv_evt.end(), compare_pair);
    // remove detvar duplicates
    auto duplicate = std::adjacent_find(dv_evt.begin(), dv_evt.end());
    int duplicates_removed = 0;
    while ( duplicate != dv_evt.end() ) {
      //std::cout << "duplicate found at: " << duplicate->first << " " << duplicate->second << std::endl;
      dv_evt.erase(duplicate);
      duplicate = std::adjacent_find(dv_evt.begin(), dv_evt.end());
      duplicates_removed++;
    }
    std::cout << Form("removed %i duplicates",duplicates_removed) << std::endl;
    // remove from CV events not in detvar
    std::vector<std::pair<int,int>> cv_evt_copy = cv_evt;
    int found_counter = 0;
    int deleted = 0;
    for ( auto p: cv_evt_copy ) {
      auto it = std::find(dv_evt.begin(),dv_evt.end(),p);
      if ( it != dv_evt.end() ) {
	found_counter++;
      } else {
	auto it2 = std::find(cv_evt.begin(),cv_evt.end(),p);
	cv_evt.erase(it2);
	deleted++;
      }
    }
    std::cout << "Deleted: " << deleted << std::endl;
    ofstream dvfile;
    dvfile.open(Form("runevt_%s.txt",detvar.Data()));
    for ( auto e: dv_evt ) {
      dvfile << e.first << " " << e.second << std::endl;
    }
    dvfile.close();
    std::cout << "CV     vector size: " << cv_evt.size() << std::endl;
    std::cout << "detvar vector size: " << dv_evt.size() << std::endl;
    std::cout << "found CV in detvar: " << found_counter << std::endl;
    
  }
  ofstream cvfile;
  cvfile.open(Form("runevt_CV.txt"));
  for ( auto e: cv_evt ) {
    cvfile << e.first << " " << e.second << std::endl;
  }
  cvfile.close();

  
  // v find
  // v erase
  
  std::cout << cv_evt.size() << std::endl;

  // now filter spacepoint files
  // for CV then detvar files
  
  TString cvout_filename = Form("systematics/detvar_%s_%imev_v08_00_00_61/spacepoints_CV.root",
				tag.Data(),mass);
  TFile *f_cvout = new TFile(cvout_filename,"recreate");
  f_cvout->mkdir("ana");
  f_cvout->cd("ana");
  TTree * oldtree;
  f_cv->GetObject("ana/t",oldtree);
  TTree * oldtree_pot;
  f_cv->GetObject("ana/pottree",oldtree_pot);
  auto newtree_pot = oldtree_pot->CloneTree();
  newtree_pot->Fill();

  auto newtree = oldtree->CloneTree(0);
  newtree->SetBranchStatus("*",1);
  oldtree->SetBranchStatus("*",1);

  int write_counter = 0;
  for ( int i = 0; i < oldtree->GetEntries(); i++ ) {
    oldtree->GetEntry(i);
    std::pair<int,int> evt_run {*run, *eve};
    auto it = std::find(cv_evt.begin(),cv_evt.end(),evt_run);
    if ( it != cv_evt.end() ) {
      write_counter++;
      newtree->Fill();
    }
  }
  f_cvout->Write();


  std::cout << Form("wrote %i events",write_counter) << std::endl;

  for ( auto detvar: detvar_v ) {
    std::vector<std::pair<int,int>> local_evt = cv_evt;
    std::cout << "--------------------" << std::endl;
    std::cout << "detvar : " << detvar << std::endl;
    TString detvar_filename = Form("systematics/detvar_%s_%imev_v08_00_00_61/spacepoints_%s_raw.root",
				   tag.Data(),mass, detvar.Data());
    TFile *f_detvar = new TFile(detvar_filename);
    TString dvout_filename = Form("systematics/detvar_%s_%imev_v08_00_00_61/spacepoints_%s.root",
				  tag.Data(),mass, detvar.Data());
    TFile *f_dvout = new TFile(dvout_filename,"recreate");
    f_dvout->mkdir("ana");
    f_dvout->cd("ana");
    TTree * oldtree_dv;
    f_detvar->GetObject("ana/t",oldtree_dv);
    TTree * oldtree_dv_pot;
    f_cv->GetObject("ana/pottree",oldtree_dv_pot);
    auto newtree_dv_pot = oldtree_dv_pot->CloneTree();
    newtree_dv_pot->Fill();
    
    auto newtree_dv = oldtree_dv->CloneTree(0);
    newtree_dv->SetBranchStatus("*",1);
    oldtree_dv->SetBranchStatus("*",1);
    
    int write_counter = 0;
    std::cout << "oldtree entries: " << oldtree_dv->GetEntries() << std::endl;
    int run_dv;
    int eve_dv;
    oldtree_dv->SetBranchAddress("run",&run_dv);
    oldtree_dv->SetBranchAddress("evt",&eve_dv);
    
    for ( int i = 0; i < oldtree_dv->GetEntries(); i++ ) {
      oldtree_dv->GetEntry(i);
      std::pair<int,int> evt_run {run_dv, eve_dv};
      auto it = std::find(local_evt.begin(),local_evt.end(),evt_run);
      if ( it != local_evt.end() ) {
	write_counter++;
	newtree_dv->Fill();
	local_evt.erase(it);
      }
    }
    std::cout << "local evt size (should be zero): " << local_evt.size() << std::endl;
    f_dvout->Write();
    std::cout << Form("wrote %i events",write_counter) << std::endl;

  }
  
  
}
