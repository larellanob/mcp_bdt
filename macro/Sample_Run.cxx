class Sample
{
private:
  TString name = "";
  TString type = ""; // ext, genie, sig, data
  TFile *f_in;
  TFile *f_out;
  TTree *t_preselected;
  double pot = 0;
  double trig = 1;
  double scaling_factor = 1.0;
  TString weight_factor = "1";
  bool scaled = false; // to prevent or at least warn of scaling multiple times
  TString legend_title = "";
  Color_t legend_color = kBlack;

public:
  Sample(){};
  Sample(TString s_name,TString s_type, TString s_filename);
  ~Sample(){
    delete f_in;
  };
  TH1F * get_th1(TString);
  void set_bdt(TString s_filename) {
    f_out = new TFile(s_filename);
    std::cout << Form("Sample %s\n\toutput file: %s ",name.Data(),s_filename.Data()) << std::endl;
  }
  TTree * get_preselection();
  TString get_name() { return name; };
  double get_pot() { return pot; };
  void set_ext_triggers(double triggers) {
    trig = triggers;
  };
  void set_scaling_factor(double s) {
    scaling_factor = s;
    scaled = true;
  };
  void set_weight_factor(TString st) { weight_factor = st; };
  TString get_weight_factor() { return weight_factor; };
  double get_scaling_factor() { return scaling_factor;};
  bool previously_scaled() { return scaled; };
  TString get_type() { return type; };
  Sample operator+(Sample &s);
};

class Run
{
private:
  TString name = "Run0"; // name used for file output and such
  TString s_run = "0"; // run number (e.g. 1, 2, 3, etc.)
  TString s_horn = "XHC"; // horn current (FHC RHC or FULL (RHC+FHC))
  std::vector<Sample*> v_samples;
  double beamon_pot;
  double beamon_triggers;
  double beamoff_triggers;
public:
  Run(){};
  Run(TString s_name,
      TString run,
      TString horn
      ){
    name = s_name;
    s_run = run;
    s_horn = horn;
  };
  ~Run(){};
  void add_sample(Sample *);
  void set_pot(double pot, double trig_on, double trig_off) {
    beamon_pot = pot;
    beamon_triggers = trig_on;
    beamoff_triggers = trig_off;
  };
  void set_samples_scaling();
  double get_pot() { return beamon_pot; };
  TString get_name() { return name; };
  TString get_run() { return s_run; };
  TString get_horn() { return s_horn; };
  void print_samples_scaling();
  std::vector<Sample*> get_samples() {
    return v_samples;
  };
  friend class RunPlot;
  void add_signal_samples();
};



// implementation of member functions


  
Sample::Sample(TString s_name,TString s_type, TString s_filename) {
  name = s_name;
  type = s_type;
  f_in = new TFile(s_filename);
  std::cout << Form("Sample %s\n\tinput file: %s ",name.Data(),s_filename.Data()) << std::endl;

  // we set pot
  TTree *T_pot = (TTree*)f_in->Get("wcpselection/T_pot");
  TH1 *htemp = new TH1D("htemp","",1,0,1);
  T_pot->Draw("0.5>>htemp","pot_tor875","goff");
  const double mcpot = htemp->Integral();
  delete htemp;
  pot = mcpot;
  std::cout << "\tsample pot: " << pot << std::endl;
}

TTree * Sample::get_preselection()
{
  TTree *T_eval = (TTree*)f_in->Get("wcpselection/T_eval");
  T_eval->AddFriend((TTree*)f_in->Get("wcpselection/T_PFeval"));
  T_eval->AddFriend((TTree*)f_in->Get("wcpselection/T_KINEvars"));
  T_eval->AddFriend((TTree*)f_in->Get("wcpselection/T_BDTvars"));

  //TTree *tbdt = (TTree*)f_out->Get("tree");
  //std::cout << "tbdt at " << tbdt << std::endl;
  T_eval->AddFriend("tree",f_out);
  
  T_eval->SetAlias("flag_cosmic","match_completeness_energy<=0.1*truth_energyInside");
  T_eval->SetAlias("flag_sig","match_completeness_energy>0.1*truth_energyInside");
  T_eval->SetAlias("flag_kdar","truth_nuEnergy==235.531814575195312");
  //T_eval->SetAlias("flag_e_elas","intType==1098");
  T_eval->SetAlias("flag_presel", "match_found_asInt == 1 && stm_eventtype != 0 && stm_lowenergy ==0 && stm_LM ==0 && stm_TGM ==0 && stm_STM==0 && stm_FullDead == 0 && stm_clusterlength >0");
  
  //T_eval->SetAlias("flag_FC","match_isFC");

  T_eval->SetAlias("full_selection","flag_presel && kine_reco_Enu > 0.0 && match_isFC");
  return T_eval;
  
}

  

void Run::add_sample(Sample *new_sample)
{
  bool do_add = true;
  for ( Sample *s: v_samples ) {
    if ( new_sample->get_name() == s->get_name() ) {
      std::cout << "WARNING: sample with name " << s->get_name()
		<< " already in " << name
		<< ". It will not be added, to avoid duplicates." <<std::endl;
      do_add = false;
      break;
    }
  }
  if ( do_add ) {
    v_samples.push_back(new_sample);
  }
}


void Run::set_samples_scaling()
{
  for ( Sample *s: v_samples ) {
    if ( s->previously_scaled() ) {
      std::cout << "WARNING: you tried to scale sample " << s->get_name() <<
	", which had already been scaled before. Will skip it\n";
      continue;
    }
    if ( s->get_type() == "ext" ) {
      s->set_scaling_factor(beamon_triggers/beamoff_triggers);
    } else if ( s->get_type() == "genie" || s->get_type() == "sig" || s->get_type() == "dirt" ) {
      if ( s->get_pot() == 0 ) {
	std::cout << "WARNING: scaling sample with 0 pot" << std::endl;
      }
      s->set_scaling_factor(beamon_pot/s->get_pot());
    } else if ( s->get_type() == "data" ) {
      ; // data shouldn't be scaled 
    } else {
      std::cout << "WARNING: can't identify type of scaling to apply (sample is not ext, genie, sig, or data).\n";
    }
  }
}



void Run::print_samples_scaling()
{
  std::cout << "Printing scalings for " << name << std::endl;
  for ( Sample *s: v_samples ) {
    std::cout << Form("\tSample: %s\tScaling %.5f\n",s->get_name().Data(),s->get_scaling_factor());
  }
}


void Run::add_signal_samples()
{
  Sample *sig15mev = new Sample("signal_15mev",
				 "sig",
				 "/exp/uboone/data/users/arellano/thesis_checks/train_sample/combinedmasses/checkout_testing_20000kev_1hits_15mev.root");
  Sample *sig20mev = new Sample("signal_20mev",
				 "sig",
				 "/exp/uboone/data/users/arellano/thesis_checks/train_sample/combinedmasses/checkout_testing_20000kev_1hits_20mev.root");
  Sample *sig30mev = new Sample("signal_30mev",
				 "sig",
				 "/exp/uboone/data/users/arellano/thesis_checks/train_sample/combinedmasses/checkout_testing_20000kev_1hits_30mev.root");
  Sample *sig50mev = new Sample("signal_50mev",
				 "sig",
				 "/exp/uboone/data/users/arellano/thesis_checks/train_sample/combinedmasses/checkout_testing_20000kev_1hits_50mev.root");
  Sample *sig80mev = new Sample("signal_80mev",
				 "sig",
				 "/exp/uboone/data/users/arellano/thesis_checks/train_sample/combinedmasses/checkout_testing_20000kev_1hits_80mev.root");
  Sample *sig100mev = new Sample("signal_100mev",
				 "sig",
				 "/exp/uboone/data/users/arellano/thesis_checks/train_sample/combinedmasses/checkout_testing_20000kev_1hits_100mev.root");
  Sample *sig150mev = new Sample("signal_150mev",
				 "sig",
				 "/exp/uboone/data/users/arellano/thesis_checks/train_sample/combinedmasses/checkout_testing_20000kev_1hits_150mev.root");
  Sample *sig200mev = new Sample("signal_200mev",
				 "sig",
				 "/exp/uboone/data/users/arellano/thesis_checks/train_sample/combinedmasses/checkout_testing_20000kev_1hits_200mev.root");
  Sample *sig250mev = new Sample("signal_250mev",
				 "sig",
				 "/exp/uboone/data/users/arellano/thesis_checks/train_sample/combinedmasses/checkout_testing_20000kev_1hits_250mev.root");
  Sample *sig300mev = new Sample("signal_300mev",
				 "sig",
				 "/exp/uboone/data/users/arellano/thesis_checks/train_sample/combinedmasses/checkout_testing_20000kev_1hits_300mev.root");
  Sample *sig350mev = new Sample("signal_350mev",
				 "sig",
				 "/exp/uboone/data/users/arellano/thesis_checks/train_sample/combinedmasses/checkout_testing_20000kev_1hits_350mev.root");
  Sample *sig400mev = new Sample("signal_400mev",
				 "sig",
				 "/exp/uboone/data/users/arellano/thesis_checks/train_sample/combinedmasses/checkout_testing_20000kev_1hits_400mev.root");
  
  sig15mev->set_bdt("outputs/20000kev_1hits_combinedmasses/gpu_checkout_testing_20000kev_1hits_15mev.root");
  sig20mev->set_bdt("outputs/20000kev_1hits_combinedmasses/gpu_checkout_testing_20000kev_1hits_20mev.root");
  sig30mev->set_bdt("outputs/20000kev_1hits_combinedmasses/gpu_checkout_testing_20000kev_1hits_30mev.root");
  sig50mev->set_bdt("outputs/20000kev_1hits_combinedmasses/gpu_checkout_testing_20000kev_1hits_50mev.root");
  sig80mev->set_bdt("outputs/20000kev_1hits_combinedmasses/gpu_checkout_testing_20000kev_1hits_80mev.root");
  sig100mev->set_bdt("outputs/20000kev_1hits_combinedmasses/gpu_checkout_testing_20000kev_1hits_100mev.root");
  sig150mev->set_bdt("outputs/20000kev_1hits_combinedmasses/gpu_checkout_testing_20000kev_1hits_150mev.root");
  sig200mev->set_bdt("outputs/20000kev_1hits_combinedmasses/gpu_checkout_testing_20000kev_1hits_200mev.root");
  sig250mev->set_bdt("outputs/20000kev_1hits_combinedmasses/gpu_checkout_testing_20000kev_1hits_250mev.root");
  sig300mev->set_bdt("outputs/20000kev_1hits_combinedmasses/gpu_checkout_testing_20000kev_1hits_300mev.root");
  sig350mev->set_bdt("outputs/20000kev_1hits_combinedmasses/gpu_checkout_testing_20000kev_1hits_350mev.root");
  sig400mev->set_bdt("outputs/20000kev_1hits_combinedmasses/gpu_checkout_testing_20000kev_1hits_400mev.root");

  add_sample(sig15mev);
  add_sample(sig20mev);
  add_sample(sig30mev);
  add_sample(sig50mev);
  add_sample(sig80mev);
  add_sample(sig100mev);
  add_sample(sig150mev);
  add_sample(sig200mev);
  add_sample(sig250mev);
  add_sample(sig300mev);
  add_sample(sig350mev);
  add_sample(sig400mev);  
}

