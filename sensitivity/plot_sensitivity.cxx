class limit_curve {
private:
  ifstream dat_file;
  TString dat_filename;
  void fill_tgraph();
  void read_file();
  std::vector<double> x_points;
  std::vector<double> y_points;
  int npoints;
  TGraph *limit;
public:
  limit_curve();
  limit_curve(TString exp, Color_t col = kBlack, Bool_t pub = true);
  ~limit_curve() { delete limit;};
  TString experiment {"Experiment"};
  Color_t color {kBlack};
  Bool_t published {true}; // published by others as opposed to your own
  Bool_t filled {false};
  void set_exp(TString exp) {experiment = exp;};
  void set_color(Color_t col) {color = col; limit->SetLineColor(color);};
  void set_fill(Bool_t fill) {filled = fill;};
  TGraph * get_tgraph() { return limit; };
  void print_tgraph();
  // overloads from ROOT
  void SetFillColorAlpha(Color_t col, float alpha) {
    this->get_tgraph()->SetFillColorAlpha(col, alpha);
  }
  void Draw(const char * kOptions = "" ) {
    this->get_tgraph()->Draw(kOptions);
  }
  void SetTitle(const char * title = "") {
    this->get_tgraph()->SetTitle(title);
  }
  void SetLineStyle(int l) {
    this->get_tgraph()->SetLineStyle(l);
  }
};

// main
void plot_sensitivity()
{
  // base canvas
  TCanvas *c1 = new TCanvas();
  //TH2 * h2 = new TH2F("h2",";Mass m_{#chi} (MeV);Millicharge #epsilon = Q/e",10,10,10000,10,0.5e-4,2.5e-1);
  // zoomed in to make it look more like sensei limits figure
  TH2 * h2 = new TH2F("h2","MicroBooNE preliminary: MC + overlay;Mass m_{#chi} (MeV);Millicharge #epsilon = Q/e",10,10,4000,10,9.5e-5,0.045);
  h2->GetXaxis()->CenterTitle();
  h2->GetYaxis()->CenterTitle();
  c1->SetLogx();
  c1->SetLogy();
  h2->Draw();

  // published limits
  limit_curve argo("ArgoNeuT",kTeal);
  limit_curve mini("MiniBooNE",kViolet);
  limit_curve mill("milliQan demonstrator",kBlue);
  limit_curve lsnd("LSND",kBlack);
  limit_curve lhcc("LHC",kGray);
  limit_curve slac("SLAC",kRed);
  limit_curve sens("SENSEI",kOrange);
  limit_curve supk("Super K PB+MD",kMagenta);
  limit_curve bebc("BEBC",kAzure-5);
  limit_curve chrm("CHARM II",kOrange+3);

  mini.Draw("same");
  mill.Draw("same");
  lsnd.Draw("same");
  lhcc.Draw("same");
  slac.Draw("same");
  argo.Draw("same");
  sens.Draw("same");
  supk.Draw("same");
  bebc.Draw("same");
  chrm.Draw("same");
  

  // own limits (microboone)
  /*
  limit_curve syst("def-th_no-syst",kGreen+1,false);
  syst.SetLineStyle(1);
  syst.SetTitle("#muBooNE normal thresholds (no systematics)");
  syst.Draw("same LF");
  
  limit_curve no_syst("def-th_full-syst",kGreen+2,false);
  no_syst.SetLineStyle(2);
  no_syst.SetTitle("#muBooNE normal thresholds (with full systematics)");
  no_syst.Draw("same LF");
  */
  limit_curve gamma3d("gamma3d",kGreen+1,false);
  gamma3d.SetLineStyle(1);
  gamma3d.SetTitle("#muBooNE normal thresholds (gamma3d)");
  //gamma3d.Draw("same LF");
  limit_curve blip("blip",kGreen+2,false);
  blip.SetLineStyle(2);
  blip.SetTitle("#muBooNE normal thresholds (blipreco 30 step)");
  //blip.Draw("same LF");
  limit_curve blip10("blip_10rounds",kGreen+3,false);
  blip10.SetLineStyle(3);
  blip10.SetTitle("#muBooNE normal thresholds (blipreco 10 step)");
  //blip10.Draw("same LF");
  limit_curve blip50("blip_50rounds",kGreen+4,false);
  blip50.SetLineStyle(4);
  blip50.SetTitle("#muBooNE normal thresholds (blipreco 50 step)");
  //blip50.Draw("same LF");
  limit_curve blip100("blip_100rounds",kGreen-1,false);
  blip100.SetLineStyle(5);
  blip100.SetTitle("#muBooNE normal thresholds (blipreco 100 step)");
  //  blip100.Draw("same LF");
  limit_curve blipunbug("blip_unbug",kGreen-2,false);
  blipunbug.SetLineStyle(6);
  blipunbug.SetTitle("#muBooNE normal thresholds (blipreco unbugged)");
  //blipunbug.Draw("same LF");
  limit_curve blipunbug_100("blip_unbug_100",kGreen-3,false);
  blipunbug_100.SetLineStyle(7);
  blipunbug_100.SetTitle("MicroBooNE Blipreco3D");
  //blipunbug_100.Draw("same LF");
  limit_curve blipunbug_150("blip_unbug_150",kGreen-4,false);
  blipunbug_150.SetLineStyle(8);
  blipunbug_150.SetTitle("MicroBooNE 2-blips");
  blipunbug_150.Draw("same LF");

  
  // legend
  TLegend *myleg = new TLegend(0.6,0.15,0.9,0.5);
  myleg->AddEntry(argo.get_tgraph());
  myleg->AddEntry(mini.get_tgraph());
  myleg->AddEntry(mill.get_tgraph());
  myleg->AddEntry(lsnd.get_tgraph());
  myleg->AddEntry(lhcc.get_tgraph());
  myleg->AddEntry(slac.get_tgraph());
  myleg->AddEntry(sens.get_tgraph());
  myleg->AddEntry(supk.get_tgraph());
  myleg->AddEntry(bebc.get_tgraph());
  myleg->AddEntry(chrm.get_tgraph());
  //myleg->AddEntry(syst.get_tgraph());
  //myleg->AddEntry(no_syst.get_tgraph());
  /*
  myleg->AddEntry(gamma3d.get_tgraph());
  myleg->AddEntry(blip.get_tgraph());
  myleg->AddEntry(blip10.get_tgraph());
  myleg->AddEntry(blip50.get_tgraph());
  */
  //myleg->AddEntry(blip100.get_tgraph());
  //myleg->AddEntry(blipunbug.get_tgraph());
  //myleg->AddEntry(blipunbug_100.get_tgraph());
  myleg->AddEntry(blipunbug_150.get_tgraph());
  myleg->Draw("same");

  gStyle->SetOptStat(0);
  c1->SaveAs("img/limits_nov23.pdf");

  
}

// class member implementation
limit_curve::limit_curve(){};
limit_curve::limit_curve(TString exp, Color_t col = kBlack, Bool_t pub = true)
{
  experiment = exp;
  published = pub;
  TString filename;
  if ( published ) {
    filename = "published/"+exp+".dat";
  } else {
    filename = "limits_nov23/"+exp+".txt";
  }
  dat_file.open(filename);
  read_file();
  set_color(col);
  fill_tgraph();
}

void limit_curve::read_file()
{
  double x,y;
  while (1) {
    dat_file >> x >> y;
    x_points.push_back(x);
    y_points.push_back(y);
    if (!dat_file.good()) break;
  }
  if ( x_points.size() != y_points.size() ) {
    std::cout << Form("ERROR: sizes of x and y don't match, not filling exp %s",experiment.Data()) << std::endl;
  } else {
    npoints = x_points.size();
  }
  limit = new TGraph(npoints);
  if ( x_points.size() <2 ) {
    std::cout << Form("WARNING: Empty data for exp %s",experiment.Data()) << std::endl;
  }
}


void limit_curve::fill_tgraph()
{
  for ( int i = 0; i < x_points.size(); i++ ) {
    limit->SetPoint(i,x_points[i],y_points[i]);
  }
  limit->SetLineColor(color);
  limit->SetLineWidth(3);
  limit->SetTitle(experiment);
  if ( !published ) { // own limits filled
    limit->SetFillColorAlpha(color,0.2);
  }
}

void limit_curve::print_tgraph()
{
  for ( int i = 0; i < npoints; i++ ) {
    std::cout << limit->GetPointX(i) << "  " << limit->GetPointY(i) << std::endl;
  }
}

