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
  limit_curve(TString exp, TString title, Color_t col = kBlack, Bool_t pub = true);
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
  void SetMarkerStyle(int l) {
    this->get_tgraph()->SetMarkerStyle(l);
  }
  void SetMarkerColor(int l) {
    this->get_tgraph()->SetMarkerColor(l);
  }
};

void plot_vert_line(double x, TString label = "",double y = 1.4e-4)
{
  TLine *l1 = new TLine(x,9.5e-5,x,0.045);
  l1->SetLineStyle(2);
  l1->Draw();

  TLatex *tl = new TLatex();
  tl->SetTextFont(132);
  tl->DrawLatex(x-x*0.09,y,label);
}

// main
void plot_sensitivity()
{
  // base canvas
  TCanvas *c1 = new TCanvas("c1","c1",900,600);
  //TH2 * h2 = new TH2F("h2",";Mass m_{#chi} (MeV);Millicharge #epsilon = Q/e",10,10,10000,10,0.5e-4,2.5e-1);
  // zoomed in to make it look more like sensei limits figure
  TH2 * h2 = new TH2F("h2","Millicharge sensitivity;Mass m_{#chi} (MeV);Millicharge #epsilon = Q/e",10,10,4000,10,9.5e-5,0.045);
  h2->GetXaxis()->CenterTitle();
  h2->GetYaxis()->CenterTitle();
  c1->SetLogx();
  c1->SetLogy();
  h2->Draw();

  // published limits
  limit_curve argo("ArgoNeuT","ArgoNeuT",kTeal);
  limit_curve mini("MiniBooNE","MiniBooNE",kViolet);
  limit_curve mill("milliQan demonstrator","milliQan demonstrator",kBlue);
  limit_curve lsnd("LSND","LSND",kBlack);
  limit_curve lhcc("LHC","LHC",kGray);
  limit_curve slac("SLAC","SLAC",kRed);
  limit_curve sens("SENSEI","SENSEI",kOrange);
  limit_curve supk("Super K PB+MD","Super K PB+MD",kMagenta);
  limit_curve bebc("BEBC","BEBC",kAzure-5);
  limit_curve chrm("CHARM II","CHARM II",kOrange+3);

  TLegend *myleg = new TLegend(0.64,0.12,0.89,0.55);
  myleg->SetBorderSize(0);

  //mini.Draw("same");  myleg->AddEntry(mini.get_tgraph());
  //mill.Draw("same");  myleg->AddEntry(mill.get_tgraph());
  lsnd.Draw("same");  myleg->AddEntry(lsnd.get_tgraph());
  //lhcc.Draw("same");  myleg->AddEntry(lhcc.get_tgraph());
  slac.Draw("same"); myleg->AddEntry(slac.get_tgraph());
  argo.Draw("same");  myleg->AddEntry(argo.get_tgraph());
  sens.Draw("same");  myleg->AddEntry(sens.get_tgraph());
  //supk.Draw("same");  myleg->AddEntry(supk.get_tgraph());
  bebc.Draw("same");  myleg->AddEntry(bebc.get_tgraph());
  //chrm.Draw("same");  myleg->AddEntry(chrm.get_tgraph());
  


  limit_curve th_5000("15mev_5000kev","Wirecell 1-interaction 5000kev min. recoil",kGreen-2,false);
  th_5000.SetLineStyle(8);
  //th_5000.Draw("same LF"); myleg->AddEntry(th_5000.get_tgraph());

  limit_curve th_30000("15mev_30000kev","Wirecell 1-interaction 30000kev min. recoil",kRed-2,false);
  th_30000.SetLineStyle(8);
  //th_30000.Draw("same LF"); myleg->AddEntry(th_30000.get_tgraph());


  limit_curve th_20000("15mev_20000kev_pawelhistts","#splitline{Wirecell 1-interaction}{20000kev min. recoil}",kBlue-2,false);
  th_20000.SetLineStyle(8);
  //th_20000.Draw("same LF"); myleg->AddEntry(th_20000.get_tgraph());

  limit_curve th_20000_cm("20000kev_combinedmasses_syst",
			  //"#splitline{Wirecell 1-interaction}{All parents}",kMagenta,false);
			  "Wirecell 1-interaction with systematics",kMagenta,false);
  th_20000_cm.SetLineStyle(1); th_20000_cm.SetMarkerStyle(kFullTriangleUp); 
  th_20000_cm.Draw("same LP"); myleg->AddEntry(th_20000_cm.get_tgraph());

  limit_curve th_20000_noether("20000kev_combinedmasses_noether_240730",
			  //"#splitline{Wirecell 1-interaction}{All parents}",kMagenta,false);
			  "Wirecell 1-interaction",kMagenta-5,false);
  th_20000_noether.SetLineStyle(1); th_20000_noether.SetMarkerStyle(kFullTriangleDown); 
  th_20000_noether.Draw("same LP"); myleg->AddEntry(th_20000_noether.get_tgraph());

  
  /*
  limit_curve th_20000_cm_t70("20000kev_combinedmasses_truth70",
			  "Truth 70%",kBlue-2,false);
  th_20000_cm_t70.SetLineStyle(8); th_20000_cm_t70.SetMarkerStyle(kFullTriangleUp); 
  th_20000_cm_t70.Draw("same LFP"); myleg->AddEntry(th_20000_cm_t70.get_tgraph());
  */




  /*
    // per parent
  limit_curve pi0("20000kev_combinedmasses_per_parent_pi0",
		  "#pi^{0} only",kOrange,false);
  pi0.SetLineStyle(8); pi0.SetMarkerStyle(kFullDiamond); pi0.SetMarkerColor(kOrange);
  pi0.Draw("same LFP"); myleg->AddEntry(pi0.get_tgraph());

  limit_curve eta("20000kev_combinedmasses_per_parent_eta",
		  "#eta only",kRed,false);
  eta.SetLineStyle(8); eta.SetMarkerStyle(kFullTriangleUp); eta.SetMarkerColor(kRed); 
  eta.Draw("same LFP"); myleg->AddEntry(eta.get_tgraph());

  limit_curve etp("20000kev_combinedmasses_per_parent_etp",
		  "#eta' only",kBlue,false);
  etp.SetLineStyle(8); etp.SetMarkerStyle(kFullTriangleDown); etp.SetMarkerColor(kBlue); 
  etp.Draw("same LFP"); myleg->AddEntry(etp.get_tgraph());

  limit_curve rho("20000kev_combinedmasses_per_parent_rho",
		  "#rho only",kGreen,false);
  rho.SetLineStyle(8); rho.SetMarkerStyle(kFullCircle);rho.SetMarkerColor(kGreen); 
  rho.Draw("same LFP"); myleg->AddEntry(rho.get_tgraph());
  */

  
  
  /*
  // vertical lines for meson death
  plot_vert_line(67.5,  "#frac{#it{M}_{#pi^{0}}}{2}");
  plot_vert_line(279.,"#frac{#it{M}_{#eta}}{2}");
  plot_vert_line(385.0,"#frac{#it{M}_{#rho}}{2}",0.0007);
  plot_vert_line(478.89,"#frac{#it{M}_{#eta'}}{2}",0.004);
  */
  myleg->Draw("same");

  gStyle->SetOptStat(0);
  gSystem->Exec(Form("mkdir -p img/"));
  TDatime dt;
  c1->SaveAs(Form("img/sens_d%i_t%i.pdf",dt.GetDate(),dt.GetTime()));
  c1->SaveAs(Form("img/sens_d%i_t%i.png",dt.GetDate(),dt.GetTime()));

  
}

// class member implementation
limit_curve::limit_curve(){};
limit_curve::limit_curve(TString exp, TString title, Color_t col = kBlack, Bool_t pub = true)
{
  experiment = exp;
  published = pub;
  TString filename;
  if ( published ) {
    filename = "published/"+exp+".dat";
  } else {
    filename = "mysensitivities/"+exp+".txt";
  }
  dat_file.open(filename);
  read_file();
  set_color(col);
  fill_tgraph();
  this->SetTitle(title);
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
    limit->SetFillColorAlpha(color,0.1);
  }
}

void limit_curve::print_tgraph()
{
  for ( int i = 0; i < npoints; i++ ) {
    std::cout << limit->GetPointX(i) << "  " << limit->GetPointY(i) << std::endl;
  }
}

