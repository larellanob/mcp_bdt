void fill_tgraph_vector(TString dir,   std::vector<TGraph *> &tgraphs);

double get_epsilon(double v0, int nhits)
{
  double exponent = 2*nhits+2;
  double powe0 = pow(0.001,exponent); // epsilon ^ 2nhits+2
  //double nevents {15.032}; ///////////////////// <<<<< 15 EVENT LIMIT
  double nevents = 1;
  return pow(powe0*nevents*v0,1./(exponent));
  //return v0;
}

double get_epsilon(double v0, int nhits, double mass)
{
  double exponent = 2*nhits+2;
  double powe0 = pow(0.001,exponent); // epsilon ^ 2nhits+2
  double nevents {15.032}; ///////////////////// <<<<< 15 EVENT LIMIT
  return pow(powe0*nevents*v0,1./(exponent));   
}

void plot_sensitivity(int threshold = 1000)
{

  // published data
  const char* dir = "/home/luciano/Physics/manchester/mcp_bdt/sensitivity/published/";
  std::vector<TGraph *> tgraphs;
  std::cout << "fill 1 " << std::endl;
  fill_tgraph_vector(dir,tgraphs);  

  // own limits
  //const char* dir_sim = Form("/uboone/app/users/arellano/bdt_mcp/mass_scan/sensitivity_2206_with_milliqan/simulated/%ikev/",threshold);
  const char* dir_sim = Form("./limits_apr_23/");
  std::vector<TGraph *> tgraphs_sim;
  std::cout << "fill 2 " << std::endl;
  fill_tgraph_vector(dir_sim,tgraphs_sim);
  
  auto * c1 = new TCanvas();
  c1->SetLogx();
  c1->SetLogy();

  
  //tgraphs[1]->SetLineColor(kBlue);
  //tgraphs[1]->SetLineStyle(2);
  //tgraphs[1]->SetLineWidth(5);


  tgraphs[5]->SetLineColor(kTeal);
  tgraphs[5]->SetLineWidth(5);

  tgraphs[4]->SetLineColor(kGray);
  tgraphs[4]->SetLineWidth(3);

  tgraphs[1]->SetLineColor(kRed);
  tgraphs[1]->SetLineWidth(3);

  tgraphs[0]->SetLineColor(kViolet);
  tgraphs[0]->SetLineWidth(3);

  tgraphs[2]->SetLineColor(kBlack);
  tgraphs[2]->SetLineWidth(3);

  tgraphs[3]->SetLineColor(kBlue);
  tgraphs[3]->SetLineWidth(3);

  tgraphs[1]->SetTitle("SLAC");
  //tgraphs[1]->SetTitle("MicroBooNE (estimate using BNB)");
  tgraphs[5]->SetTitle("ArgoNeuT");
  tgraphs[4]->SetTitle("LHC");
  tgraphs[0]->SetTitle("MiniBooNE");
  tgraphs[2]->SetTitle("LSND");
  tgraphs[3]->SetTitle("milliQan demonstrator");

  gStyle->SetOptStat(0);

  for ( size_t gr{0}; gr < tgraphs_sim.size(); gr++ ) {
    tgraphs_sim[gr]->SetLineColor(kGreen+gr+1);
    tgraphs_sim[gr]->SetLineWidth(3);
    tgraphs_sim[gr]->SetLineStyle(3+gr);
    tgraphs_sim[gr]->SetFillColorAlpha(kGreen+gr+1,0.2);
  }

  
  TH2 * h2 = new TH2F("h2",";Mass m_{#chi} (MeV);Millicharge #epsilon = Q/e",10,10,10000,10,0.5e-4,2.5e-1);
  h2->GetXaxis()->CenterTitle();
  h2->GetYaxis()->CenterTitle();
  h2->Draw();
  for ( size_t gr{0}; gr < tgraphs.size(); gr++ ) {
    tgraphs[gr]->Draw("same");
  }
  for ( size_t gr{0}; gr < tgraphs_sim.size(); gr++ ) {
    tgraphs_sim[gr]->Draw("same FL"); // F = filled, L = line
  }

  TLegend *myleg = new TLegend(0.6,0.15,0.9,0.45);
  myleg->AddEntry(tgraphs[0]);
  myleg->AddEntry(tgraphs[1]); // pheno estimate
  myleg->AddEntry(tgraphs[2]);
  myleg->AddEntry(tgraphs[3]);
  myleg->AddEntry(tgraphs[4]);
  myleg->AddEntry(tgraphs[5]);
  for ( int i = 0; i < tgraphs_sim.size(); i++ ) {
    myleg->AddEntry(tgraphs_sim[i]);
  }
  myleg->Draw();

  c1->SaveAs(Form("img/15ev_plot_%ikev.pdf",threshold));
  delete c1;
  
}
void fill_tgraph_vector(TString dir,   std::vector<TGraph *> &tgraphs)
{
  double nhits {1.};
  //const char* dir = "/uboone/app/users/arellano/bdt_mcp/mass_scan/sensitivity/published/";
  void *dirpointer = gSystem->OpenDirectory(dir);
  std::vector<TString> full_files;
  std::vector<TString> files;
  std::vector<TString> experiments;
  const char *entry;
  TString dir_iteration;
  while ( ( entry = (char*)gSystem->GetDirEntry(dirpointer))) {
    dir_iteration = entry;
    if ( dir_iteration.EndsWith(".") || dir_iteration == "README.txt" ) {
      continue;
    }
    full_files.push_back(dir+dir_iteration);
    files.push_back(dir_iteration);
    experiments.push_back(dir_iteration.ReplaceAll(".dat",""));
    std::cout << dir+dir_iteration << std::endl;
  }
  std::cout << "Found limits for the following experiments:" << std::endl;
  std::cout << experiments.size() << std::endl;
  for ( int i = 0; i < experiments.size(); i++ ) {
    std::cout << experiments[i] << std::endl;
  }
  std::vector<std::vector<Double_t>> experimentsx;
  std::vector<std::vector<Double_t>> experimentsy;
  std::vector<Double_t> datax;
  std::vector<Double_t> datay;
    for ( int i = 0; i < files.size(); i++ ) {
    Double_t x,y;
    std::cout << "opening file " << full_files[i] << std::endl;
    ifstream in;
    in.open(full_files[i]);
    datax.clear();
    datay.clear();
    while (1) {
      in >> x >> y;
      //std::cout << x << " " << y << std::endl;
      datax.push_back(x);
      datay.push_back(y);
      if (!in.good()) break;
    }
    experimentsx.push_back(datax);
    experimentsy.push_back(datay);
  }
  std::cout << experimentsx[0].size() << std::endl;

  for ( int exp = 0; exp < experiments.size(); exp++ ) {
    TGraph *gr = new TGraph(experimentsx[exp].size());
    tgraphs.push_back(gr);
  }

  for ( int exp = 0; exp < experiments.size(); exp++ ) {
    std::cout << "Graph " << experiments[exp] << std::endl;
    if ( experiments[exp].Contains("2hits") ) { nhits = 2; }
    if ( experiments[exp].Contains("3hits") ) { nhits = 3; }
    if ( experiments[exp].Contains("4hits") ) { nhits = 4; }
    if ( nhits != 1 ) { tgraphs[exp]->SetTitle(Form("%.0f hits",nhits)); }
    if ( experiments[exp].Contains("no-syst") ) { tgraphs[exp]->SetTitle("#muBooNE normal thresholds (no systematics)"); }
    if ( experiments[exp].Contains("with-bkg-syst") ) { tgraphs[exp]->SetTitle("#muBooNE normal thresholds (with bkg quadrature"); }
    if ( experiments[exp].Contains("limit") ) { tgraphs[exp]->SetTitle(Form("%.0f hits (limit)",nhits)); }
    std::cout << Form("using %.0f hits", nhits) << std::endl;
    for ( int point = 0; point < experimentsx[exp].size(); point++ ){
      tgraphs[exp]->SetPoint(point, experimentsx[exp][point], experimentsy[exp][point]);
      std::cout << tgraphs[exp]->GetMinimum() << " " << tgraphs[exp]->GetMaximum() << std::endl;
    }
  }
}
