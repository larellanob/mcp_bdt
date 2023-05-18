void change_units()
{
  ifstream infile("published/BEBC gev.dat");
  ofstream outfile("published/BEBC.dat");
  while ( !infile.eof() )
    {
      double x,y;
      infile >> x >> y;
      std::cout << 1000*x << " " << sqrt(y) << std::endl;
      outfile << 1000*x << " " << y << std::endl;
    }
}
