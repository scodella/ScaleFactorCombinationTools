TString LegComp = "SF Run2016 4/fb";
 
float xComp[7]  = {40., 60., 85., 120., 170., 250., 485.};
float exComp[7] = {10., 10., 15.,  20.,  30.,  50., 185.};



double funSFb_Comp_CSVv2L(double xx) {
  double SFb = 0.752528*((1.+(0.459115*x))/(1.+(0.367048*x)));
  return SFb;
};  

float SFb_Comp_error_CSVv2L[] = {
  0.018552,
  0.932649,
  0.018653,
  0.0188339,
  0.0274406,
  0.0456386,
  0.0604846,
  0.0767694
};

double funSFb_Comp_CSVv2M(double xx) {
  double SFb = 0.918081;
  return SFb;
}
float SFb_Comp_error_CSVv2M[] = {
  0.0177431,
  0.017958,
  0.018256,
  0.0185481,
  0.0183675,
  0.0194489,
  0.0240062
};


double funSFb_Comp_CSVv2T(double xx) {
  double SFb = 0.13067*((1.+(1.84506*x))/(1.+(0.268922*x)));
  return SFb;
}

float SFb_Comp_error_CSVv2T[] = {
  0.831369,
  0.0185548,
  0.0190936,
  0.0431566,
  0.0576586,
  0.0798476,
  0.0894715
};

