TString LegComp = "SF 13 TeV No S8";
 
float xComp[7]  = {40., 60., 85., 120., 170., 250., 485.};
float exComp[7] = {10., 10., 15.,  20.,  30.,  50., 185.};

double funSFb_Comp_CSVv2L(double xx) {
  double SFb = 0.967342*((1.+(0.000360684*xx))/(1.+(0.000268909*xx)));
  return SFb;
};  

float SFb_Comp_error_CSVv2L[] = {
  0.0212899,
  0.0237228,
  0.0216037,
  0.0224728,
  0.0206988,
  0.0196368,
  0.0439392
};

double funSFb_Comp_CSVv2M(double xx) {
  double SFb = 0.916291*((1.+(0.0118628*xx))/(1.+(0.0108104*xx)));
  return SFb;
}
float SFb_Comp_error_CSVv2M[] = {
  0.0249721,
  0.0339905,
  0.0349183,
  0.0222068,
  0.0174197,
  0.0277083,
  0.028926
};


double funSFb_Comp_CSVv2T(double xx) {
  double SFb = 0.834695*((1.+(0.00887379*xx))/(1.+(0.00715742*xx)));
  return SFb;
}

float SFb_Comp_error_CSVv2T[] = {
  0.035379,
  0.0403344,
  0.0409835,
  0.0318152,
  0.0306925,
  0.0277734,
  0.0442292
};

