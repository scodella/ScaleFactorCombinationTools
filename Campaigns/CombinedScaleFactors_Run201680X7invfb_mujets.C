TString LegComp = "SF Run2016 7.6/fb mujets";
 
float xComp[7]  = {40., 60., 85., 120., 170., 250., 485.};
float exComp[7] = {10., 10., 15.,  20.,  30.,  50., 185.};


double funSFb_Comp_CSVv2L(double xx) {
  double SFb = 0.956768*((1.+(0.107522*xx))/(1.+(0.110892*xx)));
  return SFb;
};  

float SFb_Comp_error_CSVv2L[] = {
  0.00975533,
  0.0128151,
  0.0119232,
  0.0171376,
  0.01996,
  0.0197942,
  0.0447977
};

double funSFb_Comp_CSVv2M(double xx) {
  double SFb = 0.628408*((1.+(0.576276*xx))/(1.+(0.396304*xx)));
  return SFb;
}
float SFb_Comp_error_CSVv2M[] = {
  0.0146231,
  0.0158532,
  0.0151635,
  0.0230163,
  0.0297416,
  0.0329956,
  0.0370984
};


double funSFb_Comp_CSVv2T(double xx) {
  double SFb = 0.850182*((1.+(-0.000611408*xx))/(1.+(-0.000654205*xx)));

  return SFb;
}

float SFb_Comp_error_CSVv2T[] = {
  0.0216494,
  0.0206873,
  0.0278909,
  0.0345319,
  0.0504979,
  0.0480927,
  0.0440766
};

