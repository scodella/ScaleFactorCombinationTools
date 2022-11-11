TString LegComp = "SF Run2016 7.6/fb comb";
 
float xComp[7]  = {40., 60., 85., 120., 170., 250., 485.};
float exComp[7] = {10., 10., 15.,  20.,  30.,  50., 185.};


double funSFb_Comp_CSVv2L(double xx) {
  double SFb = 0.631151*((1.+(1.12333*xx))/(1.+(0.754135*xx)));
  return SFb;
};  

float SFb_Comp_error_CSVv2L[] = {
  0.00733156,
  0.00729381,
  0.00663374,
  0.0122059,
  0.0103633,
  0.0196129,
  0.0441719
};

double funSFb_Comp_CSVv2M(double xx) {
  double SFb = 0.671462*((1.+(0.548773*xx))/(1.+(0.404254*xx)));
  return SFb;
}
float SFb_Comp_error_CSVv2M[] = {
  0.00987782,
  0.00908035,
  0.00838405,
  0.0142745,
  0.0144621,
  0.0325831,
  0.0365897
};


double funSFb_Comp_CSVv2T(double xx) {
  double SFb = 0.859096+(3.11204e-05*xx);

  return SFb;
}

float SFb_Comp_error_CSVv2T[] = {
  0.0127854,
  0.0115238,
  0.0118904,
  0.0175575,
  0.020734,
  0.0475482,
  0.0422291
};

