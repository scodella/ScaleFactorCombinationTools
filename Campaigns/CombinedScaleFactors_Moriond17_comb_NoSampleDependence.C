TString LegComp = "comb (no sample dependence)";
TString CompTitle = "_CompNoSampleDependence";
 
float xComp[9]  = {25., 40., 60., 85., 120., 170., 250., 450., 800.};
float exComp[9] = { 5., 10., 10., 15.,  20.,  30.,  50., 150., 200.};

double funSFb_Comp_CSVv2L(double xx) {
  double SFb = SFb = 0.887973*((1.+(0.0523821*xx))/(1.+(0.0460876*xx)));
  return SFb;
};  

float SFb_Comp_error_CSVv2L[] = {
  0.0229587, 
  0.00768633, 
  0.00581377, 
  0.00521596, 
  0.00413102, 
  0.0041927, 
  0.00890183, 
  0.0142091, 
  0.0197221
};

double funSFb_Comp_CSVv2M(double xx) {
  double SFb = 0.561694*((1.+(0.31439*xx))/(1.+(0.17756*xx)));
  return SFb;
}
float SFb_Comp_error_CSVv2M[] = {
  0.0370744, 
  0.00995909, 
  0.00735735, 
  0.00726552, 
  0.00564967, 
  0.00652676, 
  0.0111144, 
  0.0194543, 
  0.0321448,
};


double funSFb_Comp_CSVv2T(double xx) {
  double SFb = 0.817647*((1.+(0.038703*xx))/(1.+(0.0312388*xx)));

  return SFb;
}

float SFb_Comp_error_CSVv2T[] = {
  0.0326888, 
  0.0118847, 
  0.00910059, 
  0.0092945, 
  0.00867169, 
  0.0095481, 
  0.0167115, 
  0.0292375, 
  0.0530721
};


double funSFb_Comp_DeepCSVL(double xx) {
  double SFb = 0.954392*((1.+(0.00724505*xx))/(1.+(0.00666699*xx)));
  return SFb;
};  

float SFb_Comp_error_DeepCSVL[] = {
  0.0173998, 
  0.00571747, 
  0.00584571, 
  0.00487424, 
  0.00361892, 
  0.00417187, 
  0.00812666, 
  0.0154203, 
  0.0258248
};

double funSFb_Comp_DeepCSVM(double xx) {
  double SFb =  0.703522*((1.+(0.214774*xx))/(1.+(0.152657*xx)));

  return SFb;
}
float SFb_Comp_error_DeepCSVM[] = {
  0.0174673, 
  0.00864863, 
  0.0071851, 
  0.00628848, 
  0.00555223, 
  0.00533402, 
  0.011194, 
  0.0213006, 
  0.0285493
};


double funSFb_Comp_DeepCSVT(double xx) {
  double SFb = 0.534877*((1.+(0.305766*xx))/(1.+(0.16437*xx)));

  return SFb;
}

float SFb_Comp_error_DeepCSVT[] = {
  0.0251679, 
  0.0104582, 
  0.0083786, 
  0.00741188, 
  0.00785357, 
  0.00879909, 
  0.0161075, 
  0.0299195, 
  0.0369527
};

