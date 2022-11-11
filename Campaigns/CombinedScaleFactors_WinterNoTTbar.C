TString LegComp = "SF Winter no ttbar";

 double funSFb_Comp_TCHPT(double x) {
   double SFb = 0.671344*((1.+(0.111537*x))/(1.+(0.0796576*x)));
   return SFb; 
 }
 double funSFb_Comp_JPL(double x) {
   double SFb = 1.65936*((1.+(0.546853*x))/(1.+(0.93234*x)));
   return SFb; 
 }
 double funSFb_Comp_JPM(double x) {
   double SFb = 0.896758*((1.+(0.11259*x))/(1.+(0.105615*x)));
   return SFb; 
 }
 double funSFb_Comp_JPT(double x) {
   double SFb = 0.802128*((1.+(0.024399*x))/(1.+(0.0215834*x)));
   return SFb; 
 }
 double funSFb_Comp_CSVL(double x) {
   double SFb = 1.00572*((1.+(0.013676*x))/(1.+(0.0143279*x)));
   return SFb; 
 }
 double funSFb_Comp_CSVM(double x) {
   double SFb = (0.939158+(0.000158694*x))+(-2.53962e-07*(x*x));
   return SFb; 
 }
 double funSFb_Comp_CSVT(double x) {
   double SFb = (0.9203+(-3.32421e-05*x))+(-7.74664e-08*(x*x));
   return SFb; 
 }
 double funSFb_Comp_CSVV1L(double x) {
   double SFb = 1.27429*((1.+(0.385903*x))/(1.+(0.504169*x)));
   return SFb; 
 }
 double funSFb_Comp_CSVV1M(double x) {
   double SFb = 0.951164+(-2.13656e-05*x);
   return SFb; 
 }
 double funSFb_Comp_CSVV1T(double x) {
   double SFb = (0.903734+(0.000110656*x))+(-1.84386e-07*(x*x));
   return SFb; 
 }
 double funSFb_Comp_CSVSLV1L(double x) {
   double SFb = 0.986055*((1.+(0.000242583*x))/(1.+(0.000206821*x)));
   return SFb; 
 }
 double funSFb_Comp_CSVSLV1M(double x) {
   double SFb = (0.945328+(0.000207679*x))+(-5.22182e-07*(x*x));
   return SFb; 
 }
 double funSFb_Comp_CSVSLV1T(double x) {
   double SFb = (0.919663+(9.18523e-05*x))+(-3.72057e-07*(x*x));
   return SFb; 
 }  
