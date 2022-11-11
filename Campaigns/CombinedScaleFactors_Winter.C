TString LegComp = "SF Winter";

 double funSFbWinter_TCHPT(double x) {
   double SFb = 0.703389*((1.+(0.088358*x))/(1.+(0.0660291*x)));
   return SFb; 
 }
 double funSFbWinter_CSVL(double x) {
    double SFb = 0.997942*((1.+(0.00923753*x))/(1.+(0.0096119*x)));
    return SFb; 
 }
 double funSFbWinter_CSVM(double x) {
   double SFb = (0.938887+(0.00017124*x))+(-2.76366e-07*(x*x));
   return SFb; 
 }
 double funSFbWinter_CSVT(double x) {
   double SFb = (0.927563+(1.55479e-05*x))+(-1.90666e-07*(x*x));
   return SFb; 
 }
 double funSFbWinter_CSVV1L(double x) {
   double SFb = 1.7586*((1.+(0.799078*x))/(1.+(1.44245*x)));
   return SFb; 
 }
 double funSFbWinter_CSVV1M(double x) {
   double SFb = 0.952067+(-2.00037e-05*x);
   return SFb; 
 }
 double funSFbWinter_CSVV1T(double x) {
   double SFb = (0.912578+(0.000115164*x))+(-2.24429e-07*(x*x));
   return SFb; 
 }
 double funSFbWinter_CSVSLV1L(double x) {
   double SFb = 0.970168*((1.+(0.00266812*x))/(1.+(0.00250852*x)));
   return SFb; 
 }
 double funSFbWinter_CSVSLV1M(double x) {
   double SFb = ((0.939238+(0.000278928*x))+(-7.49693e-07*(x*x)))+(2.04822e-10*(x*(x*x)));
   return SFb; 
 }
 double funSFbWinter_CSVSLV1T(double x) {
   double SFb = (0.928257+(9.3526e-05*x))+(-4.1568e-07*(x*x));
   return SFb; 
 }
