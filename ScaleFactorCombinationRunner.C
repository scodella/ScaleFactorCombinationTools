void ScaleFactorCombinationRunner(bool Compile, TString Algo, TString WP, TString MeasType) {
  
  if (Compile)
    gROOT->ProcessLine(".L ScaleFactorCombination.C++");
  else 
    gSystem->Load("ScaleFactorCombination_C.so");

  if (WP=="Loose"  || WP=="All") gROOT->ProcessLine(".x ScaleFactorCombinationExec.C(\"" + Algo + "\", \"Loose\",  \"" + MeasType + "\")");
  if (WP=="Medium" || WP=="All") gROOT->ProcessLine(".x ScaleFactorCombinationExec.C(\"" + Algo + "\", \"Medium\", \"" + MeasType + "\")");
  if (WP=="Tight"  || WP=="All") gROOT->ProcessLine(".x ScaleFactorCombinationExec.C(\"" + Algo + "\", \"Tight\",  \"" + MeasType + "\")");

}
  
