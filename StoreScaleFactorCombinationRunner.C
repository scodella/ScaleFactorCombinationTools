void StoreScaleFactorCombinationRunner(bool Compile, string BTagger, string OP = "All", string MeasType = "mujets", string JetFlavour = "6", string StoreSystematicBreakdown = "NONE") {
  
  if (Compile)
    gROOT->ProcessLine(".L ScaleFactorCombination.C++");
  else 
    gSystem->Load("ScaleFactorCombination_C.so");

  const int nAlgorithms = 2;
  string Algorithms[nAlgorithms] = { "DeepCSV", "DeepJet" };

  const int nWorkingPoints = 3;
  string WorkingPoints[nWorkingPoints] = { "Loose", "Medium", "Tight" };

  const int nMeasTypes = 2;
  string MeasTypes[nWorkingPoints] = { "comb", "mujets" };

  for (int algo = 0; algo<nAlgorithms; algo++) {
    if (BTagger=="All" || BTagger==Algorithms[algo]) {
      for (int wp = 0; wp<nWorkingPoints; wp++) { 
        if (OP=="All" || OP==WorkingPoints[wp]) { 
          for (int mt = 0; mt<nMeasTypes; mt++) {
            if (MeasType=="All" || MeasType==MeasTypes[mt]) {

               TString StoreCommand =  "StoreScaleFactorCombination(\"" + Algorithms[algo] + "\", \"" + WorkingPoints[wp] + "\", \"" + MeasTypes[mt] + "\", " + JetFlavour + ", \"" + StoreSystematicBreakdown + "\")";  
               gROOT->ProcessLine(StoreCommand);

            }
          } 
        }
      }
    }
  } 

}

