// General data info
string CampaignNameString = "Moriond18CDE";
TString CampaignName = CampaignNameString;
TString CampaignLuminosity = "X.X fb^{-1}"; 
TString CenterOfMassEnergy = " (13 TeV, 2017)";

// Taggers
int const nTaggingAlgorithms = 3;
string TaggingAlgorithmName[nTaggingAlgorithms] = {"CSVv2", "DeepCSV", "DeepFlavour"};

// Pt bins
int const nBinsCampaign = 9;
float MinPtCampaign, MaxPtCampaign;
float xPt[nBinsCampaign]  = {25., 40., 60., 85., 120., 170., 250., 450., 800.};
float exPt[nBinsCampaign] = { 5., 10., 10., 15.,  20.,  30.,  50., 150., 200.};

// Some general uncertainty options
bool StatisticalCorrelation = true;
bool InflateStatistic       = false;
bool AddSampleDependence    = true;
float SampleDependence      = 0.01; 
float cJetsInflationFactor[3] = {2.5, 3.0, 3.5};

// Measurements
int const nMeasurements = 11;
string MeasurementName[nMeasurements] = {"PtRel",  "System8",  "LT",     "TagCount", "Kin",   "TnP",   "IF",    "TnPH",  "TnPL",  "mujets",  "ttbar"};
string MeasurementFlag[nMeasurements] = {"ptrel",  "system8",  "lt",     "TagCount", "kin",   "TnP",   "if",    "TnPH",  "TnPL",  "mujets",  "ttbar"};
string MeasurementPlot[nMeasurements] = {"_PT",    "_S8",      "_LT",    "_TC",      "_Kin",  "_TnP",  "_IF",   "_TnPH", "_TnPL", "_mujets", "_ttbar"};
string TypeMeasurement[nMeasurements] = {"mujets", "mujets",   "mujets", "ttbar",    "ttbar", "ttbar", "ttbar", "ttbar", "ttbar", "mujets",  "ttbar"};
bool UseThisMeasurement[nTaggingAlgorithms][nMeasurements] = { { true,  true,  true, false,  true, false,  false, false, false, false, false},
							       { true,  true,  true, false,  true, false,  false, false, false, false, false},
							       {false, false, false, false, false, false,  false, false, false, false, false} };
bool VetoedMeasurement [nTaggingAlgorithms][nMeasurements] = { {false, false, false, false, false, false,  false, false, false, false, false}, 
							       {false, false, false, false, false, false,  false, false, false, false, false}, 
							       {false, false, false, false, false, false,  false, false, false, false, false} };
bool PlotTheMeasurement[nTaggingAlgorithms][nMeasurements] = { { true,  true,  true, false,  true, false,  false, false, false, false, false},
							       { true,  true,  true, false,  true, false,  false, false, false, false, false},
							       {false, false, false, false, false, false,  false, false, false, false, false} };

// Plot parameters
int GraphXBins[nMeasurements] = {9, 5, 9, 6, 7, 6, 3, 6, 6, 9, 7};
int GraphStyle[nMeasurements] = {3005, -1, -1, -1, -1, -1, -1, -1, -1, -1};
int GraphColor[nMeasurements] = {kBlue, kGreen, kRed, kOrange, kBlack, 30, 46, 28, 30, kBlue, kRed};
int GraphMarker[nMeasurements] = {22, 29, 20, 21, 33, 23, 24, 23, 23, 21, 20};
float GraphSize[nMeasurements] = {1.3, 1.3, 1.3, 1., 1., 1., 1., 1., 1., 1., 1.};
int GraphWidth[nMeasurements] = {6, 6, 4, 4, 4, 4, 4, 4, 4, 4, 4};
int MeasurementOrder[nMeasurements] = {3, 2, 5, 7, 4, 6, 8, 9, 10, 0, 1};

// Uncertainties
int const nUncertainties = 43;
string UncertaintyName[nUncertainties] = {"", "_statistic", "_pileup", "_mupt", "_gluonsplitting", "_bfragmentation", "_cfragmentation", "_jetaway", "_mudr", "_cb", "_jes", "_ttbarmodelling", "_l2c", "_ptrel", "_ipbias", "_ltothers", "_jer", "_tt", "_kt", "_tct", "_ntrk", "_njet", "_jeta", "_dmux", "_ksl", "_ntrkgen", "_btempcorr", "_ltempcorr", "_flavFrac", "_bkg", "_scale1", "_ps", "_met", "_hpp", "_topmass", "_hdamp", "_qcdscale", "_sel", "_trig", "_tWth", "_isr", "_fsr", "_mistag"};

bool IsForBreakdown[nUncertainties] = {false, true,         false,     true,    true,              true,              false,             true,       false,   false, true,   false,             true,   true,     false,     false,       false,  false, false, false,  false,   false,   false,   false,   false,  false,      false,        false,         true,       true,     false,     false,  true,  true,   true,       false,    false,       false,  false,   false, false, false, false};
bool IsPtCorrelated[nUncertainties] = {false, false,        true,      true,    false,             true,              false,             false,      true,     true, false,  false,             false,  false,    false,     true,        false,  false, false, false,  false,   false,   false,   false,   false,  false,      false,        false,         true,       false,    false,     false,  false, false,  true,       false,    false,       false,  false,   false,  true,  true,  true};

int MeasurementSystematic[nMeasurements][nUncertainties] = { {1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
							     {1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
							     {1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},  
							     {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},    
							     {1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1},  
							     {1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0},  
							     {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},  
							     {1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0},  
							     {1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0},  
							     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},  
							     {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} };

// !!! Fraction of events with statistical correlation: without weighting
float FracPTJP[] = {0.09, 0.096491, 0.128465, 0.162212, 0.182612, 0.203323, 0.235041, 0.239796, 0.24};
float FracS8JP[] = {0.20, 0.206224, 0.230399, 0.310995, 0.331386, 0.356070, 0.403861, 0.432426, 0.43};
float FracPTS8[] = {0.46, 0.467892, 0.488014, 0.413078, 0.434579, 0.44793,  0.401653, 0.413572, 0.41};

// TTbar global measurement
TString Chi2Strategy = "Overall";

// This is for avereging the results on ttbar spectrum:
float TTbarPt[nBinsCampaign] = {0., 85., 50., 36., 26., 14., 2.1, 0., 0.};
float TTbarSpectrumScale;
