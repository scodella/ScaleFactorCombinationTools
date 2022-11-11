TString CampaignName = "Run201525ns";
string CampaignNameString = "Run201525ns";
TString CampaignLuminosity = "2.4 fb^{-1}"; 
TString CenterOfMassEnergy = "13 TeV";

// Taggers
int const nTaggingAlgorithms = 2;
string TaggingAlgorithmName[nTaggingAlgorithms] = {"CSVv2", "JP"};

// Pt bins
int const nBinsCampaign = 7;
float MinPtCampaign, MaxPtCampaign;

float xPt[nBinsCampaign]  = {40., 60., 85., 120., 170., 250., 485.};
float exPt[nBinsCampaign] = {10., 10., 15.,  20.,  30.,  50., 185.};

// 
bool StatisticalCorrelation = true;
bool InflateStatistic       = false;
bool AddSampleDependence    = true;
float SampleDependence      = 0.01; 

// Measurements
int const nMeasurements = 8;
string MeasurementName[nMeasurements] = {"PtRel", "System8", "IP3D", "LT", "LT J/#Psi", "LT top", "KIN",  "TagCountTT"};
string MeasurementFlag[nMeasurements] = {"ptrel", "system8", "ip3d", "lt", "ltjpsi",    "lttop",  "kin",  "tctt"};
string MeasurementPlot[nMeasurements] = {"_PT",   "_S8",     "_IP", "_LT", "_PSI",      "_TT",    "_KIN", "_TCT"};
bool UseThisMeasurement[nTaggingAlgorithms][nMeasurements] = { {true, false, false, true,  false, false, false, false},
							       {true, true,  false, false, false, false, false, false} };
// Plot parameters
int GraphXBins[nMeasurements] = {7, 4, 3, 7, 9, 6, 5, 6};
int GraphStyle[nMeasurements] = {3005, -1, -1, -1, -1, -1, -1, -1};
int GraphColor[nMeasurements] = {kBlue, kGreen, kBlue, kRed, kOrange, kBlack, 28, kBlack};
int GraphMarker[nMeasurements] = {22, 29, 23, 20, 33, 21, 25, 21};
float GraphSize[nMeasurements] = {1.3, 1.3, 1.3, 1.3, 1., 1., 1., 1.};
int GraphWidth[nMeasurements] = {6, 6, 4, 4, 4, 4, 4, 4};
int MeasurementOrder[nMeasurements] = {4, 1, 2, 3, 5, 0, 6};

// Uncertainties
int const nUncertainties = 27;
string UncertaintyName[nUncertainties] = {"", "_statistic", "_pileup", "_mupt", "_gluonsplitting", "_bfragmentation", "_cfragmentation", "_jetaway", "_mudr", "_cb", "_jes", "_ttbarmodelling", "_l2c", "_ptrel", "_ipbias", "_ltothers", "_jer", "_tt", "_kt", "_tct", "_ntrk", "_njet", "_jeta", "_dmux", "_ksl", "_ntrkgen", "_bcorr"};

//bool IsForBreakdown[nUncertainties] = {false, true,         true,     true,    true,              true,              true,             true,       true,   true, true,   true,             true,   true,    true,     true,       true,  true, true, true,  true,   true,   true,   true,   true,  true,      true};
bool IsForBreakdown[nUncertainties] = {false, true,         false,     true,    true,              true,              false,             true,       false,   false, true,   false,             true,   true,     false,     false,       false,  false, false, false,  false,   false,   false,   false,   false,  false,      false };
bool IsPtCorrelated[nUncertainties] = {false, false,        true,      true,    false,             true,              false,             false,      true,     true, false,  false,             false,  false,    false,     true,        false,  false, false, false,  false,   false,   false,   false,   false,  false,      false};

int MeasurementSystematic[nMeasurements][nUncertainties] = { {1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1},  
							     {1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
							     {1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
							     {1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0},  
							     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},  
							     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},  
							     {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},  
							     {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0} };

// !!! Fraction of events with statistical correlation: without weighting
float FracPTJP[] = {0.096491, 0.128465, 0.162212, 0.182612, 0.203323, 0.235041, 0.239796};
float FracS8JP[] = {0.206224, 0.230399, 0.310995, 0.331386, 0.356070, 0.403861, 0.432426};
float FracPTS8[] = {0.467892, 0.488014, 0.413078, 0.434579, 0.44793,  0.401653, 0.413572};

// TTbar global measurement
float TTbarScale = 100.;
TString Chi2Strategy = "Overall";

// This is for avereging the results on ttbar spectrum:
float TTbarPt[nBinsCampaign] = {9.5, 42.5, 76., 56., 26.5, 4.5, 0.};
float TTbarSpectrumScale = 215.;
