TString CampaignName = "Run2015B";
string CampaignNameString = "Run2015B";
TString CampaignLuminosity = "40 pb^{-1}"; 
TString CenterOfMassEnergy = "13 TeV";


int const nTaggers = 3;
TString TaggerName[nTaggers] = {"CSVv2L", "CSVv2M", "CSVv2T"};

int const nBinsCampaign = 7;
float MinPtCampaign, MaxPtCampaign;

float xPt[nBinsCampaign]  = {40., 60., 85., 120., 170., 250., 485.};
float exPt[nBinsCampaign] = {10., 10., 15.,  20.,  30.,  50., 185.};

bool StatisticalCorrelation = true;
bool InflateStatistic       = false;
bool AddSampleDependence    = false;
float SampleDependence      = 0.00; 
float cJetsInflation[3]     = {2., 2., 2.};

int const nUncertainties = 26;
string UncertaintyName[nUncertainties] = {"", "_statistic", "_pileup", "_mupt", "_gluonsplitting", "_bfragmentation", "_cfragmentation", "_jetaway","_mudr", "_cb", "_jes", "_ttbarmodelling", "_l2c", "_ptrel", "_ipbias", "_ltothers", "_jer", "_tt", "_kt", "_tct", "_ntrk", "_njet", "_jeta", "_dmux", "_ksl", "_ntrkgen"};

int const nMeasurements = 8;
string MeasurementName[nMeasurements] = {"PtRel", "System8", "IP3D", "LT", "LT J/#Psi", "LT top", "Kin TT", "TagCountTT"};
string MeasurementFlag[nMeasurements] = {"ptrel", "system8", "ip3d", "lt", "ltjpsi", "lttop", "kintt", "tctt"};
int MeasurementSystematic[nMeasurements][nUncertainties] = { {1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},  
							     {1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
							     {1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},   
							     {1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1},  
							     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},  
							     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},  
							     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},  
							     {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0} };

bool UseThisMeasurement[nMeasurements] = {true, true, true, true, true, true, true, true};

int GraphXBins[nMeasurements] = {4, 4, 3, 7, 9, 6, 6, 6};
int GraphStyle[nMeasurements] = {3005, -1, -1, -1, -1, -1, -1, -1};
int GraphColor[nMeasurements] = {kBlue, kGreen, kBlue, kRed, kOrange, kBlack, 28, kBlack};
int GraphMarker[nMeasurements] = {22, 29, 23, 20, 33, 21, 25, 21};
float GraphSize[nMeasurements] = {1.3, 1.3, 1.3, 1.3, 1., 1., 1., 1.};
int GraphWidth[nMeasurements] = {6, 6, 4, 4, 4, 4, 4, 4};

// Fraction of events with statistical correlation
//float frac_PTJP[NBINS], frac_S8JP[NBINS], frac_PTS8[NBINS];
// Fraction of events with statistical correlation: without weighting


// Fraction of events with statistical correlation: without weighting
float frac_PTJP[] = {0.096491, 0.128465, 0.162212, 0.182612, 0.203323, 0.235041, 0.239796};
float frac_S8JP[] = {0.206224, 0.230399, 0.310995, 0.331386, 0.356070, 0.403861, 0.432426};
float frac_PTS8[] = {0.467892, 0.488014, 0.413078, 0.434579, 0.44793,  0.401653, 0.413572};
//float frac_PTJP[] = {0.467892, 0.488014, 0.413078, 0.434579, 0.44793,  0.401653, 0.413572};
//float frac_S8JP[] = {1., 1., 1., 1., 1., 1., 1.};
//float frac_PTS8[] = {0.467892, 0.488014, 0.413078, 0.434579, 0.44793,  0.401653, 0.413572};

// // Fraction of events with statistical correlation: with weighting
// frac_PTJP[0] = 0.214657; frac_SYJP[0] = 0.280682; frac_PTSY[0] = 0.727274;
// frac_PTJP[1] = 0.177045; frac_SYJP[1] = 0.293238; frac_PTSY[1] = 0.579063;
// frac_PTJP[2] = 0.205043; frac_SYJP[2] = 0.320358; frac_PTSY[2] = 0.610821;
// frac_PTJP[3] = 0.236862; frac_SYJP[3] = 0.414303; frac_PTSY[3] = 0.544984;
// frac_PTJP[4] = 0.259465; frac_SYJP[4] = 0.433469; frac_PTSY[4] = 0.570806;
// frac_PTJP[5] = 0.280611; frac_SYJP[5] = 0.455601; frac_PTSY[5] = 0.587821;
// frac_PTJP[6] = 0.282864; frac_SYJP[6] = 0.494757; frac_PTSY[6] = 0.546719;
// frac_PTJP[7] = 0.298543; frac_SYJP[7] = 0.516355; frac_PTSY[7] = 0.555483;
// frac_PTJP[8] = 0.293440; frac_SYJP[8] = 0.293440; frac_PTSY[8] = 0.517299;
// frac_PTJP[9] = 0.302531;
// frac_PTJP[10] = 0.315165;
// frac_PTJP[11] = 0.318102;
// frac_PTJP[12] = 0.227775;
// frac_PTJP[13] = 0.242758;
// frac_PTJP[14] = 0.239796;

// correlations (order:     PT/SY, PT/IP, PT/JP, PT/PSI, PT/TT, PT/KT, PT/TCT, SY/IP, SY/JP, SY/PSI, SY/TT, SY/KT/, SY/TCT, IP/JP, IP/PSI, IP/TT, IP/KT, IP/TCT, JP/PSI, JP/TT, JP/KT, JP/TCT, PSI/TT, PSI/KT, PSI/TCT, TT/KT, TT/TCT, KT/TCT)
//
float corrStat[]          = { 1,     1,     1,     0,      0,     0,     0,      1,     1,     0,      0,     0,      0,      0,     1,      0,     0,     0,      0,      0,     0,     0,       0,      0,      0,      1,     1,      1 };
float corrSystPU[] =    {-1,    -1,    -1,    -1,     -1,    -1,    -1,      1,     1,     1,      1,     1,      1,      1,     1,      1,     1,     1,      1,      1,      1,    1,       1,      1,      1,      1,     1,      1 };
float corrSystMuPt[]      =  {-1,    -1,     1,     0,      0,     0,     0,      1,    -1,     0,      0,     0,      0,     -1,     0,      0,     0,     0,      0,      0,      0,    0,       0,      0,      0,      0,     0,      0 };
float corrSystGluon[]     = { 1,     1,     1,     0,      0,     0,     0,      1,    -1,     0,      0,     0,      0,     -1,     0,      0,     0,     0,      0,      0,      0,    0,       0,      0,      0,      0,     0,      0 };
float corrSystBfrag[]     = {-1,    -1,     1,     1,      1,     0,     0,      1,    -1,    -1,     -1,     0,      0,     -1,    -1,     -1,     0,     0,      1,      1 ,     0,    0,       1,      0,      0,      0,     0,      0 };
float corrSystAway[]      =  { 1,     1,     0,     0,      0,     0,     0,      1,     0,     0,      0,     0,      0,      0,     0,      0,     0,     0,      0,      0,      0,    0,       0,      0,      0,      0,     0,     0 };
float corrSystMudr[]      = {-1,     1,     0,     0,      0,     0,     0,     -1,     0,     0,      0,     0,      0,      0,     0,      0,     0,     0,      0,      0,      0,    0,       0,      0,      0,      0,     0,     0 };
float corrSystLife[]      = { 0,     0,     0,     0,      0,     0,     0,      0,     0,     0,      0,     0,      0,      0,     0,      0,     0,     0,      1,      1,      0,    0,       1,      0,      0,      0,     0,     0 };
float corrSystJES[]       = { 0,     0,     0,     0,      0,     0,     0,      0,     0,     0,      0,     0,      0,      0,     0,      0,     0,     0,      0,      1,      1,    1,       0,      0,      0,      1,     1,     1 };
float corrSystTTModel[]   = { 0,     0,     0,     0,      0,     0,     0,      0,     0,     0,      0,     0,      0,      0,     0,      0,     0,     0,      0,      0,      0,    0,       0,      0,      0,      1,     1,     1 };

float TTbarScale = 100.;
TString Chi2str = "Overall";

// This is for avereging the results on ttbar spectrum:
float ttbar_pt[nBinsCampaign] = {9.5, 42.5, 76., 56., 26.5, 4.5, 0.};
float ttbar_tot = 215.;
