TString CampaignName = "7TeVLegacy";
string CampaignNameString = "7TeVLegacy";
TString CampaignLuminosity = "5.1 fb^{-1}"; 
TString CenterOfMassEnergy = "7 TeV";


int const nTaggers = 4;
TString TaggerName[nTaggers] = {"TCHPT", "CSVL", "CSVM", "CSVT"};

int const nBinsCampaign = 15;
float MinPtCampaign, MaxPtCampaign;

float xPt[nBinsCampaign] = {25., 35., 45., 55., 65., 75., 90., 110., 140., 185., 235., 290., 360., 450., 585.};
float exPt[nBinsCampaign] = {5.,  5.,  5.,  5.,  5.,  5., 10.,  10.,  20.,  25.,  25.,  30.,  40.,  50.,  85.};

bool StatisticalCorrelation = true;
bool InflateStatistic       = false;
bool AddSampleDependence    = false;
float SampleDependence      = 0.00; 
float cJetsInflation[3]     = {2., 2., 2.};

//int const nUncertainties = 13;
//TString UncertaintyName[nUncertainties] = {"", "Statistic", "PileUp", "MuPt", "GluonSplitting", "BFragmentation", "JetAway","DeltaR", "LT-Bias", "LT-Cb", "JES", "ttbar-modelling", "Uncorrelated"};
//int const nUncertainties = 27;
//TString UncertaintyName[nUncertainties] = {"", "_statistic", "_pileup", "_mupt", "_gluonsplitting", "_bfragmentation", "_cfragmentation", "_jetaway","_mudr", "_ltbias", "_cb", "_jes", "_ttbarmodelling", "_l2c", "_ptrel", "_ipbias", "_ltothers", "_jer", "_tt", "_kt", "_tct", "_ntrk", "_njet", "_jeta", "_dmux", "_ksl", "_ntrkgen"};

/*
enum { NBINS = 20 };

int const nMeasurements = 8;
TString MeasurementName[nMeasurements] = {"PtRel", "System8", "IP3D", "LT", "LT J/#Psi", "LT top", "Kin TT", "Tag Count TT"};
bool UseThisMeasurement[nMeasurements] = {true, true, true, true, false, true, true, true};
float MeasuredScaleFactor[nMeasurements][NBINS], MeasuredScaleFactorError[nMeasurements][NBINS], MeasuredScaleFactorStatistic[nMeasurements][NBINS]; 

float xpt[nMeasurements+1][NBINS], expt[nMeasurements+1][NBINS]; 

float sf[NBINS], sf_stat[NBINS], sf_syst[NBINS], sf_eror[NBINS];

float eff_PT[NBINS], eff_stat_PT[NBINS], effMC_PT[NBINS], effMC_stat_PT[NBINS];
float sf_PT[NBINS], sf_stat_PT[NBINS], sf_syst_PT[NBINS], sf_eror_PT[NBINS];
float mupt_PT[NBINS], gluon_PT[NBINS], pu_PT[NBINS], total_syst_PT[NBINS];
float away_PT[NBINS], lc_PT[NBINS], bf_PT[NBINS], mudr_PT[NBINS];
float eff_PT_etaLT12[NBINS], eff_stat_PT_etaLT12[NBINS], effMC_PT_etaLT12[NBINS], effMC_stat_PT_etaLT12[NBINS], sf_PT_etaLT12[NBINS], sf_stat_PT_etaLT12[NBINS];
float eff_PT_etaGT12[NBINS], eff_stat_PT_etaGT12[NBINS], effMC_PT_etaGT12[NBINS], effMC_stat_PT_etaGT12[NBINS], sf_PT_etaGT12[NBINS], sf_stat_PT_etaGT12[NBINS];

float eff_S8[NBINS], eff_stat_S8[NBINS], effMC_S8[NBINS], effMC_stat_S8[NBINS];
float sf_S8[NBINS], sf_stat_S8[NBINS], sf_syst_S8[NBINS], sf_eror_S8[NBINS];
float mupt_S8[NBINS], gluon_S8[NBINS], pu_S8[NBINS], total_S8[NBINS];
float away_S8[NBINS], bef_S8[NBINS], mudr_S8[NBINS];
float ptrel_S8[NBINS], clos_S8[NBINS];
float eff_S8_etaLT12[NBINS], eff_stat_S8_etaLT12[NBINS], effMC_S8_etaLT12[NBINS], effMC_stat_S8_etaLT12[NBINS], sf_S8_etaLT12[NBINS], sf_stat_S8_etaLT12[NBINS];
float eff_S8_etaGT12[NBINS], eff_stat_S8_etaGT12[NBINS], effMC_S8_etaGT12[NBINS], effMC_stat_S8_etaGT12[NBINS], sf_S8_etaGT12[NBINS], sf_stat_S8_etaGT12[NBINS];

float eff_IP[NBINS], eff_stat_IP[NBINS], effMC_IP[NBINS], effMC_stat_IP[NBINS];
float sf_IP[NBINS], sf_stat_IP[NBINS], sf_syst_IP[NBINS], sf_eror_IP[NBINS];
float mupt_IP[NBINS], gluon_IP[NBINS], pu_IP[NBINS], total_syst_IP[NBINS];
float away_IP[NBINS], lc_IP[NBINS], bf_IP[NBINS], mudr_IP[NBINS];
float eff_IP_etaLT12[NBINS], eff_stat_IP_etaLT12[NBINS], effMC_IP_etaLT12[NBINS], effMC_stat_IP_etaLT12[NBINS], sf_IP_etaLT12[NBINS], sf_stat_IP_etaLT12[NBINS];
float eff_IP_etaGT12[NBINS], eff_stat_IP_etaGT12[NBINS], effMC_IP_etaGT12[NBINS], effMC_stat_IP_etaGT12[NBINS], sf_IP_etaGT12[NBINS], sf_stat_IP_etaGT12[NBINS];

float eff_JP[NBINS], eff_stat_JP[NBINS], effMC_JP[NBINS], effMC_stat_JP[NBINS];
float sf_JP[NBINS], sf_stat_JP[NBINS], sf_syst_JP[NBINS], sf_eror_JP[NBINS];
float mupt_JP[NBINS], gluon_JP[NBINS], bf_JP[NBINS], pu_JP[NBINS], total_syst_JP[NBINS];
float cor_JP[NBINS], inc_JP[NBINS], others_JP[NBINS], bias_JP[NBINS], away_JP[NBINS], jes_JP[NBINS], jer_JP[NBINS];

float eff_PSI[NBINS], eff_stat_PSI[NBINS], effMC_PSI[NBINS], effMC_stat_PSI[NBINS]; 
float sf_PSI[NBINS], sf_stat_PSI[NBINS], sf_syst_PSI[NBINS], sf_eror_PSI[NBINS];
float cosalpha_PSI[NBINS], deltaR_PSI[NBINS], deltaRprobe_PSI[NBINS]; 
float lxys_PSI[NBINS], mass_PSI[NBINS], ptsum_PSI[NBINS]; 
float z0_PSI[NBINS], pu_PSI[NBINS], bias_PSI[NBINS];
float cor_PSI[NBINS], bf_PSI[NBINS], other_PSI[NBINS];

float eff_TT[NBINS], eff_stat_TT[NBINS], effMC_TT[NBINS], effMC_stat_TT[NBINS];
float sf_TT[NBINS], sf_stat_TT[NBINS], sf_syst_TT[NBINS], sf_eror_TT[NBINS];
float pu_TT[NBINS], bf_TT[NBINS], total_syst_TT[NBINS];
float model_TT[NBINS], fit_TT[NBINS], jes_TT[NBINS];
float cor_TT[NBINS];

double Chi2Normal[NBINS];

*/

// Fraction of events with statistical correlation
//float frac_PTJP[NBINS], frac_S8JP[NBINS], frac_PTS8[NBINS];
// Fraction of events with statistical correlation: without weighting


// Fraction of events with statistical correlation: without weighting
float frac_PTJP[] = {0.117946, 0.096491, 0.112438, 0.128465, 0.144013, 0.159494, 
		     0.162212, 0.178839, 0.182612, 0.203323, 0.223828, 0.235041, 
		     0.227775, 0.242758, 0.239796};
float frac_S8JP[] = {0.193411, 0.206224, 0.230399, 0.310995, 0.331386, 0.356070,
		     0.403861, 0.432426, 0.476164};
float frac_PTS8[] = {0.609818, 0.467892, 0.488014, 0.413078, 0.434579, 0.44793,
		     0.401653, 0.413572, 0.383506};

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

// This is for avereging the results on ttbar spectrum (from 7 TeV combination code):
float ttbar_pt[15] = {0., 803050., 896000., 901200., 835650., 730950.,
		      1106000., 697600., 641000., 241100., 72950., 30670., 14380., 6183., 2761.};
float ttbar_tot = 6981000.;
