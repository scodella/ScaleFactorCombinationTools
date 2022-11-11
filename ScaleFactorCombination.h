#include "Measurements/UL16APV/CampaignInfo.h"
//#include "Measurements/Prompt18/CampaignInfo.h"
//#include "Measurements/2016Legacy/CampaignInfo.h"
//#include "Measurements/DeepFlavour2017/CampaignInfo.h"
//#include "Measurements/Moriond18/CampaignInfo.h"

bool PrintPNG  =  true;
bool PrintPDF  =  true;
bool PrintC    =  true;
bool PrintRoot =  true;

bool etamin = false;
bool etamax = false;

bool PrintComparison             =  false;
bool SampleCombinationComparison =  false;
bool PrintOnlyTypeMeasurements   =   true;
bool PrintPtFit                  =   true;

bool StoreFittedScaleFactors =  true;
bool SystematicBreakdown     = false;
bool CategoryBreakdown       =  true;
bool PtCorrelationBreakdown  = false;
bool StatisticBreakdown      = false; 

int TaggingAlgorithm = -1;

enum { NBINS = 20 };

float SystematicCorrelation[nUncertainties][NBINS][nMeasurements*(nMeasurements-1)/2];

float MeasuredScaleFactor[nMeasurements][NBINS], MeasuredScaleFactorUncertainty[nMeasurements][NBINS][nUncertainties];

float xpt[nMeasurements][NBINS], expt[nMeasurements][NBINS]; 
float MeasuredScaleFactorValue[nMeasurements][NBINS], MeasuredScaleFactorError[nMeasurements][NBINS], MeasuredScaleFactorStatistic[nMeasurements][NBINS];

bool isUsedSystematic[nUncertainties];

int nTotalMeasurements;
int MeasurementBinIndex[1000];
int CampaignBinIndex[1000];
int nBinsCampaignForFit;
int FirstBinCampaignForFit;
int LastBinCampaignForFit;
//int nBinsCampaignForPlot;

TString TSfun1;

float sf[NBINS], sf_stat[NBINS], sf_error[NBINS], sf_uncerbreak[NBINS][50];
float fun_val[NBINS], fun_err[NBINS], fun_unc[NBINS][50];

double NormalizedChi2;

