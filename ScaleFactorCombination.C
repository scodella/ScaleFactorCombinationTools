#include <iostream>
#include <fstream>
#include <TROOT.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TGraphErrors.h>
#include <TH2F.h>
#include <TF1.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TText.h>
#include <TString.h>
#include <TLatex.h>

#include "TMatrixDSym.h"
#include "TVectorD.h"

#include "ScaleFactorCombination.h"
//#include "BTagCalibrationStandalone.cc"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

float TakeMaximum(float up, float down) {

  float ThisMaximum = fabs(up);
  if (fabs(down)>ThisMaximum) ThisMaximum = fabs(down);
  if (up<0.) ThisMaximum *= -1.;

  return ThisMaximum;

}

float TakeAverage(float up, float down) {

  float ThisAverage = (fabs(up) + fabs(down))/2.;
  if (up<0.) ThisAverage *= -1.;
  
  return ThisAverage;

}

//#include "Campaigns/CombinedScaleFactors_Run201680X7invfb_mujets.C"
//#include "Campaigns/CombinedScaleFactors_Run201680X7invfb_comb.C"
#include "Campaigns/CombinedScaleFactors_Moriond17_comb_NoSampleDependence.C"

TMatrixD BuildCovarianceMatrix(TString ErrorCategory = "", bool isFit = true) {

  TMatrixD MatCov(nTotalMeasurements, nTotalMeasurements);

  MatCov *= 0.;

  int RowIndex = 0;
  for (int im1 = 0; im1<nTotalMeasurements; im1++) {

      int ColumnIndex = 0;
      for (int im2 = 0; im2<nTotalMeasurements; im2++) {

	int jm1 = MeasurementBinIndex[im1];
	int jm2 = MeasurementBinIndex[im2];
	
	int bpt1 = CampaignBinIndex[im1];
	int bpt2 = CampaignBinIndex[im2];
	
	// This term contains systematics specific of a method, which are not correlated with other methods, for which we can ignore the bin-to-bin correlation, and for which we do not need to compute the contribution to the total uncertainty
	float SpecificUncertainty = 0.;
	if (jm1==jm2) SpecificUncertainty = pow(MeasuredScaleFactorUncertainty[jm1][bpt1][0], 2);
	
	for (int is = 1; is<nUncertainties; is++) {
	  
	  // "" if we are making the fit
	  // Specifying a systematic to compute its contribution to the total uncertatiny
	  if ((ErrorCategory=="" || ErrorCategory==UncertaintyName[is]) && isUsedSystematic[is]) {
	    
	    float ThisMatrixElement = 0.;
	    
	    float ThisMatrixTerm = MeasuredScaleFactorUncertainty[jm1][bpt1][is]*MeasuredScaleFactorUncertainty[jm2][bpt2][is];

	    if (bpt1==bpt2 // Always correlate a systematic in the same pt bin
		|| IsPtCorrelated[is] || !isFit) { // Systematic to be correlated across the pt bins
	      
	      // !!! Statistical uncertainty needs special treatment for correlation across
	      // the methods
	      if (UncertaintyName[is]=="_statistic") {
		
		if (jm1==jm2) 
		  ThisMatrixElement = fabs(ThisMatrixTerm);
		else {
		  
		  if ((MeasurementName[jm1]=="PtRel" && MeasurementName[jm2]=="System8") ||
		      (MeasurementName[jm2]=="PtRel" && MeasurementName[jm1]=="System8")) {
		    
		    ThisMatrixElement = sqrt(FracPTS8[bpt1])*fabs(ThisMatrixTerm);
		    
		  } else if ((MeasurementName[jm1]=="PtRel" && MeasurementName[jm2]=="IP3D") ||
			     (MeasurementName[jm2]=="PtRel" && MeasurementName[jm1]=="IP3D")) {
		    
		    ThisMatrixElement = fabs(ThisMatrixTerm);
		    
		  } else if ((MeasurementName[jm1]=="PtRel" && MeasurementName[jm2]=="LT") ||
			     (MeasurementName[jm2]=="PtRel" && MeasurementName[jm1]=="LT")) {
		    
		    ThisMatrixElement = sqrt(FracPTS8[bpt1])*FracS8JP[bpt1]*fabs(ThisMatrixTerm);
		    
		  } else if ((MeasurementName[jm1]=="System8" && MeasurementName[jm2]=="IP3D") ||
			     (MeasurementName[jm2]=="System8" && MeasurementName[jm1]=="IP3D")) {
		    
		    ThisMatrixElement = sqrt(FracPTS8[bpt1])*fabs(ThisMatrixTerm);
		    
		  } else if ((MeasurementName[jm1]=="LT" && MeasurementName[jm2]=="IP3D") ||
			     (MeasurementName[jm2]=="LT" && MeasurementName[jm1]=="IP3D")) {
		    
		    ThisMatrixElement = sqrt(FracPTS8[bpt1]*FracS8JP[bpt1])*fabs(ThisMatrixTerm);
		    
		  } else if ((MeasurementName[jm1]=="System8" && MeasurementName[jm2]=="LT") ||
			     (MeasurementName[jm2]=="System8" && MeasurementName[jm1]=="LT")) {
		    
		    ThisMatrixElement = sqrt(FracS8JP[bpt1])*fabs(ThisMatrixTerm);
		    
		  } 
		  
		}
		
	      }
	      // Now look at the other systematics. The following list must contains:
	      // - systematics correlated across the pt bin (or between measurements)
	      // - systematics for which we want to provide the specific contribution to the total uncertainty
	      // Any other systematic can be left to the SpecificUncertainty term
	      else if (ThisMatrixTerm!=0. && (jm1!=jm2 || IsPtCorrelated[is] || IsForBreakdown[is] || ErrorCategory!="")) {
		
		float Coefficient = 1.;
		if (UncertaintyName[is]=="_ltothers" && bpt1!=bpt2) Coefficient = 0.5; // !!! Special for LT systematics
		
		//if ((UncertaintyName[is]!="_jetaway" && UncertaintyName[is]!="_jes") || jm1==jm2 || !isFit) // !!! Special case for jetaway systematic (as in Run1)
		if (UncertaintyName[is]!="_jetaway" || jm1==jm2 || !isFit) // !!! Special case for jetaway systematic (as in Run1)
		  ThisMatrixElement = Coefficient*ThisMatrixTerm;
		
	      } 
	      
	      MatCov(RowIndex, ColumnIndex) += ThisMatrixElement;
	      
	      if (ThisMatrixElement!=0. && jm1==jm2 && bpt1==bpt2) 
		SpecificUncertainty -= ThisMatrixTerm;
	      
	    }
	    
	  } 
	  
	}
	
	if (ErrorCategory=="") {
	  if (jm1==jm2 && bpt1==bpt2) {
	    if (SpecificUncertainty>0.00001) {  
	      MatCov(RowIndex, ColumnIndex) += SpecificUncertainty;
	    } else {
	      cout << " Warning: SpecificUncertainty too small (" << SpecificUncertainty<< ") for " << MeasurementName[jm1] << " in " << xPt[bpt1] << endl;
	    }
	  }
	}

	ColumnIndex++;
	
      }
      
      RowIndex++;
      
  }
  
  return MatCov;

}

void ScaleFactorsFit() {
  
  cout << "  Build MatU" << endl;

  TMatrixD MatU(nTotalMeasurements, nBinsCampaignForFit);

  MatU *= 0.;

  int RowIndex = 0;
  for (int im = 0; im<nTotalMeasurements; im++) {
    
    int bpt = CampaignBinIndex[im] - FirstBinCampaignForFit;

    int jm = MeasurementBinIndex[im];
    if (MeasurementName[jm]!="TagCountTT") MatU(RowIndex, bpt) = 1.;
    else {
      
      for (int bpt2 = 0; bpt2<nBinsCampaignForFit; bpt2++) 
	MatU(RowIndex, bpt2) = TTbarPt[FirstBinCampaignForFit+bpt2]/TTbarSpectrumScale;

    }
    
    RowIndex++;
    
  }
  
  cout << "  Build MatCov" << endl;

  TMatrixD MatCov = BuildCovarianceMatrix();
  TMatrixD MatCovStatistics = BuildCovarianceMatrix("_statistic");
  
  cout << "  Build VecSF" << endl;
  
  TMatrixD VecSF(nTotalMeasurements, 1);

  VecSF *= 0.;

  RowIndex = 0;
  for (int im = 0; im<nTotalMeasurements; im++) {

    int jm = MeasurementBinIndex[im];
    int bpt = CampaignBinIndex[im];
    
    VecSF(RowIndex, 0) = MeasuredScaleFactor[jm][bpt];
    
    RowIndex++;
    
  }
  
  cout << "  Make fit" << endl;
  
  TMatrixD MatUTr(nBinsCampaignForFit, nTotalMeasurements); 
  MatUTr.Transpose(MatU);
  
  double Determinant;
  TMatrixD MatCovInv(MatCov);
  MatCovInv.Invert(&Determinant);
  if (Determinant==0.)
    cout << "      MatCov not invertible" << endl;
  
  TMatrixD Aux1(nBinsCampaignForFit, nTotalMeasurements);
  Aux1.Mult(MatUTr, MatCovInv);
  
  TMatrixD Aux2(nBinsCampaignForFit, nBinsCampaignForFit);
  Aux2.Mult(Aux1, MatU); 
  Aux2.Invert(&Determinant);

  if (Determinant==0.)
    cout << "    Aux2 not invertible" << endl;
  
  TMatrixD MatCoefficients(nBinsCampaignForFit, nTotalMeasurements);
  MatCoefficients.Mult(Aux2, Aux1);

  cout << "  Get SFs" << endl;

  TMatrixD MatSF(nBinsCampaignForFit, 1);
  MatSF.Mult(MatCoefficients, VecSF); 

  cout << "  Get SF errors" << endl;

  TMatrixD MatCoefficientsTr(nTotalMeasurements, nBinsCampaignForFit);
  MatCoefficientsTr.Transpose(MatCoefficients);

  TMatrixD Aux3(nBinsCampaignForFit, nTotalMeasurements);
  Aux3.Mult(MatCoefficients, MatCov);

  TMatrixD MatSFErr(nBinsCampaignForFit, nBinsCampaignForFit);
  MatSFErr.Mult(Aux3, MatCoefficientsTr);

  TMatrixD Aux3Statistics(nBinsCampaignForFit, nTotalMeasurements);
  Aux3Statistics.Mult(MatCoefficients, MatCovStatistics);

  TMatrixD MatSFErrStatistics(nBinsCampaignForFit, nBinsCampaignForFit);
  MatSFErrStatistics.Mult(Aux3Statistics, MatCoefficientsTr);

  cout << "  Fill fit results" << endl;

  for (int bpt = FirstBinCampaignForFit; bpt<=LastBinCampaignForFit; bpt++) {

    int MatrixBin = bpt - FirstBinCampaignForFit;

    sf[bpt] = MatSF(MatrixBin, 0);

    sf_error[bpt] = sqrt(MatSFErr(MatrixBin, MatrixBin)); 
    sf_stat[bpt] = sqrt(MatSFErrStatistics(MatrixBin, MatrixBin));

    sf_uncerbreak[bpt][0] = sf_error[bpt];
    sf_uncerbreak[bpt][1] = sf_stat[bpt];

  }
 
  if (SystematicBreakdown || CategoryBreakdown) {

    for (int is = 2; is<nUncertainties; is++) 
      if (isUsedSystematic[is] && ((SystematicBreakdown && IsForBreakdown[is]) || (CategoryBreakdown && IsForCategoryBreakdown[is]))) {

	TMatrixD MatCovSystematic = BuildCovarianceMatrix(UncertaintyName[is], true);
	
	TMatrixD Aux3Systematic(nBinsCampaignForFit, nTotalMeasurements);
	Aux3Systematic.Mult(MatCoefficients, MatCovSystematic);
		
	TMatrixD MatSFErrSystematic(nBinsCampaignForFit, nBinsCampaignForFit);
	MatSFErrSystematic.Mult(Aux3Systematic, MatCoefficientsTr);

	for (int bpt = FirstBinCampaignForFit; bpt<=LastBinCampaignForFit; bpt++) {

	  /*if (UncertaintyName[is]=="_jer") {
	    double CC = 0.;
	    for (int k = 0 ; k<nTotalMeasurements; k++)
	      for (int r = 0 ; r<nTotalMeasurements; r++) {
		double CCC = MatCoefficients(bpt - FirstBinCampaignForFit, k)*MatCoefficients(bpt - FirstBinCampaignForFit, r)*MatCovSystematic(k,r);
		if (CCC!=0.) {
		  cout << " MMM " << MeasurementName[MeasurementBinIndex[k]] << " " << xPt[CampaignBinIndex[k]] << " " << MeasurementName[MeasurementBinIndex[r]] << " " << xPt[CampaignBinIndex[r]] << " " << CCC  << endl;
		  CC += CCC;
		}
	      }
	    cout << " CC " << CC << endl;
	    }*/
  
	  if (MatSFErrSystematic(bpt - FirstBinCampaignForFit, bpt - FirstBinCampaignForFit)<0.) 
	    cout << " Warning: negative covariance matrix element " << MatSFErrSystematic(bpt - FirstBinCampaignForFit, bpt - FirstBinCampaignForFit) << endl;
	  
	  double unsigned_sf_uncerbreak = sqrt(MatSFErrSystematic(bpt - FirstBinCampaignForFit, bpt - FirstBinCampaignForFit));
	  
	  double signed_sf_uncerbreak = 0.;
	  for (int im = 0; im<nTotalMeasurements; im++) {

	      int jm = MeasurementBinIndex[im];
	      int ptm  = CampaignBinIndex[im];
	      double sigma = MeasuredScaleFactorUncertainty[jm][ptm][is];
	      double alpha = MatCoefficients[bpt - FirstBinCampaignForFit][im];
	      signed_sf_uncerbreak += alpha*sigma;

	  }

	  if (IsPtCorrelated[is]) {
	    sf_uncerbreak[bpt][is] = signed_sf_uncerbreak;
	    if (abs(unsigned_sf_uncerbreak-abs(signed_sf_uncerbreak))>1.0e-08)
	      cout << " WWW " << UncertaintyName[is] << " " << xPt[bpt] << " " << unsigned_sf_uncerbreak << " " << signed_sf_uncerbreak << " " << unsigned_sf_uncerbreak-abs(signed_sf_uncerbreak) << endl;
	  } else {
	    /*cout << " ### " << UncertaintyName[is] << " " << xPt[bpt] << " " << unsigned_sf_uncerbreak << " " << signed_sf_uncerbreak << " " << unsigned_sf_uncerbreak-abs(signed_sf_uncerbreak) << " " << sf_error[bpt] << " " << FirstBinCampaignForFit << endl;
	    for (int im = 0; im<nTotalMeasurements; im++) {
	      if (CampaignBinIndex[im]==bpt) {
	      int jm = MeasurementBinIndex[im];
	      if (MeasuredScaleFactorUncertainty[jm][bpt][is]>0.)
		cout << "               " << MeasurementName[jm] << " " << MeasuredScaleFactorUncertainty[jm][bpt][is] << "     " << MatCoefficients[bpt - FirstBinCampaignForFit][im] << endl;
	      else if (MeasuredScaleFactorUncertainty[jm][bpt][is]<0.)
		cout << "        --->   " << MeasurementName[jm] << " " << MeasuredScaleFactorUncertainty[jm][bpt][is] <<  "    " << MatCoefficients[bpt - FirstBinCampaignForFit][im] << endl;
	      }
	      }*/
	    sf_uncerbreak[bpt][is] = unsigned_sf_uncerbreak;
	  }

	}
	
      }
    
  }

  if (PtCorrelationBreakdown) {

    for (int is = 2; is<nUncertainties; is++) 
      if (isUsedSystematic[is] && IsPtCorrelated[is]) {
	
	TMatrixD MatCovSystematic = BuildCovarianceMatrix(UncertaintyName[is]);
	
	TMatrixD Aux3Systematic(nBinsCampaignForFit, nTotalMeasurements);
	Aux3Systematic.Mult(MatCoefficients, MatCovSystematic);
		
	TMatrixD MatSFErrSystematic(nBinsCampaignForFit, nBinsCampaignForFit);
	MatSFErrSystematic.Mult(Aux3Systematic, MatCoefficientsTr);

	for (int bpt = FirstBinCampaignForFit; bpt<=LastBinCampaignForFit; bpt++) {
  
	  sf_uncerbreak[bpt][is] = sqrt(MatSFErrSystematic(bpt - FirstBinCampaignForFit, bpt));
	  
	}
	
      }
    
  }

  cout << "  Compute Chi2" << endl;

  NormalizedChi2 = 0.;
  for (int im1 = 0; im1<nTotalMeasurements; im1++)
    for (int im2 = 0; im2<nTotalMeasurements; im2++) 
      for (int bpt1 = FirstBinCampaignForFit; bpt1<=LastBinCampaignForFit; bpt1++) 
	for (int bpt2 = FirstBinCampaignForFit; bpt2<=LastBinCampaignForFit; bpt2++) {
	  
	  int jm1 = MeasurementBinIndex[im1];
	  int jm2 = MeasurementBinIndex[im2];
	  
	  int MatIdx1 = bpt1 - FirstBinCampaignForFit;
	  int MatIdx2 = bpt2 - FirstBinCampaignForFit;

	  if (MeasurementName[jm1]!="TagCountTT" && MeasurementName[jm2]!="TagCountTT") {
	    
	    if (bpt1==CampaignBinIndex[im1] && bpt2==CampaignBinIndex[im2]) 
	      NormalizedChi2 += MatU(im1, MatIdx1)*(MeasuredScaleFactor[jm1][bpt1] - sf[bpt1])*MatCovInv(im1, im2)*MatU(im2, MatIdx2)*(MeasuredScaleFactor[jm2][bpt2] - sf[bpt2]);
	    
	  } else if (MeasurementName[jm1]=="TagCountTT" && MeasurementName[jm2]=="TagCountTT") {

	    NormalizedChi2 += MatU(im1, MatIdx1)*(MeasuredScaleFactor[jm1][0] - sf[bpt1])*MatCovInv(im1, im2)*MatU(im2, MatIdx2)*(MeasuredScaleFactor[jm2][0] - sf[bpt2]);
	    
	  } else if (MeasurementName[jm1]=="TagCountTT" && MeasurementName[jm2]!="TagCountTT") {

	    if (bpt2==CampaignBinIndex[im2]) 
	      NormalizedChi2 += MatU(im1, MatIdx1)*(MeasuredScaleFactor[jm1][0] - sf[bpt1])*MatCovInv(im1, im2)*MatU(im2, MatIdx2)*(MeasuredScaleFactor[jm2][bpt2] - sf[bpt2]);
	    
	  } else if (MeasurementName[jm1]!="TagCountTT" && MeasurementName[jm2]=="TagCountTT") {

	    if (bpt2==CampaignBinIndex[im1]) 
	      NormalizedChi2 += MatU(im1, MatIdx1)*(MeasuredScaleFactor[jm1][bpt1] - sf[bpt1])*MatCovInv(im1, im2)*MatU(im2, MatIdx2)*(MeasuredScaleFactor[jm2][0] - sf[bpt2]);
	    
	  }
	  
	}
  
  NormalizedChi2 /= (nTotalMeasurements - nBinsCampaignForFit);
  cout << "Normalized Chi2 = " << NormalizedChi2 << endl;
  
}

TGraphErrors *MakeTGraphErrors(int nBins, float VectX[NBINS], float VectY[NBINS], float VectXError[NBINS], float VectYError[NBINS], int ThisStyle, int ThisColor, int ThisMarker, float ThisSize, int ThisWidth) {
  
  TGraphErrors *ThisTGraphErrors = new TGraphErrors(nBins, VectX, VectY, VectXError, VectYError);

  if (ThisStyle>=0) ThisTGraphErrors->SetFillStyle(ThisStyle);
  ThisTGraphErrors->SetFillColor(ThisColor);
  ThisTGraphErrors->SetMarkerColor(ThisColor);
  ThisTGraphErrors->SetLineColor(ThisColor);
  ThisTGraphErrors->SetLineStyle(1);
  ThisTGraphErrors->SetMarkerStyle(ThisMarker);
  ThisTGraphErrors->SetMarkerSize(ThisSize);
  ThisTGraphErrors->SetLineWidth(ThisWidth);

  return ThisTGraphErrors;

}

void ReadMeasurements(string  BTagger, TString OP, TString MeasType) {

  BTagEntry::OperatingPoint WP;
  if (OP=="Loose") WP = BTagEntry::OP_LOOSE;
  if (OP=="Medium") WP = BTagEntry::OP_MEDIUM;
  if (OP=="Tight") WP = BTagEntry::OP_TIGHT;

  BTagEntry::JetFlavor JF = BTagEntry::FLAV_B;

  nTotalMeasurements = 0; nBinsCampaignForFit = 0; FirstBinCampaignForFit = -1; //nBinsCampaignForPlot = 0;

  for (int is = 0; is<nUncertainties; is++) isUsedSystematic[is] = false;

  for (int nm = 0; nm<nMeasurements; nm++) 
    if (UseThisMeasurement[TaggingAlgorithm][nm]) {

      string CSVFileName = "./Measurements/" + CampaignNameString + "/" + BTagger + "_" + MeasurementName[nm] + ".csv";
      BTagCalibration calib(BTagger, CSVFileName);
      
      //BTagCalibrationReader *CentralReader = new BTagCalibrationReader(&calib, WP, MeasurementFlag[nm], "central");
      std::vector<std::string> v_syst;
      for (int is = 0; is<nUncertainties; is++) {
        if (MeasurementSystematic[nm][is]==1) {
          string thisSyst = UncertaintyName[is];
          v_syst.push_back("up"+thisSyst);
          v_syst.push_back("down"+thisSyst);
        }
      }  

      BTagCalibrationReader *CentralReader = new BTagCalibrationReader(WP, "central", v_syst);      
      CentralReader->load(calib, JF, MeasurementFlag[nm]);      

      int nMeasBins = -1; float LastMeasuredSF = 0., LowBinEdge = 0.;

      for (int bpt = 0; bpt<nBinsCampaign; bpt++) {

        MeasuredScaleFactor[nm][bpt] = CentralReader->eval(JF, 0., xPt[bpt]);
	
	if (MeasuredScaleFactor[nm][bpt]!=LastMeasuredSF) {
	  
	  if (nMeasBins>=0) {
	    
	    xpt[nm][nMeasBins] = (LowBinEdge + xPt[bpt] - exPt[bpt])/2.;
	    expt[nm][nMeasBins] = xpt[nm][nMeasBins] - LowBinEdge;
	    MeasuredScaleFactorValue[nm][nMeasBins] = LastMeasuredSF;
	    
	  }
	  
	  nMeasBins++;

	  if (!VetoedMeasurement[TaggingAlgorithm][nm] && MeasuredScaleFactor[nm][bpt]!=0. && 
	      (MeasType=="comb" || TypeMeasurement[nm]==MeasType) && 
	      !(MeasurementName[nm]=="System8" && xPt[bpt]>100. && CampaignName.Contains("Moriond18")) && // Carefully here...
	      !(MeasurementName[nm]=="PtRel" && xPt[bpt]>200. && CampaignName.Contains("Moriond18") && CampaignName!="Moriond18B") && // Carefully here...
	      !(MeasurementName[nm]=="PtRel" && xPt[bpt]>600. && CampaignName=="UL17")    &&  // Carefully here... 
              !(MeasurementName[nm]=="PtRel" && xPt[bpt]>600. && CampaignName=="UL18")    &&  // Carefully here...  
              //!(MeasurementName[nm]=="PtRel" && xPt[bpt]>600. && CampaignName=="UL16APV") &&  // Carefully here...
              !(MeasurementName[nm]=="PtRel" && xPt[bpt]>600. && CampaignName=="UL16")        // Carefully here...     
	      ) {

	    for (int is = 0; is<nUncertainties; is++) 
	      if (MeasurementSystematic[nm][is]!=0)
		isUsedSystematic[is] = true;
	    
	    MeasurementBinIndex[nTotalMeasurements] = nm;
	    CampaignBinIndex[nTotalMeasurements] = bpt;
	    nTotalMeasurements++;
	    
	    if (FirstBinCampaignForFit<0) FirstBinCampaignForFit = bpt;

	  }

	  LowBinEdge = xPt[bpt] - exPt[bpt];

	  LastMeasuredSF = MeasuredScaleFactor[nm][bpt];
	  
	}      

      }

      if (LastMeasuredSF>0.) {

	xpt[nm][nMeasBins] = (LowBinEdge + xPt[nBinsCampaign-1] + exPt[nBinsCampaign-1])/2.;
	expt[nm][nMeasBins] = xpt[nm][nMeasBins] - LowBinEdge;
	MeasuredScaleFactorValue[nm][nMeasBins] = LastMeasuredSF;

	nMeasBins++;

      }

      if (nMeasBins>nBinsCampaignForFit && !VetoedMeasurement[TaggingAlgorithm][nm] && 
	  (MeasType=="comb" || TypeMeasurement[nm]==MeasType)) 
	nBinsCampaignForFit = nMeasBins;

      //if (PlotTheMeasurement[TaggingAlgorithm][nm])
      //if (nMeasBins>nBinsCampaignForPlot) 
      //  nBinsCampaignForPlot = nMeasBins;

      float TotalErrorForFit[50], TotalCSVErrorForFit[50];
      for (int bpt = 0; bpt<nBinsCampaign; bpt++) {
	
	TotalErrorForFit[bpt] = 0.;
	TotalCSVErrorForFit[bpt] = 0.;
	
      }
      
      float TotalError[50], TotalCSVError[50];
      for (int bpt = 0; bpt<nMeasBins; bpt++) {
	
	TotalError[bpt] = 0.;
	TotalCSVError[bpt] = 0.;
	
      }
      
      bool hasSystematicScaling = false;

      for (int is = 0; is<nUncertainties; is++) { 

	for (int bpt = 0; bpt<nBinsCampaign; bpt++)
	  MeasuredScaleFactorUncertainty[nm][bpt][is] = 0.;
	
	if (MeasurementSystematic[nm][is]!=0) {

	  string UpFlag = "up" + UncertaintyName[is], DownFlag = "down" + UncertaintyName[is];
	  //BTagCalibrationReader *SystUpReader = new BTagCalibrationReader(&calib, WP, MeasurementFlag[nm], UpFlag);
	  //BTagCalibrationReader *SystDownReader = new BTagCalibrationReader(&calib, WP, MeasurementFlag[nm], DownFlag);
          
	  for (int bpt = 0; bpt<nBinsCampaign; bpt++) {

            if (MeasuredScaleFactor[nm][bpt]==0.) continue;

	    //float UpSF = SystUpReader->eval(JF, 0., xPt[bpt]); UpSF -= MeasuredScaleFactor[nm][bpt];
	    //float DownSF = SystDownReader->eval(JF, 0., xPt[bpt]); DownSF -= MeasuredScaleFactor[nm][bpt];
            float UpSF = CentralReader->eval_auto_bounds(UpFlag, JF, 0.1, xPt[bpt]); UpSF -= MeasuredScaleFactor[nm][bpt];
            float DownSF = CentralReader->eval_auto_bounds(DownFlag, JF, 0.1, xPt[bpt]); DownSF -= MeasuredScaleFactor[nm][bpt];
	    
            //if (UpSF<0. || DownSF>0.) {
	    //  if (IsForCategoryBreakdown[is]) {
            //    if (IsPtCorrelated[is])
	    //      cout << " FFF      " << MeasurementName[nm] << " " << UncertaintyName[is] << " " << bpt << " " << UpSF << " " << DownSF << endl;
	    //	  else
	    //	    cout << " FFF ---> " << MeasurementName[nm] << " " << UncertaintyName[is] << " " << bpt << " " << UpSF << " " << DownSF << endl;
	    //  }
	    //}
	    if (fabs(UpSF)>0.5 || fabs(DownSF)>0.5) {
              cout << " FFF      " << BTagger << " " << OP << " " << MeasurementName[nm] << " " << UncertaintyName[is] << " " << bpt << " " << MeasuredScaleFactor[nm][bpt] << " " << UpSF << " " << DownSF << endl; 
            }

	    if (CampaignName=="Run201680X7invfb" && UncertaintyName[is]=="" && MeasurementFlag[nm]=="kin") {
	      
	      float FakeSystematic = MeasuredScaleFactor[nm][bpt]*0.016;

	      float PowUpSF = UpSF*UpSF;
	      PowUpSF += FakeSystematic*FakeSystematic;
	      UpSF = sqrt(PowUpSF)*UpSF/fabs(UpSF);
	      
	      float PowDownSF = DownSF*DownSF;
	      PowDownSF += FakeSystematic*FakeSystematic;
	      DownSF = sqrt(PowDownSF)*DownSF/fabs(DownSF);
	      
	      MeasuredScaleFactorUncertainty[nm][bpt][28] = FakeSystematic;
	      
	    } else if (CampaignName=="Moriond17paper" && UncertaintyName[is]=="_fsr" && 
		       (MeasurementFlag[nm]=="TnP" || MeasurementFlag[nm]=="TnPH" || MeasurementFlag[nm]=="TnPL") ) {
	      UpSF *= sqrt(2.); DownSF *= sqrt(2.);
	      hasSystematicScaling = true;
	    } 
	    
	    float ThisError =  TakeAverage(UpSF, DownSF);
	    MeasuredScaleFactorUncertainty[nm][bpt][is] = ThisError;
	    if (UncertaintyName[is]=="") TotalCSVErrorForFit[bpt] = ThisError;
	    else TotalErrorForFit[bpt] += ThisError*ThisError;

	  }

	  for (int bpt = 0; bpt<nMeasBins; bpt++) {
	    
	    //float UpSF = SystUpReader->eval(JF, 0., xpt[nm][bpt]); UpSF -= MeasuredScaleFactorValue[nm][bpt];
	    //float DownSF = SystDownReader->eval(JF, 0., xpt[nm][bpt]); DownSF -= MeasuredScaleFactorValue[nm][bpt];
            float UpSF = CentralReader->eval_auto_bounds(UpFlag, JF, 0.1,  xpt[nm][bpt]); UpSF -= MeasuredScaleFactorValue[nm][bpt];
            float DownSF = CentralReader->eval_auto_bounds(DownFlag, JF, 0.1,  xpt[nm][bpt]); DownSF -= MeasuredScaleFactorValue[nm][bpt];

	    if (CampaignName=="Run201680X7invfb" && UncertaintyName[is]=="" && MeasurementFlag[nm]=="kin") {
	      
	      float FakeSystematic = MeasuredScaleFactor[nm][bpt]*0.016;

	      float PowUpSF = UpSF*UpSF;
	      PowUpSF += FakeSystematic*FakeSystematic;
	      UpSF = sqrt(PowUpSF)*UpSF/fabs(UpSF);
	      
	      float PowDownSF = DownSF*DownSF;
	      PowDownSF += FakeSystematic*FakeSystematic;
	      DownSF = sqrt(PowDownSF)*DownSF/fabs(DownSF);

	    } else if (CampaignName=="Moriond17paper" && UncertaintyName[is]=="_fsr" && 
		       (MeasurementFlag[nm]=="TnP" || MeasurementFlag[nm]=="TnPH" || MeasurementFlag[nm]=="TnPL") ) {
	      UpSF *= sqrt(2.); DownSF *= sqrt(2.);
	      hasSystematicScaling = true;
	    }
	   
	    float ThisError = TakeAverage(UpSF, DownSF);
	    if (UncertaintyName[is]=="_statistic") MeasuredScaleFactorStatistic[nm][bpt] = ThisError;
	    if (UncertaintyName[is]=="") TotalCSVError[bpt] = ThisError;
	    else TotalError[bpt] += ThisError*ThisError;

 	  }

	}
	
      }

      //if (MeasurementName[nm]=="LT") hasSystematicScaling = true;
      //if (MeasurementName[nm]=="TnPH") hasSystematicScaling = true;
      //if (MeasurementName[nm]=="TnPL") hasSystematicScaling = true;
      //if (MeasurementName[nm]=="Kin") hasSystematicScaling = true;

      for (int bpt = 0; bpt<nBinsCampaign; bpt++) {
	
	if (TotalCSVErrorForFit[bpt]!=0. && !hasSystematicScaling)
	  MeasuredScaleFactorUncertainty[nm][bpt][0] = TotalCSVErrorForFit[bpt];
	else 
	  MeasuredScaleFactorUncertainty[nm][bpt][0] = sqrt(TotalErrorForFit[bpt]);

	if (TotalCSVErrorForFit[bpt]!=0. && !hasSystematicScaling && sqrt(TotalErrorForFit[bpt])>TotalCSVErrorForFit[bpt]) {
	  int toFix = 0;
	  for (int is = 2; is<nUncertainties; is++) 
	    if (abs(MeasuredScaleFactorUncertainty[nm][bpt][is])>0.7)
	      toFix++;
	  if (toFix==1)
	    for (int is = 2; is<nUncertainties; is++) 
	      if (abs(MeasuredScaleFactorUncertainty[nm][bpt][is])>0.7) {
		cout << " toFix " << MeasuredScaleFactorUncertainty[nm][bpt][is] << endl;
		MeasuredScaleFactorUncertainty[nm][bpt][is] = TotalCSVErrorForFit[bpt]*TotalCSVErrorForFit[bpt];
		for (int is2 = 2; is2<nUncertainties; is2++)
		  if (is2!=is)
		    MeasuredScaleFactorUncertainty[nm][bpt][is] -= MeasuredScaleFactorUncertainty[nm][bpt][is2]*MeasuredScaleFactorUncertainty[nm][bpt][is2];
		MeasuredScaleFactorUncertainty[nm][bpt][is] = sqrt(MeasuredScaleFactorUncertainty[nm][bpt][is]);
		cout << " toFix " << MeasuredScaleFactorUncertainty[nm][bpt][is] << endl;
	      }
	}

      }

      for (int bpt = 0; bpt<nMeasBins; bpt++) {
	
	if (TotalCSVError[bpt]!=0. && !hasSystematicScaling)// && MeasurementName[nm]!="PtRel" && MeasurementName[nm]!="System8") 
	  MeasuredScaleFactorError[nm][bpt] = TotalCSVError[bpt];
	else 
	  MeasuredScaleFactorError[nm][bpt] = sqrt(TotalError[bpt]);

      }
      
      if (InflateStatistic) { // Run1 stuff ...

	for (int bpt = 0; bpt<nMeasBins; bpt++) {

	  int BinMultiplicity = 0;
	  for (int cpt = 0; cpt<nBinsCampaign; cpt++) 
	    if (xPt[cpt]>xpt[nm][bpt]-expt[nm][bpt] && xPt[cpt]<xpt[nm][bpt]+expt[nm][bpt])
	      BinMultiplicity++;
  
	  for (int cpt = 0; cpt<nBinsCampaign; cpt++) 
	    if (xPt[cpt]>xpt[nm][bpt]-expt[nm][bpt] && xPt[cpt]<xpt[nm][bpt]+expt[nm][bpt]) 
	      MeasuredScaleFactorUncertainty[nm][cpt][1] *= sqrt(BinMultiplicity);

	}

      }

      for (int bpt = 0; bpt<nMeasBins; bpt++) {

	float xPtShift = 0.;
	if (nm==1) xPtShift =  1.;
	if (nm==2) xPtShift = -1.;
	if (nm==3) xPtShift =  2.;
	if (nm==4) xPtShift = -2.;
	if (MeasurementName[nm]=="mujets") xPtShift = 1;
	if (MeasurementName[nm]=="ttbar") xPtShift = -1;
	if (expt[nm][bpt]>=15.) xPtShift *= 2.;
	else if (xpt[nm][bpt]>160.) xPtShift *= 3.;
	xpt[nm][bpt] += xPtShift;
 
      }

    }

  LastBinCampaignForFit = FirstBinCampaignForFit + nBinsCampaignForFit - 1;
  
  MinPtCampaign = xPt[FirstBinCampaignForFit] - exPt[FirstBinCampaignForFit];
  MaxPtCampaign = xPt[LastBinCampaignForFit] + exPt[LastBinCampaignForFit];

  TTbarSpectrumScale = 0.;
  for (int bpt = FirstBinCampaignForFit; bpt<=LastBinCampaignForFit; bpt++)
    TTbarSpectrumScale += TTbarPt[bpt];

}
  
void SetAlgorithmToUse(string BTagger) {

  TaggingAlgorithm = -1;

  for (int tg = 0; tg<nTaggingAlgorithms; tg++) 
    if (BTagger==TaggingAlgorithmName[tg]) TaggingAlgorithm = tg;
  
}

TStyle* PlotStyle() {

  TStyle *plotStyle = new TStyle("PLOT","Plot style CMS");
  
  plotStyle->SetOptTitle(0);  
  plotStyle->SetOptStat(0); 
  plotStyle->SetOptFit(0);
  
  plotStyle->SetPaperSize(20.,26.);
  
  plotStyle->SetEndErrorSize(2);
  plotStyle->SetErrorX(0.);
  
  plotStyle->SetFrameBorderMode(0);
  plotStyle->SetFrameFillColor(0);
  plotStyle->SetFrameFillStyle(0);
  plotStyle->SetCanvasBorderMode(0);
  plotStyle->SetFillColor(0);
  plotStyle->SetCanvasColor(0);
  plotStyle->SetCanvasBorderSize(2);
  
  plotStyle->SetPadBorderMode(0);
  plotStyle->SetPadColor(0);
  plotStyle->SetPadGridX(false);
  plotStyle->SetPadGridY(false);
  plotStyle->SetGridColor(0);
  plotStyle->SetGridStyle(3);
  plotStyle->SetGridWidth(1);
  
  plotStyle->SetCanvasDefX(0);
  plotStyle->SetCanvasDefY(0);
  
  plotStyle->SetHistLineColor(1);
  plotStyle->SetHistLineStyle(0);
  plotStyle->SetHistLineWidth(1);
  
  plotStyle->SetPadTickX(1);
  plotStyle->SetPadTickY(1);
  
  plotStyle->SetPadLeftMargin(0.16);
  plotStyle->SetPadRightMargin(0.02);
  plotStyle->SetPadTopMargin(0.06);
  plotStyle->SetPadBottomMargin(0.13);
  
  plotStyle->SetLabelFont(42,"x");
  plotStyle->SetTitleFont(42,"x");
  plotStyle->SetLabelFont(42,"y");
  plotStyle->SetTitleFont(42,"y");
  plotStyle->SetLabelFont(42,"z");
  plotStyle->SetTitleFont(42,"z");
  
  plotStyle->SetLabelSize(0.05,"x");
  plotStyle->SetTitleSize(0.06,"x");
  plotStyle->SetLabelSize(0.05,"y");
  plotStyle->SetTitleSize(0.06,"y");
  plotStyle->SetLabelSize(0.035,"z");
  plotStyle->SetTitleSize(0.035,"z");
  
  plotStyle->SetLabelOffset(0.007,"XYZ");
  
  plotStyle->SetNdivisions(510,"XYZ");
  plotStyle->SetStripDecimals(kTRUE);
  plotStyle->SetTickLength(0.03, "XYZ");
  
  plotStyle->SetTitleXOffset(0.9);
  plotStyle->SetTitleYOffset(1.25);
  
  return plotStyle;

}

void ScaleFactorCombination(string BTagger, TString OP, TString MeasType) {

  // @@@ Define some general parameter and read the measurements
  
  SetAlgorithmToUse(BTagger);

  TString title = BTagger;
  if (OP=="Loose")  title += "L"; 
  if (OP=="Medium") title += "M"; 
  if (OP=="Tight")  title += "T"; 

  float exstat[NBINS];
  for (int bpt = 0; bpt<nBinsCampaign; bpt++) exstat[bpt] = 0.0001;
    
  if (!StatisticalCorrelation) 
    for (int bpt = 0; bpt<nBinsCampaign; bpt++) {
      FracPTJP[bpt] = 1.; 
      FracS8JP[bpt] = 1.; 
      FracPTS8[bpt] = 1.;
    }
  
  ReadMeasurements(BTagger, OP, MeasType);

  // @@@ Performe the combination fit
  
  cout << "NOW STARTING THE FITS " << endl;

  ScaleFactorsFit();

  cout << "FITS DONE, NOW MAKING PLOTS" << endl;

  float xPt_fit[NBINS], exPt_fit[NBINS], sf_fit[NBINS], sf_error_fit[NBINS];
  for (int bpt = 0; bpt<nBinsCampaign; bpt++) {
    xPt_fit[bpt] = xPt[bpt];
    exPt_fit[bpt] = exPt[bpt];
    sf_fit[bpt] = sf[bpt];
    sf_error_fit[bpt] = sf_error[bpt];
  }
  
  // Some campaign related adjustement
  if (CampaignName=="Run201680X4invfb") {

    if  (OP!="Medium" && BTagger=="CSVv2") {
    
      for (int bpt = 3; bpt<nBinsCampaign; bpt++) {
	
	sf[bpt] = sf[1];
	if (fabs(sf[1] - MeasuredScaleFactorValue[0][bpt])>sf_error[bpt])
	  sf_error[bpt] = fabs(sf[1] - MeasuredScaleFactorValue[0][bpt]);
	if (fabs(sf[1] - MeasuredScaleFactorValue[2][bpt])>sf_error[bpt])
	  sf_error[bpt] = fabs(sf[1] - MeasuredScaleFactorValue[2][bpt]);
	
      }
      
    }

    for (int bpt = 0; bpt<nBinsCampaign; bpt++) {
      if (sf_error[bpt]/sf[bpt]<0.02) sf_error[bpt] = 0.02*sf[bpt];
      sf_error_fit[bpt] = sf_error[bpt];
    }

  } else if (CampaignName.Contains("Moriond18")) {

    bool DowngradeHightPt = false;/*true;
    if (CampaignName=="Moriond18CDE" && OP=="Loose" && MeasType=="comb") DowngradeHightPt = false;
    if (CampaignName=="Moriond18B"   && OP=="Loose") DowngradeHightPt = false;
    if (CampaignName=="Moriond18EF"  && OP=="Loose") DowngradeHightPt = false;
    if (CampaignName=="Moriond18"  && BTagger=="CSVv2" && OP!="Medium") DowngradeHightPt = false;
    if (CampaignName=="Moriond18B" && BTagger=="CSVv2" && OP!="Medium") DowngradeHightPt = false;
    if (CampaignName=="Moriond18CDE" && BTagger=="CSVv2" && OP!="Medium") DowngradeHightPt = false;
    if (CampaignName=="Moriond18EF" && BTagger=="CSVv2" && OP!="Medium") DowngradeHightPt = false;*/
    //if (CampaignName=="Moriond18") DowngradeHightPt = true;
    if (BTagger=="CSVv2" && OP=="Loose") DowngradeHightPt = true;
    //if (BTagger=="DeepCSV" && OP=="Tight") DowngradeHightPt = true;
    if (BTagger=="DeepCSV") DowngradeHightPt = true;

    if (DowngradeHightPt) {
      //xPt_fit[nBinsCampaign-1] = xPt_fit[nBinsCampaign-1] + exPt_fit[nBinsCampaign-1];
      //exPt_fit[nBinsCampaign-1] *= 2.;
      sf_fit[nBinsCampaign-1] = sf_fit[nBinsCampaign-2];
      sf_error_fit[nBinsCampaign-1] *= 2.;
      //if (CampaignName!="Moriond18B" && CampaignName!="Moriond18CDE" && OP=="Loose" && MeasType=="mujets") sf_error_fit[nBinsCampaign-2] *= 2.;
    }

  } else if (CampaignName.Contains("2016Legacy")) {

    if (OP=="Loose" && !(BTagger=="DeepJet" && MeasType=="comb")) 
      sf_error_fit[nBinsCampaign-1] *= 1.5;

  } else if (CampaignName.Contains("DeepFlavour2017")) {

    if (OP!="Loose" && MeasType=="mujets") {
      sf_error_fit[nBinsCampaign-1] *= 2.5;
    }

  } else if (CampaignName.Contains("Prompt18")) {

    if (BTagger=="DeepCSV" && OP=="Loose"  && MeasType=="mujets")
      sf_error_fit[nBinsCampaign-1] *= 1.;
    else if (BTagger=="DeepJet" && OP=="Medium" && MeasType=="mujets")
      sf_error_fit[nBinsCampaign-1] *= 20.;
    else if (BTagger=="DeepCSV" && OP=="Medium" && MeasType=="comb") {
      sf_fit[nBinsCampaign-1] = sf_fit[nBinsCampaign-2];
      sf_error_fit[nBinsCampaign-1] *= 200.;
      sf_error_fit[nBinsCampaign-2] *= 200.;
    } else if (MeasType=="comb") {
      sf_error_fit[nBinsCampaign-1] *= 2.;
    } else      
      sf_error_fit[nBinsCampaign-1] *= 200.;
    
  }

  // TGraph for pt dependence fit
  TGraphErrors* sff = new TGraphErrors(nBinsCampaign,xPt,sf_fit,exPt_fit,sf_error_fit);

  // Do we have sample dependence?
  if (AddSampleDependence) {

    for (int bpt = 0; bpt<nBinsCampaign; bpt++)
      sf_error[bpt] = sqrt(sf_error[bpt]*sf_error[bpt] + SampleDependence*SampleDependence);

  } else if (NormalizedChi2>1.) {

    for (int bpt = 0; bpt<nBinsCampaign; bpt++)
      sf_error[bpt] *= sqrt(NormalizedChi2);

  }

  // @@@ Study the dependence on jet pt of the combined scale factors

  // Define fitting function
  TString FittingFunction;
  if (CampaignName=="Winter13") {
    if ( title == "CSVL") FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)";
    else if ( title == "CSVM" || title == "CSVT") 
      FittingFunction = "[0]+[1]*x+[2]*x*x";
  } else if (CampaignName=="7TeVLegacy") {
    if ( title == "CSVL") FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)";
    else if ( title == "CSVM" ) 
      FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)";
    else if ( title == "CSVT" ) 
      FittingFunction = "[0]+[1]*log(x+19)*log(x+18)*(3-[2]*log(x+18))";
    else if ( title == "TCHPT" ) 
      FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)";
  } else if (CampaignName=="Run2015B") {
    FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)";
  } else if (CampaignName=="Run201525ns") {
    if (title.Contains("CSVv2") || title!="JPL")
      FittingFunction = "[0]+[1]*log(x+[3])*log(x+[3])*(3-[2]*log(x+[3]))";
    else 
      FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)"; 
  } else if (CampaignName=="Run201576X") {
    FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)"; 
  } else if (CampaignName=="Run201680X4invfb") {
    if (OP!="Medium" || BTagger=="cMVAv2") FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)"; 
    else FittingFunction = "[0]"; 
  } else if (CampaignName=="Run201680X7invfb") {
    if (MeasType=="mujets") 
      FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)"; 
    else 
      FittingFunction = "[0]+[1]*x";
  } else if (CampaignName=="Moriond17") {
    FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)";
    if (OP=="Loose" && BTagger=="cMVAv2" && MeasType=="ttbar")
      FittingFunction = "[0]+[1]*x"; 
    //if (OP=="Tight" && BTagger=="CSVv2" && MeasType=="mujets")
      //FittingFunction = "[0]+[1]*x";
  /*} else if (CampaignName=="Moriond18B") {
    FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)";
    if (OP=="Medium") FittingFunction = "[0]+[1]*log(x+19)*log(x+18)*(3-[2]*log(x+18))";
    if (BTagger=="CSVv2") FittingFunction = "[0]+[1]*log(x+19)*log(x+18)*(3-[2]*log(x+18))"; 
  } else if (CampaignName=="Moriond18CDE") {
    FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)";
    if (OP=="Loose" || (OP=="Medium" && MeasType=="comb")) FittingFunction = "[0]+[1]*log(x+19)*log(x+18)*(3-[2]*log(x+18))";
    if (OP=="Medium" && BTagger=="CSVv2") FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)";
  } else if (CampaignName=="Moriond18EF") {
    FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)";
    if (OP=="Medium" && MeasType=="comb") FittingFunction = "[0]+[1]*log(x+19)*log(x+18)*(3-[2]*log(x+18))";
    if (OP=="Loose") FittingFunction = "[0]+[1]*log(x+19)*log(x+18)*(3-[2]*log(x+18))";
    if (BTagger=="CSVv2") FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)";*/
  } else if (CampaignName.Contains("Moriond18")) {
    FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)";
    //if (OP=="Loose" && MeasType=="mujets") FittingFunction = "[0]+[1]*log(x+19)*log(x+18)*(3-[2]*log(x+18))";
    //if (OP!="Tight") FittingFunction = "[0]+[1]*log(x+19)*log(x+18)*(3-[2]*log(x+18))";
    if (OP=="Loose") FittingFunction = "[0]+[1]*log(x+19)*log(x+18)*(3-[2]*log(x+18))";
  } else if (CampaignName.Contains("DeepFlavour2017")) {
    if (MeasType=="mujets") {
      FittingFunction = "[0]+[1]*log(x+19)*log(x+18)*(3-[2]*log(x+18))";
    } else {
      FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)";
    }
    //if (OP=="Loose" && MeasType=="mujets") FittingFunction = "[0]+[1]*log(x+19)*log(x+18)*(3-[2]*log(x+18))";
    //if (OP!="Tight") FittingFunction = "[0]+[1]*log(x+19)*log(x+18)*(3-[2]*log(x+18))";
    //if (OP=="Loose") FittingFunction = "[0]+[1]*log(x+19)*log(x+18)*(3-[2]*log(x+18))";
  } else if (CampaignName.Contains("Prompt18")) {
    FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)";
    if (MeasType=="comb")
      FittingFunction = "[0]+[1]*log(x+19)*log(x+18)*(3-[2]*log(x+18))";
  } else if (CampaignName.Contains("2016Legacy")) {
    FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)";
  } else if (CampaignName.Contains("2018Ultimate")) {
    //FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)";
    //if (MeasType=="comb")
      FittingFunction = "[0]+[1]*log(x+19)*log(x+18)*(3-[2]*log(x+18))";
  } else if (CampaignName.Contains("UL17")) {
    FittingFunction = "[0]+[1]*log(x+19)*log(x+18)*(3-[2]*log(x+18))"; 
  } else if (CampaignName.Contains("UL18")) {
    FittingFunction = "[0]+[1]*log(x+19)*log(x+18)*(3-[2]*log(x+18))";
  } else if (CampaignName.Contains("UL16APV")) {
    if (BTagger=="DeepCSV" && (OP=="Tight" || (OP=="Medium" && MeasType=="comb"))) {    
      FittingFunction = "[0]+[1]*log(x+19)*log(x+18)*(3-[2]*log(x+18))";
    } else {//if (OP=="Loose" || BTagger=="DeepJet") {
      FittingFunction = "[0]*(1.+[1]*x)/(1.+[2]*x)";
    }
  } else if (CampaignName.Contains("UL16")) {
    FittingFunction = "[0]+[1]*log(x+19)*log(x+18)*(3-[2]*log(x+18))";
  }

  TF1* fun1 = new TF1("fun1", FittingFunction, MinPtCampaign, MaxPtCampaign);

  // Initialise some parameter in the fitting function
  float MaxPtFit = MaxPtCampaign, MinPtFit = MinPtCampaign;
  if (CampaignName=="Winter13") {
    if ( title == "CSVM" || title == "CSVT")
      fun1->SetParameters(0.938887,0.00017124,-2.76366e-07);
  } else if (CampaignName=="7TeVLegacy") {
    if ( title == "CSVM") 
      fun1->SetParameters(0.899, 0.0565, 0.0606);
    else if ( title == "CSVT")
      fun1->SetParameters(0.9, 0.07/log(670)/log(670), 2./log(670));    
  } else if (CampaignName=="Run2015B") {
    fun1->SetParameters(0.853, 0.0527, 0.453);
    if (OP!="Loose") {
      fun1->FixParameter(1, 0.);
      fun1->FixParameter(2, 0.);
    }
  } else if (CampaignName=="Run201576X") {
    fun1->SetParameter(0, 0.916293);
    fun1->SetParameter(1, 0.0118628);
    fun1->SetParameter(2, 0.0108104);
  } else if (CampaignName=="Run201680X4invfb") {
    if (OP!="Medium") {
      fun1->SetParameter(0, 0.916293);
      fun1->SetParameter(1, 0.0118628);
      fun1->SetParameter(2, 0.0108104);
    } 
  } else if (CampaignName=="Run201680X7invfb") {
    if (MeasType=="mujets") {
      fun1->SetParameter(0,  0.92934);
      fun1->SetParameter(1, -0.00205);
      fun1->SetParameter(2,  0.00208);
    } 
  } else if (CampaignName=="Moriond17") {//.0.846867*((1.+(0.0190562*x))/(1.+(0.0157319
    if (OP=="Tight" && BTagger=="CSVv2" && MeasType=="mujets") {
      fun1->SetParameter(0,  0.846867);
      fun1->SetParameter(1,  0.0190562);
      fun1->SetParameter(2,  0.0157319);
    }
  } else if (CampaignName.Contains("Moriond18")) {
    if (OP=="Loose") {
      fun1->SetParameter(0,  1.1378);
      fun1->SetParameter(1,  0.00599383);
      fun1->SetParameter(2,  0.326135);
    } 		
    if (OP=="Medium") {
      fun1->SetParameter(0,  1.09079);
      fun1->SetParameter(1,  0.180763);
      fun1->SetParameter(2,  0.216798);
    } 	
    if (OP=="Tight") {
      fun1->SetParameter(0,  0.931126);
      fun1->SetParameter(1,  0.00108871);
      fun1->SetParameter(2,  0.00103758);
    } 
  } else if (CampaignName.Contains("Prompt18")) {
    if (BTagger=="DeepCSV" && OP=="Loose" && MeasType=="comb") {
      fun1->SetParameter(0,  0.976211); 
      fun1->SetParameter(1, -0.000122163);
      fun1->SetParameter(2, -4.2984e-05);
    }
    if (BTagger=="DeepCSV" && OP=="Medium" && MeasType=="comb") {
      fun1->SetParameter(0,  0.909339); 
      fun1->SetParameter(1,  0.00354);
      fun1->SetParameter(2,  0.471623);
    }
  } else if (CampaignName.Contains("2018Ultimate")) {
    if (BTagger=="DeepCSV" && OP=="Tight" && MeasType=="comb") {
      fun1->SetParameter(0,  0.828145); 
      fun1->SetParameter(1,  0.00612709);
      fun1->SetParameter(2,  0.39278);
    }
  } /*else if (CampaignName.Contains("UL16APV")) {
    if (BTagger=="DeepCSV" && OP=="Medium" && MeasType=="comb") {
      fun1->SetParameter(0,  1.02033);
      fun1->SetParameter(1,  0.00044452);
      fun1->SetParameter(2,  0.00070354);
    }
  }*/

  // Fit the pt dependence of the combination result
  fun1->SetLineColor(kBlack);
  fun1->SetLineWidth(2);
  fun1->SetLineStyle(1);
  sff->Fit("fun1","0","",MinPtFit,MaxPtFit);
  sff->Fit("fun1","rve0","",MinPtFit,MaxPtFit);
  //$$  sff->Fit("fun1","rvee0"); 
  fun1->Draw("same"); 

  // build fun2 from binned scale factors;
  TString BinnedScaleFactors = "";
  for (int bpt = FirstBinCampaignForFit; bpt<=LastBinCampaignForFit; bpt++) {

    BinnedScaleFactors += "(x>=";
    BinnedScaleFactors += xPt[bpt] - exPt[bpt];
    BinnedScaleFactors += " && x<=";
    BinnedScaleFactors += xPt[bpt] + exPt[bpt];
    BinnedScaleFactors += ") ? ";
    BinnedScaleFactors += sf[bpt];
    BinnedScaleFactors += " : ";
  }
  BinnedScaleFactors += " 0. ";

  // @@@ Fill the final function value and uncertainties
  float fit_val[NBINS], fit_err[NBINS];

  // Compute the function value for the pt bins
  std::cout << std::endl;
  for (int bpt = 0; bpt<nBinsCampaign; bpt++) {

    fit_val[bpt] = fun1->Integral(xPt[bpt]-1.,xPt[bpt]+1.) / 2.; // Run1 keep
    //fit_val[bpt] = sf[bpt];

    //fit_err[bpt] = fun1->IntegralError(xPt[bpt]-exPt[bpt],xPt[bpt]+exPt[bpt]) / exPt[bpt];
    fit_err[bpt] = sf_error[bpt]*fit_val[bpt]/sf[bpt]; // Run1 keep

    if (StoreFittedScaleFactors) {

      fun_val[bpt] = fit_val[bpt];
      fun_err[bpt] = fit_err[bpt];
      
    } else {

      fun_val[bpt] = sf[bpt];
      fun_err[bpt] = sf_error[bpt];

    }

    if (SystematicBreakdown || PtCorrelationBreakdown || StatisticBreakdown || CategoryBreakdown)
      for (int uu = 0; uu<nUncertainties; uu++)	 
	fun_unc[bpt][uu] = sf_uncerbreak[bpt][uu]*fun_val[bpt]/sf[bpt];

  }

  // @@@ Print the combination results

  cout << "Print fit results" << endl << endl; 
  cout << "  Normalized Chi2 = " << NormalizedChi2 << endl << endl;
  TSfun1 = fun1->GetExpFormula("p");
  std::cout << " Tagger: " << title << " within " << MinPtCampaign 
	    << " < pt < " << MaxPtCampaign << " GeV, abs(eta) < 2.5, x = pt" << std::endl;
  std::cout << "  SFb = " << TSfun1 << ";" << std::endl << std::endl;
  std::cout << " Fit details" << std::endl;
  //for (int bpt = FirstBinCampaignForFit; bpt<=LastBinCampaignForFit; bpt++) 
  //std::cout << "   " << xPt[bpt]-exPt[bpt] << "-" << xPt[bpt]+exPt[bpt] << ": SF(bin) = " << sf[bpt] << " +/- " << sf_error[bpt] << "; SF(fit) = " << fit_val[bpt] << " +/- " << fit_err[bpt] << " (" << 100.*fit_err[bpt]/fit_val[bpt] << "%); Stat: " << 100.*sf_stat[bpt]/sf[bpt] << "% (" << 100.*sf_stat[bpt]/sf_error[bpt] << "%)" << std::endl;
  std::cout << " SFb errors" << std::endl;
  for (int bpt = FirstBinCampaignForFit; bpt<=LastBinCampaignForFit; bpt++) 
  std::cout << "  " << sf_error[bpt] << ", " << std::endl;
  std::cout << std::endl;
  
  cout << " Comparison with average ttbar-based measurements" << endl;
   float ScaleFactorTTbarFit = 0., ScaleFactorTTbarFitError = 0.;
   float ScaleFactorTTbarBin = 0., ScaleFactorTTbarBinError = 0.;
   for (int bpt = FirstBinCampaignForFit; bpt<=LastBinCampaignForFit; bpt++) {
     ScaleFactorTTbarFit += fit_val[bpt]*TTbarPt[bpt]/TTbarSpectrumScale; 
     ScaleFactorTTbarFitError += fit_err[bpt]*TTbarPt[bpt]/TTbarSpectrumScale;
     ScaleFactorTTbarBin += sf[bpt]*TTbarPt[bpt]/TTbarSpectrumScale; 
     ScaleFactorTTbarBinError += sf_error[bpt]*TTbarPt[bpt]/TTbarSpectrumScale;
   }
   std::cout << "  ttbar-like SF (fit ) = " << ScaleFactorTTbarFit << " +/- " << ScaleFactorTTbarFitError << std::endl;
   std::cout << "  ttbar-like SF (bins) = " << ScaleFactorTTbarBin << " +/- " << ScaleFactorTTbarBinError << std::endl << std::endl;

   // Take final function for storing scale factors
   if (!StoreFittedScaleFactors) TSfun1 = BinnedScaleFactors;

   // @@@ Producing the plots

   if (PrintPNG || PrintPDF || PrintC || PrintRoot) {
     
     // set the plot style
     TCanvas *c1; TPad *pad1, *pad2;
       
     gROOT->ProcessLine( "gErrorIgnoreLevel = 1001;");
     
     /* Kirill's Style
	static TStyle* plotStyle = PlotStyle(); 
	gROOT->SetStyle("PLOT");   
	gROOT->ForceStyle(); 
	
	pad1->Draw();
	pad2->Draw();
     */// End Kirill's Style 
       
     // Fill vector with results from previous campaign in case a comparison is wished
     float SFb_Comp[nBinsCampaign], SFb_Comp_error[nBinsCampaign];
     
     if (PrintComparison) {
  
       for (int bpt = 0; bpt<nBinsCampaign; bpt++) {
	 
	 double xx = xPt[bpt];
	 
	 if ( title == "CSVv2L" ) {
	   SFb_Comp[bpt]      = funSFb_Comp_CSVv2L(xx);
	   SFb_Comp_error[bpt] = SFb_Comp_error_CSVv2L[bpt];
	 }
	 if ( title == "CSVv2M" ) { 
	   SFb_Comp[bpt]      = funSFb_Comp_CSVv2M(xx);
	   SFb_Comp_error[bpt] = SFb_Comp_error_CSVv2M[bpt];
	 }
	 if ( title == "CSVv2T" ) { 
	   SFb_Comp[bpt]      = funSFb_Comp_CSVv2T(xx);
	   SFb_Comp_error[bpt] = SFb_Comp_error_CSVv2T[bpt];
	 }
	 if ( title == "DeepCSVL" ) {
	   SFb_Comp[bpt]      = funSFb_Comp_DeepCSVL(xx);
	   SFb_Comp_error[bpt] = SFb_Comp_error_DeepCSVL[bpt];
	 }
	 if ( title == "DeepCSVM" ) { 
	   SFb_Comp[bpt]      = funSFb_Comp_DeepCSVM(xx);
	   SFb_Comp_error[bpt] = SFb_Comp_error_DeepCSVM[bpt];
	 }
	 if ( title == "DeepCSVT" ) { 
	   SFb_Comp[bpt]      = funSFb_Comp_DeepCSVT(xx);
	   SFb_Comp_error[bpt] = SFb_Comp_error_DeepCSVT[bpt];
	 }
	 
       }

     }
     
     TGraphErrors* sf_Comp = new TGraphErrors(nBinsCampaign,xPt,SFb_Comp,exPt,SFb_Comp_error);

     TString PlotHeader = title;
     PlotHeader.ReplaceAll("L", " L"); 
     PlotHeader.ReplaceAll("M", " M"); 
     PlotHeader.ReplaceAll("T", " T"); 
     
     if (!SampleCombinationComparison && PrintPtFit) { // @@@ Two canvas plot for single measurement combinations

       // Build the canvas
     
       c1 = new TCanvas("c1", "plots",200,0,700,700);
       c1->SetFillColor(10);
       c1->SetFillStyle(4000);
       c1->SetBorderSize(2);
     
       pad1 = new TPad("pad1","This is pad1",0.02,0.52,0.98,0.98,21);
       pad2 = new TPad("pad2","This is pad2",0.02,0.03,0.98,0.49,21);
          
       // Run2015B Setting
       gStyle->SetOptFit(0);
       gStyle->SetOptStat(0);
       gStyle->SetOptTitle(0);
       c1->Range(0,0,1,1);
       c1->SetFillColor(10);
       c1->SetBorderMode(0);
       c1->SetBorderSize(2);
       c1->SetTickx(1);
       c1->SetTicky(1);
       c1->SetLeftMargin(0.16);
       c1->SetRightMargin(0.02);
       c1->SetTopMargin(0.05);
       c1->SetBottomMargin(0.13);
       c1->SetFrameFillColor(0);
       c1->SetFrameFillStyle(0);
       c1->SetFrameBorderMode(0);
       
       pad1->SetFillColor(0);
       pad1->SetBorderMode(0);
       pad1->SetBorderSize(2);
       //pad1->SetLogy();
       pad1->SetLogx();
       pad1->SetTickx(1);
       pad1->SetTicky(1);
       pad1->SetLeftMargin(0.16);
       pad1->SetRightMargin(0.02);
       pad1->SetTopMargin(0.065);
       pad1->SetBottomMargin(0.13);
       pad1->SetFrameFillStyle(0);
       pad1->SetFrameBorderMode(0);
       pad1->SetFrameFillStyle(0);
       pad1->SetFrameBorderMode(0);
       pad1->Draw();
       
       pad2->SetFillColor(0);
       pad2->SetBorderMode(0);
       pad2->SetBorderSize(2);
       //pad2->SetGridy();
       pad2->SetLogx();
       pad2->SetTickx(1);
       pad2->SetTicky(1);
       pad2->SetLeftMargin(0.16);
       pad2->SetRightMargin(0.02);
       //pad2->SetTopMargin(0.05);
       //pad2->SetBottomMargin(0.31);
       pad2->SetTopMargin(0.065);
       pad2->SetBottomMargin(0.13);
       pad2->SetFrameFillStyle(0);
       pad2->SetFrameBorderMode(0);
       pad2->SetFrameFillStyle(0);
       pad2->SetFrameBorderMode(0);
       pad2->Draw();
       // End Run2015B Setting 
       
       // @@@ Plot the measurements in the top pad
      
       pad1->cd();
       
       float PlotMaxPtCampaign = MaxPtCampaign;
       //if (MeasType=="ttbar") PlotMaxPtCampaign = 300.;
       float YAxisWidth = 0.4;
       //if (OP=="Tight") YAxisWidth = 0.4;
       cout <<MinPtCampaign<<" "<<PlotMaxPtCampaign << endl;
       TH2F* histo = new TH2F("histo","",58,MinPtCampaign,PlotMaxPtCampaign,100,1.-YAxisWidth,1.+YAxisWidth);
       
       histo->Draw(""); 
       histo->SetLabelSize(0.05, "XYZ");
       histo->SetTitleSize(0.06, "XYZ"); 
       histo->SetLabelFont(42, "XYZ"); 
       histo->SetTitleFont(42, "XYZ");
       //histo->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
       histo->GetXaxis()->SetTitle("p_{T} [GeV]");
       //histo->GetYaxis()->SetTitle("Data/Simulation SF_{b}");
       histo->GetYaxis()->SetTitle("SF_{b}");
       //histo->SetTitleOffset(1.1,"X"); // Ideal for .png
       histo->SetTitleOffset(0.95,"X");
       histo->SetTitleOffset(0.8,"Y");
       histo->SetTickLength(0.06,"X");
       histo->SetNdivisions(509, "XYZ");
       histo->GetXaxis()->SetMoreLogLabels();
       histo->GetXaxis()->SetNoExponent();
       
       TGraphErrors *ScaleFactorStatistic[nMeasurements];
       TGraphErrors *ScaleFactorError[nMeasurements];
               
       int nLegendMeasurements = 1; // Counting from 1 to include the combination
       for (int nm = 0; nm<nMeasurements; nm++) {

	 ScaleFactorStatistic[nm] = MakeTGraphErrors(GraphXBins[nm], xpt[nm], MeasuredScaleFactorValue[nm], exstat, MeasuredScaleFactorStatistic[nm], GraphStyle[nm], GraphColor[nm], GraphMarker[nm], GraphSize[nm], GraphWidth[nm]);
	 
	 ScaleFactorError[nm] = MakeTGraphErrors(GraphXBins[nm], xpt[nm], MeasuredScaleFactorValue[nm], expt[nm], MeasuredScaleFactorError[nm], GraphStyle[nm], GraphColor[nm], GraphMarker[nm], GraphSize[nm], 2);
	 
	 if (UseThisMeasurement[TaggingAlgorithm][nm] && PlotTheMeasurement[TaggingAlgorithm][nm] && 
	     (!PrintOnlyTypeMeasurements || TypeMeasurement[nm]==MeasType || MeasType=="comb")) { 
	   
	   ScaleFactorStatistic[nm]->Draw("P"); 
	   ScaleFactorError[nm]->Draw("P");  
	   nLegendMeasurements++;
	   
	 }
	 
       }
       
       // Run2015B Style
       TLatex *tex = new TLatex(0.2,0.88,"CMS"); 
       tex->SetNDC(); 
       tex->SetTextAlign(13);
       tex->SetTextFont(61);
       tex->SetTextSize(0.07475);
       tex->SetLineWidth(2); 
       tex->Draw();                                                                                       
       TLatex *tex2 = new TLatex(0.2,0.79,"Preliminary"); 
       tex2->SetNDC();
       tex2->SetTextAlign(13);
       tex2->SetTextFont(52);
       tex2->SetTextSize(0.05681);
       tex2->SetLineWidth(2);   
       tex2->Draw();   
       
       TLatex *text1 = new TLatex(0.98,0.95125, CampaignLuminosity + CenterOfMassEnergy); 
       text1->SetNDC();                                              
       text1->SetTextAlign(31);                          
       text1->SetTextFont(42);    
       text1->SetTextSize(0.04875);   
       text1->SetLineWidth(2);    
       text1->Draw(); 
       // End Run2015B Style

       // Print the combination result in the top pad
       TGraphErrors* sf0 = new TGraphErrors(nBinsCampaign,xPt,sf,exPt,sf_error);
       bool drawHatchedBand = true;
       if (drawHatchedBand) {
         sf0->SetFillStyle(3005); 
         sf0->SetFillColor(kGray+3);
       } else {
         sf0->SetLineStyle(2);    
         sf0->SetFillStyle(0);      
       }
       sf0->Draw("e2"); // to plot fit band

       // Add legend
       TLegend *leg = new TLegend(0.48,0.64,0.70,0.89);
       leg->SetBorderSize(0);
       leg->SetFillColor(kWhite);
       leg->SetTextFont(62);
       leg->SetHeader(PlotHeader); 
       leg->SetNColumns((nLegendMeasurements-1)/3+1);
       int nLegendEntries = 0;
       //int CombinationPosition = nLegendMeasurements-1;
       int CombinationPosition = leg->GetNColumns()*(nLegendMeasurements/leg->GetNColumns()) - 1;
       for (int nm = 0; nm<nMeasurements; nm++) 
	 if (UseThisMeasurement[TaggingAlgorithm][nm] && PlotTheMeasurement[TaggingAlgorithm][nm] && 
	     (!PrintOnlyTypeMeasurements || TypeMeasurement[nm]==MeasType || MeasType=="comb")) {
	   TString LegText = MeasurementName[nm];
           LegText.ReplaceAll("System8", "System-8");
	   leg->AddEntry(ScaleFactorError[nm], LegText, "PL");
	   nLegendEntries++;
	   if (nLegendEntries==CombinationPosition) leg->AddEntry(sf0,"weighted average","PF");
	 }
       leg->SetY1(0.89 - 0.25*leg->GetNRows()/4);
       if (leg->GetNColumns()>1) leg->SetX2(0.92);
       if (leg->GetNColumns()>2) { leg->SetX1(0.36); leg->SetX2(0.95); }
       leg->Draw();
       
       // Superimpose the single measurements on the combination and legend   
       for (int nm = 0; nm<nMeasurements; nm++)
	 if (UseThisMeasurement[TaggingAlgorithm][MeasurementOrder[nm]] && PlotTheMeasurement[TaggingAlgorithm][MeasurementOrder[nm]] && 
	     (!PrintOnlyTypeMeasurements || TypeMeasurement[MeasurementOrder[nm]]==MeasType || MeasType=="comb")) {
	   ScaleFactorStatistic[MeasurementOrder[nm]]->Draw("P");
	   ScaleFactorError[MeasurementOrder[nm]]->Draw("P");
	 }
       
       // @@@ Now plotting the combination in the bottom pad
       
       pad2->cd();
       
       histo = new TH2F("histo","",58,MinPtCampaign,PlotMaxPtCampaign,100,1.-YAxisWidth,1.+YAxisWidth);
       
       histo->Draw("");
       histo->SetLabelSize(0.05, "XYZ");
       histo->SetTitleSize(0.06, "XYZ"); 
       histo->SetLabelFont(42, "XYZ"); 
       histo->SetTitleFont(42, "XYZ");
       //histo->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
       histo->GetXaxis()->SetTitle("p_{T} [GeV]");
       //histo->GetYaxis()->SetTitle("Data/Simulation SF_{b}");
       histo->GetYaxis()->SetTitle("SF_{b}");
       //histo->SetTitleOffset(1.1,"X"); // Ideal for .png
       histo->SetTitleOffset(0.95,"X"); // 
       histo->SetTitleOffset(0.8,"Y");
       histo->SetTickLength(0.06,"X");
       histo->SetNdivisions(509, "XYZ");
       histo->GetXaxis()->SetMoreLogLabels();
       histo->GetXaxis()->SetNoExponent();

       if (PrintComparison) {
     
	 sf_Comp->SetFillColor(kYellow-9);
	 sf_Comp->Draw("e2"); 
 
       }

       // Superimpose the combination result
       //sf0->SetFillStyle(3005);
       //sf0->SetFillColor(kBlack);
       sf0->Draw("e2"); 
       
       // Add the functions values to the bottom pad
       TGraphErrors* fun0 = new TGraphErrors(16,xPt,fit_val,exPt,fit_err);
       fun0->SetMarkerColor(kRed);
       fun0->SetLineColor(kRed);
       fun0->SetLineStyle(1);
       fun0->SetMarkerStyle(24);
       fun0->SetMarkerSize(0.001);
       fun0->SetLineWidth(2);
       fun0->Draw("P"); 
       fun1->SetMarkerColor(kBlack);
       fun1->Draw("same"); 
       
       // Add legend
       leg = new TLegend(0.48,0.64,0.70,0.89);
       leg->SetBorderSize(0);
       leg->SetFillColor(kWhite);
       leg->SetTextFont(62);
       leg->SetTextSize(0.05);   
       leg->SetHeader(PlotHeader);
       leg->AddEntry(sf0,"weighted average","PF");
       leg->AddEntry(fun1,"fit","L");
       leg->AddEntry(fun0,"fit #pm (stat #oplus syst)","LE");
       if (PrintComparison) leg->AddEntry(sf_Comp,LegComp,"F");
       leg->Draw();
       
       /* // Run1 Style
	  tex = new TLatex(0.17,1,"CMS Preliminary, " + CampaignLuminosity + "  at #sqrt{s} = " + CenterOfMassEnergy);
	  tex->SetNDC();
	  tex->SetTextAlign(13);
	  tex->SetLineWidth(2);
	  tex->Draw();
       */ // End Run1 Style
       
       // Run2015B Style
       tex = new TLatex(0.2,0.88,"CMS"); 
       tex->SetNDC(); 
       tex->SetTextAlign(13);
       tex->SetTextFont(61);
       tex->SetTextSize(0.07475);
       tex->SetLineWidth(2); 
       tex->Draw();                                                                                       
       tex2 = new TLatex(0.2,0.79,"Preliminary"); 
       tex2->SetNDC();
       tex2->SetTextAlign(13);
       tex2->SetTextFont(52);
       tex2->SetTextSize(0.05681);
       tex2->SetLineWidth(2);   
       tex2->Draw();      
       
       text1 = new TLatex(0.98,0.95125, CampaignLuminosity + CenterOfMassEnergy); 
       text1->SetNDC();                                              
       text1->SetTextAlign(31);                          
       text1->SetTextFont(42);    
       text1->SetTextSize(0.04875);   
       text1->SetLineWidth(2);    
       text1->Draw();  
       // End Run2015B Style
   
     } else if (!PrintPtFit) { // @@@ One canvas plot for measurements' comaparison

       // Build the canvas
     
       c1 = new TCanvas("c1", "plots",200,0,700,350);
       c1->SetFillColor(10);
       c1->SetFillStyle(4000);
       c1->SetBorderSize(2);
     
       pad1 = new TPad("pad1","This is pad1",0.02,0.03,0.98,0.98,21);
       
       // Run2015B Setting
       gStyle->SetOptFit(0);
       gStyle->SetOptStat(0);
       gStyle->SetOptTitle(0);
       c1->Range(0,0,1,1);
       c1->SetFillColor(10);
       c1->SetBorderMode(0);
       c1->SetBorderSize(2);
       c1->SetTickx(1);
       c1->SetTicky(1);
       c1->SetLeftMargin(0.16);
       c1->SetRightMargin(0.02);
       c1->SetTopMargin(0.05);
       c1->SetBottomMargin(0.13);
       c1->SetFrameFillColor(0);
       c1->SetFrameFillStyle(0);
       c1->SetFrameBorderMode(0);
       
       pad1->SetFillColor(0);
       pad1->SetBorderMode(0);
       pad1->SetBorderSize(2);
       //pad1->SetLogy();
       pad1->SetLogx();
       pad1->SetTickx(1);
       pad1->SetTicky(1);
       pad1->SetLeftMargin(0.16);
       pad1->SetRightMargin(0.02);
       pad1->SetTopMargin(0.065);
       pad1->SetBottomMargin(0.13);
       pad1->SetFrameFillStyle(0);
       pad1->SetFrameBorderMode(0);
       pad1->SetFrameFillStyle(0);
       pad1->SetFrameBorderMode(0);
       pad1->Draw();

       // End Run2015B Setting 
       
       // @@@ Plot the measurements in the pad

       pad1->cd();
  
       float PlotMaxPtCampaign = MaxPtCampaign;
       //if (MeasType=="ttbar") PlotMaxPtCampaign = 300.;
       float YAxisWidth = 0.4;
       if (OP=="Tight") YAxisWidth = 0.4;
       TH2F* histo = new TH2F("histo","",58,MinPtCampaign,PlotMaxPtCampaign,100,1.-YAxisWidth,1.+YAxisWidth);
       
       histo->Draw(""); 
       histo->SetLabelSize(0.05, "XYZ");
       histo->SetTitleSize(0.06, "XYZ"); 
       histo->SetLabelFont(42, "XYZ"); 
       histo->SetTitleFont(42, "XYZ");
       //histo->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
       histo->GetXaxis()->SetTitle("p_{T} [GeV]");
       //histo->GetYaxis()->SetTitle("Data/Simulation SF_{b}");
       histo->GetYaxis()->SetTitle("SF_{b}");
       //histo->SetTitleOffset(1.1,"X"); // Ideal for .png
       histo->SetTitleOffset(0.95,"X");
       histo->SetTitleOffset(0.8,"Y");
       histo->SetTickLength(0.06,"X");
       histo->SetNdivisions(509, "XYZ");
       histo->GetXaxis()->SetMoreLogLabels();
       histo->GetXaxis()->SetNoExponent();
       
       TGraphErrors *ScaleFactorStatistic[nMeasurements];
       TGraphErrors *ScaleFactorError[nMeasurements];
               
       int nLegendMeasurements = 1; // Counting from 1 to include the combination
       for (int nm = 0; nm<nMeasurements; nm++) {

	 ScaleFactorStatistic[nm] = MakeTGraphErrors(GraphXBins[nm], xpt[nm], MeasuredScaleFactorValue[nm], exstat, MeasuredScaleFactorStatistic[nm], GraphStyle[nm], GraphColor[nm], GraphMarker[nm], GraphSize[nm], GraphWidth[nm]);
	 
	 ScaleFactorError[nm] = MakeTGraphErrors(GraphXBins[nm], xpt[nm], MeasuredScaleFactorValue[nm], expt[nm], MeasuredScaleFactorError[nm], GraphStyle[nm], GraphColor[nm], GraphMarker[nm], GraphSize[nm], 2);
	 
	 if (UseThisMeasurement[TaggingAlgorithm][nm] && PlotTheMeasurement[TaggingAlgorithm][nm] && 
	     (!PrintOnlyTypeMeasurements || TypeMeasurement[nm]==MeasType || MeasType=="comb")) { 
	   
	   ScaleFactorStatistic[nm]->Draw("P"); 
	   ScaleFactorError[nm]->Draw("P");  
	   nLegendMeasurements++;
	   
	 }
	 
       }
       
       // Run2015B Style
       TLatex *tex = new TLatex(0.2,0.88,"CMS"); 
       tex->SetNDC(); 
       tex->SetTextAlign(13);
       tex->SetTextFont(61);
       tex->SetTextSize(0.07475);
       tex->SetLineWidth(2); 
       tex->Draw();                                                                                       
       TLatex *tex2 = new TLatex(0.2,0.79,"Preliminary"); 
       tex2->SetNDC();
       tex2->SetTextAlign(13);
       tex2->SetTextFont(52);
       tex2->SetTextSize(0.05681);
       tex2->SetLineWidth(2);   
       tex2->Draw();   
       
       TLatex *text1 = new TLatex(0.98,0.95125, CampaignLuminosity + CenterOfMassEnergy); 
       text1->SetNDC();                                              
       text1->SetTextAlign(31);                          
       text1->SetTextFont(42);    
       text1->SetTextSize(0.04875);   
       text1->SetLineWidth(2);    
       text1->Draw(); 
       // End Run2015B Style

       // Print the combination result in the top pad
       TGraphErrors* sf0 = new TGraphErrors(nBinsCampaign,xPt,sf,exPt,sf_error);
       sf0->SetFillStyle(3005);
       sf0->SetFillColor(kGray+3);
       sf0->Draw("e2"); // to plot fit band

       // Add legend
       TLegend *leg = new TLegend(0.48,0.64,0.70,0.89);
       leg->SetBorderSize(0);
       leg->SetFillColor(kWhite);
       leg->SetTextFont(62);
       leg->SetHeader(PlotHeader); 
       leg->SetNColumns((nLegendMeasurements-1)/3+1);
       int nLegendEntries = 0;
       //int CombinationPosition = nLegendMeasurements-1;
       int CombinationPosition = leg->GetNColumns()*(nLegendMeasurements/leg->GetNColumns()) - 1;
       for (int nm = 0; nm<nMeasurements; nm++) 
	 if (UseThisMeasurement[TaggingAlgorithm][nm] && PlotTheMeasurement[TaggingAlgorithm][nm] && 
	     (!PrintOnlyTypeMeasurements || TypeMeasurement[nm]==MeasType || MeasType=="comb")) {
	   TString LegText = MeasurementName[nm];
           LegText.ReplaceAll("System8", "System-8");
	   leg->AddEntry(ScaleFactorError[nm], LegText, "PL");
	   nLegendEntries++;
	   if (nLegendEntries==CombinationPosition) leg->AddEntry(sf0,"weighted average","PF");
	 }
       leg->SetY1(0.89 - 0.25*leg->GetNRows()/4);
       if (leg->GetNColumns()>1) leg->SetX2(0.92);
       if (leg->GetNColumns()>2) { leg->SetX1(0.36); leg->SetX2(0.95); }
       leg->Draw();
       
       // Superimpose the single measurements on the combination and legend   
       for (int nm = 0; nm<nMeasurements; nm++)
	 if (UseThisMeasurement[TaggingAlgorithm][MeasurementOrder[nm]] && PlotTheMeasurement[TaggingAlgorithm][MeasurementOrder[nm]] && 
	     (!PrintOnlyTypeMeasurements || TypeMeasurement[MeasurementOrder[nm]]==MeasType || MeasType=="comb")) {
	   ScaleFactorStatistic[MeasurementOrder[nm]]->Draw("P");
	   ScaleFactorError[MeasurementOrder[nm]]->Draw("P");
	 }

     } else { // @@@ One canvas plot for combinations' comaparison

       // Build the canvas
     
       c1 = new TCanvas("c1", "plots",200,0,700,350);
       c1->SetFillColor(10);
       c1->SetFillStyle(4000);
       c1->SetBorderSize(2);
     
       pad1 = new TPad("pad1","This is pad1",0.02,0.03,0.98,0.98,21);
       
       // Run2015B Setting
       gStyle->SetOptFit(0);
       gStyle->SetOptStat(0);
       gStyle->SetOptTitle(0);
       c1->Range(0,0,1,1);
       c1->SetFillColor(10);
       c1->SetBorderMode(0);
       c1->SetBorderSize(2);
       c1->SetTickx(1);
       c1->SetTicky(1);
       c1->SetLeftMargin(0.16);
       c1->SetRightMargin(0.02);
       c1->SetTopMargin(0.05);
       c1->SetBottomMargin(0.13);
       c1->SetFrameFillColor(0);
       c1->SetFrameFillStyle(0);
       c1->SetFrameBorderMode(0);
       
       pad1->SetFillColor(0);
       pad1->SetBorderMode(0);
       pad1->SetBorderSize(2);
       //pad1->SetLogy();
       pad1->SetLogx();
       pad1->SetTickx(1);
       pad1->SetTicky(1);
       pad1->SetLeftMargin(0.16);
       pad1->SetRightMargin(0.02);
       pad1->SetTopMargin(0.065);
       pad1->SetBottomMargin(0.13);
       pad1->SetFrameFillStyle(0);
       pad1->SetFrameBorderMode(0);
       pad1->SetFrameFillStyle(0);
       pad1->SetFrameBorderMode(0);
       pad1->Draw();

       // End Run2015B Setting 

       // Plot the measurements 

       pad1->cd();
  
       float PlotMaxPtCampaign = MaxPtCampaign;
       //if (MeasType=="ttbar") PlotMaxPtCampaign = 300.;
       float YAxisWidth = 0.3;
       if (OP=="Tight") YAxisWidth = 0.4;
       TH2F* histo = new TH2F("histo","",58,MinPtCampaign,PlotMaxPtCampaign,100,1.-YAxisWidth,1.+YAxisWidth);
       
       histo->Draw(""); 
       histo->SetLabelSize(0.05, "XYZ");
       histo->SetTitleSize(0.06, "XYZ"); 
       histo->SetLabelFont(42, "XYZ"); 
       histo->SetTitleFont(42, "XYZ");
       //histo->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
       histo->GetXaxis()->SetTitle("p_{T} [GeV]");
       //histo->GetYaxis()->SetTitle("Data/Simulation SF_{b}");
       histo->GetYaxis()->SetTitle("SF_{b}");
       //histo->SetTitleOffset(1.1,"X"); // Ideal for .png
       histo->SetTitleOffset(0.95,"X");
       histo->SetTitleOffset(0.8,"Y");
       histo->SetTickLength(0.06,"X");
       histo->SetNdivisions(509, "XYZ");
       histo->GetXaxis()->SetMoreLogLabels();
       histo->GetXaxis()->SetNoExponent();
       
       TGraphErrors *ScaleFactorStatistic[nMeasurements];
       TGraphErrors *ScaleFactorError[nMeasurements];
       
       for (int nm = 0; nm<nMeasurements; nm++) {

	 if (UseThisMeasurement[TaggingAlgorithm][nm] && (MeasurementName[nm]=="mujets" || MeasurementName[nm]=="ttbar")) { 

	   ScaleFactorStatistic[nm] = MakeTGraphErrors(GraphXBins[nm], xpt[nm], MeasuredScaleFactorValue[nm], exstat, MeasuredScaleFactorStatistic[nm], GraphStyle[nm], GraphColor[nm], GraphMarker[nm], GraphSize[nm], GraphWidth[nm]);
	   
	   ScaleFactorError[nm] = MakeTGraphErrors(GraphXBins[nm], xpt[nm], MeasuredScaleFactorValue[nm], expt[nm], MeasuredScaleFactorError[nm], GraphStyle[nm], GraphColor[nm], GraphMarker[nm], GraphSize[nm], 2);
	   
	   ScaleFactorStatistic[nm]->Draw("P"); 
	   ScaleFactorError[nm]->Draw("P"); 
	   
	 }
	 
       }

       // Superimpose the combination result
       //sf0->SetFillStyle(3005);
       //sf0->SetFillColor(kBlack);
       //sf0->Draw("e2"); 
       
       // Add the functions values
       TGraphErrors* fun0 = new TGraphErrors(16,xPt,fit_val,exPt,fit_err);
       fun0->SetMarkerColor(kBlack);
       fun0->SetLineColor(kBlack);
       fun0->SetLineStyle(1);
       fun0->SetMarkerStyle(24);
       fun0->SetMarkerSize(0.001);
       fun0->SetLineWidth(2);
       //fun0->Draw("P"); 
       fun0->SetFillColor(30);
       fun0->Draw("e2"); 
       
       // Plotting the comparison with other measurements if wished 
       if (PrintComparison) {

	 sf_Comp->SetFillStyle(3005);
	 //sf_Comp->SetFillColor(30);
	 sf_Comp->SetLineStyle(1);
	 sf_Comp->SetLineWidth(2);
	 sf_Comp->Draw("e2");
 
       }

       // Add the fitted function
       fun1->SetMarkerColor(kBlack);
       fun1->Draw("same"); 
       
       // Superimpose the single measurements on the combination
       for (int nm = 0; nm<nMeasurements; nm++)
	 if (UseThisMeasurement[TaggingAlgorithm][MeasurementOrder[nm]] && 
	     (MeasurementName[MeasurementOrder[nm]]=="mujets" || MeasurementName[MeasurementOrder[nm]]=="ttbar")) {
	   ScaleFactorStatistic[MeasurementOrder[nm]]->Draw("P");
	   ScaleFactorError[MeasurementOrder[nm]]->Draw("P"); 
	 }
       
       // Add legend
       TLegend *leg = new TLegend(0.48,0.64,0.70,0.89);
       leg->SetBorderSize(0);
       leg->SetFillColor(kWhite);
       leg->SetTextFont(62);
       leg->SetTextSize(0.05);   
       leg->SetHeader(PlotHeader);
       for (int nm = 0; nm<nMeasurements; nm++) 
	 if (UseThisMeasurement[TaggingAlgorithm][nm] && (MeasurementName[nm]=="mujets" || MeasurementName[nm]=="ttbar")) {
	   TString LegText;
           if (MeasurementName[nm]=="mujets") LegText = "muon jets";
           if (MeasurementName[nm]=="ttbar")  LegText = "t#bar{t}";
	   leg->AddEntry(ScaleFactorError[nm], LegText, "PL");
	 }
       //leg->AddEntry(sf0,"weighted average","PF");
       //leg->AddEntry(fun1,"comb","L");
       leg->AddEntry(fun0,"comb #pm (stat #oplus syst)","LF");
       if (PrintComparison) leg->AddEntry(sf_Comp,LegComp,"LF");
       leg->Draw();
       
       // Run2015B Style
       TLatex *tex = new TLatex(0.2,0.88,"CMS"); 
       tex->SetNDC(); 
       tex->SetTextAlign(13);
       tex->SetTextFont(61);
       tex->SetTextSize(0.07475);
       tex->SetLineWidth(2); 
       tex->Draw();                                                                                       
       TLatex *tex2 = new TLatex(0.2,0.79,"Preliminary"); 
       tex2->SetNDC();
       tex2->SetTextAlign(13);
       tex2->SetTextFont(52);
       tex2->SetTextSize(0.05681);
       tex2->SetLineWidth(2);   
       tex2->Draw();      
       
       TLatex *text1 = new TLatex(0.98,0.95125, CampaignLuminosity + CenterOfMassEnergy); 
       text1->SetNDC();                                              
       text1->SetTextAlign(31);                          
       text1->SetTextFont(42);    
       text1->SetTextSize(0.04875);   
       text1->SetLineWidth(2);    
       text1->Draw();  
       // End Run2015B Style

     }
       
     // @@@ Save the plots

     c1->Update();
     
     TString PlotName = "./Plots/" + CampaignName + "/SFb_" + CampaignName + "_" + title;
     for (int nm = 0; nm<nMeasurements; nm++) {
       if (UseThisMeasurement[TaggingAlgorithm][nm]) {
	 TString MeasString = MeasurementPlot[nm];
	 if (VetoedMeasurement[TaggingAlgorithm][nm] || (MeasType!="comb" && TypeMeasurement[nm]!=MeasType))
	   MeasString += "nofit";
	 if ((!PlotTheMeasurement[TaggingAlgorithm][nm] && !SampleCombinationComparison) || 
	     (!SampleCombinationComparison && PrintOnlyTypeMeasurements && TypeMeasurement[nm]!=MeasType && MeasType!="comb") ||
	     (SampleCombinationComparison && MeasurementName[nm]!="mujets" && MeasurementName[nm]!="ttbar") )
	   MeasString += "noplot";
	 if (!MeasString.Contains("nofitnoplot")) PlotName += MeasString;
       }
     }
     if (PrintComparison) PlotName += CompTitle;
     if (!PrintPtFit) PlotName += "_noPtFit";
     if (PrintPNG)  c1->Print(PlotName + ".png");
     if (PrintPDF)  c1->Print(PlotName + ".pdf");
     if (PrintC)    c1->Print(PlotName + ".C");
     if (PrintRoot) c1->Print(PlotName + ".root");
   
   }

}

void StoreScaleFactorCombination(string BTagger, TString OP = "All", string MeasType = "comb", int JetFlavour = 5, TString StoreSystematicBreakdown = "NONE", TString CentralValues = "Function") {

  BTagCalibration calib(BTagger);

  const int nOperatingPoints = 3;
  string OperatingPoint[nOperatingPoints] = {"Loose", "Medium", "Tight"};
  
  const int nMeasurementTypes = 3;
  string MeasTypeName[nMeasurementTypes] = {"comb", "mujets", "ttbar"};

  if (StoreSystematicBreakdown=="Systematic") SystematicBreakdown = true;
  if (StoreSystematicBreakdown=="Category") CategoryBreakdown = true;
  if (StoreSystematicBreakdown=="PtCorrelation") PtCorrelationBreakdown = true;
  if (StoreSystematicBreakdown=="Statistic") StatisticBreakdown = true;

  float PtBinEdge[nBinsCampaign+1]; PtBinEdge[0] = xPt[0] - exPt[0];
  for (int bpt = 0; bpt<nBinsCampaign; bpt++) PtBinEdge[bpt+1] = xPt[bpt] + exPt[bpt];
  TH1F *UncertaintyHisto = new TH1F ("UncertaintyHisto", "", nBinsCampaign, PtBinEdge);  
  float DiscMin = 0., DiscMax = 1.;
  if (BTagger=="JP")   { DiscMin =  0.; DiscMax = 5.; }
  if (BTagger=="cMVA") { DiscMin = -1.; DiscMax = 1.; }

  for (int mt = 0; mt<nMeasurementTypes; mt++) 
    if((MeasType=="All" && MeasTypeName[mt]!="ttbar") || MeasType==MeasTypeName[mt]) {

      SetAlgorithmToUse(BTagger);

      for (int op = 0; op<nOperatingPoints; op++) 
	if (OP=="All" || OP.Contains(OperatingPoint[op])) {
		  
	  BTagEntry::OperatingPoint wp = BTagEntry::OP_LOOSE;
	  if (OperatingPoint[op]=="Loose")  wp = BTagEntry::OP_LOOSE;
	  if (OperatingPoint[op]=="Medium") wp = BTagEntry::OP_MEDIUM;
	  if (OperatingPoint[op]=="Tight")  wp = BTagEntry::OP_TIGHT;

	  if (CampaignName=="Run201680X4invfb" && OperatingPoint[op]=="Medium")
	    for (int nm = 0; nm<nMeasurements; nm++) 
	      if (MeasurementName[nm]=="LT") 
		VetoedMeasurement[TaggingAlgorithm][nm] = true;
	  
	  ScaleFactorCombination(BTagger, OperatingPoint[op], MeasTypeName[mt]);
	  
	  int maxUncertainty = nUncertainties;
	  if (!SystematicBreakdown && !PtCorrelationBreakdown && !CategoryBreakdown) {
	    maxUncertainty = 1;
	    if (StatisticBreakdown) maxUncertainty = 2;
	  }

	  for (int ThisFlavour = 4; ThisFlavour<=5; ThisFlavour++) {

	    if (JetFlavour==4 && ThisFlavour==5) continue;
	    if (JetFlavour==5 && ThisFlavour==4) continue;

	    float SysFactor = (ThisFlavour==5 || JetFlavour==4) ? 1. : cJetsInflationFactor[op];
	    
	    BTagEntry::JetFlavor jfl;
	    if (ThisFlavour==4) jfl = BTagEntry::FLAV_C;
	    if (ThisFlavour==5) jfl = BTagEntry::FLAV_B;
	    
	    double TotalUncertainty[nBinsCampaign];
	    double NotBrokenUncertainty[nBinsCampaign];
  
	    for (int uu = 0; uu<=maxUncertainty; uu++) {
	      if (uu<=1 || (SystematicBreakdown && IsForBreakdown[uu-1]) || (PtCorrelationBreakdown && IsPtCorrelated[uu-1]) || 
		  (StatisticBreakdown && UncertaintyName[uu-1]=="_statistic") || (CategoryBreakdown && IsForCategoryBreakdown[uu-1])) {

		if (uu>=2 && !isUsedSystematic[uu-1]) continue;

		for (int vs = 0; vs<2; vs++) {
		  
		  if (uu==0 && vs==1) continue;
		  
		  string SysFlag = "central";
		  
		  TString ThisFunction = TSfun1;
		  
		  if (uu>0) {
		    
		    if (vs==0) { SysFlag = "up"; ThisFunction += " + "; }
		    else if (vs==1) { SysFlag = "down";  ThisFunction += " - "; }
		    
		    UncertaintyHisto->Reset();
		    
		    if (uu==1) {
		      
		      for (int bpt = FirstBinCampaignForFit; bpt<=LastBinCampaignForFit; bpt++) {
			
			UncertaintyHisto->SetBinContent(bpt+1, SysFactor*fun_err[bpt]);
			if (vs==0) {
			  NotBrokenUncertainty[bpt] = pow(SysFactor*fun_err[bpt], 2);
			  TotalUncertainty[bpt] = pow(SysFactor*fun_err[bpt], 2);
			}

		      }		      
		      
		    } else { 
    
                      if (UncertaintyName[uu-1]=="_isrDef") 
                        SysFlag += "_isr";
                      else if (UncertaintyName[uu-1]=="_fsrDef")
                        SysFlag += "_fsr";
                      else 
		        SysFlag += UncertaintyName[uu-1];

		      for (int bpt = FirstBinCampaignForFit; bpt<=LastBinCampaignForFit; bpt++) {
			
			UncertaintyHisto->SetBinContent(bpt+1, SysFactor*fun_unc[bpt][uu-1]);
			if (vs==0) NotBrokenUncertainty[bpt] -= pow(SysFactor*fun_unc[bpt][uu-1], 2);

		      }
		      
		    }
		    
		  }
		  
		  /* This is not working yet ...		  
		     BTagEntry::Parameters params(wp, MeasTypeName[mt], SysFlag, jfl, -2.4, 2.4, MinPtCampaign, MaxPtCampaign, 0, 1);
		     
		     if (uu>0) {
		     
		     BTagEntry e2(UncertaintyHisto, params); 
		     ThisFunction += e2.formula;
		     cout << e2.formula << endl;
		     
		     }
		     cout << ThisFunction << endl;
		     const TF1 SFFun("SFFun", ThisFunction, MinPtCampaign, MaxPtCampaign);
		     BTagEntry e1(&SFFun, params);  
		     
		     calib.addEntry(e1);
		  */
		  
		  // ... so doing like this for the time being
		  if (uu==0) {
		    
		    BTagEntry::Parameters params(wp, MeasTypeName[mt], SysFlag, jfl, -2.5, 2.5, MinPtCampaign, MaxPtCampaign, DiscMin, DiscMax);

		    const TF1 SFFun("SFFun", ThisFunction, MinPtCampaign, MaxPtCampaign);
		    BTagEntry e1(&SFFun, params);  
		    
		    calib.addEntry(e1);
		    
		  } else if (uu==1 || SystematicBreakdown || CategoryBreakdown || StatisticBreakdown) {
		    
		    for (int bpt = FirstBinCampaignForFit; bpt<=LastBinCampaignForFit; bpt++) {
		      
		      BTagEntry::Parameters params(wp, MeasTypeName[mt], SysFlag, jfl, -2.5, 2.5, PtBinEdge[bpt], PtBinEdge[bpt+1], DiscMin, DiscMax);
		        
		      TString SystFunction = ThisFunction;
		      SystFunction += UncertaintyHisto->GetBinContent(bpt+1);
		      const TF1 SFFun("SFFun", SystFunction, PtBinEdge[bpt], PtBinEdge[bpt+1]);
		      BTagEntry e1(&SFFun, params);  
		      
		      calib.addEntry(e1);
		      
		    }
		    
		  }
		  
		}

	      }
	    }
	    
	    if (SystematicBreakdown || CategoryBreakdown) {
	      
	      for (int vs = 0; vs<2; vs++) {

		TString ThisFunction = TSfun1; string SysFlag;
		if (vs==0) { ThisFunction += " + "; SysFlag = "up_"; }
		else if (vs==1) { ThisFunction += " - "; SysFlag = "down_"; }
		SysFlag += (!CategoryBreakdown) ? "notinbreakdown" : "type3";
		
		for (int bpt = FirstBinCampaignForFit; bpt<=LastBinCampaignForFit; bpt++) {
		  
		  BTagEntry::Parameters params(wp, MeasTypeName[mt], SysFlag, jfl, -2.5, 2.5, PtBinEdge[bpt], PtBinEdge[bpt+1], DiscMin, DiscMax);
		
		  TString SystFunction = ThisFunction;
		  SystFunction += sqrt(NotBrokenUncertainty[bpt]);
		  
		  const TF1 SFFun("SFFun", SystFunction, PtBinEdge[bpt], PtBinEdge[bpt+1]);
		  BTagEntry e1(&SFFun, params);  
		  
		  calib.addEntry(e1);
		  
		}

	      }
	    
	    }

	    if (PtCorrelationBreakdown) {
	      
	      for (int vs = 0; vs<2; vs++) {

		TString ThisFunction = TSfun1; string SysFlag;
		if (vs==0) { ThisFunction += " + "; SysFlag = "up_"; }
		else if (vs==1) { ThisFunction += " - "; SysFlag = "down_"; }
		SysFlag += "ptuncorrelated";

		for (int bpt = FirstBinCampaignForFit; bpt<=LastBinCampaignForFit; bpt++) {
		
		  BTagEntry::Parameters params(wp, MeasTypeName[mt], SysFlag, jfl, -2.5, 2.5, PtBinEdge[bpt], PtBinEdge[bpt+1], DiscMin, DiscMax);

		  if (vs==0 && AddSampleDependence) 
		    NotBrokenUncertainty[bpt] -= pow(SampleDependence, 2);

		  TString SystFunction = ThisFunction;
		  SystFunction += sqrt(NotBrokenUncertainty[bpt]);
		  const TF1 SFFun("SFFun", SystFunction, PtBinEdge[bpt], PtBinEdge[bpt+1]);
		  BTagEntry e1(&SFFun, params);  
		  
		  calib.addEntry(e1);
		  
		}

		if (vs==0) SysFlag = "up_"; 
		else if (vs==1) SysFlag = "down_"; 
		SysFlag += "ptcorrelated";
		
		for (int bpt = FirstBinCampaignForFit; bpt<=LastBinCampaignForFit; bpt++) {
		
		  BTagEntry::Parameters params(wp, MeasTypeName[mt], SysFlag, jfl, -2.5, 2.5, PtBinEdge[bpt], PtBinEdge[bpt+1], DiscMin, DiscMax);
		  
		  TString SystFunction = ThisFunction;
		  SystFunction += sqrt(TotalUncertainty[bpt] - NotBrokenUncertainty[bpt]);
		  const TF1 SFFun("SFFun", SystFunction, PtBinEdge[bpt], PtBinEdge[bpt+1]);
		  BTagEntry e1(&SFFun, params);  
		  
		  calib.addEntry(e1);
		  
		}
	      
	      }
	    
	    }

	  }
	 
	  if (CampaignName=="Run201680X4invfb")
	    for (int nm = 0; nm<nMeasurements; nm++) 
	      if (MeasurementName[nm]=="LT") 
		VetoedMeasurement[TaggingAlgorithm][nm] = false; 
	  
	}
	  
    }
  
  TString CSVFileName = BTagger;
  if (OP=="Loose")  CSVFileName += "L"; 
  if (OP=="Medium") CSVFileName += "M"; 
  if (OP=="Tight")  CSVFileName += "T";
  CSVFileName += "_" + CampaignNameString;
  if (MeasType!="All") { CSVFileName += "_"; CSVFileName += MeasType; } 
  if (CategoryBreakdown) CSVFileName += "_CategoryBreakdown";
  else if (SystematicBreakdown) CSVFileName += "_SystematicBreakdown";
  if (PtCorrelationBreakdown) CSVFileName += "_PtCorrelationBreakdown";
  if (StatisticBreakdown) CSVFileName += "_StatisticBreakdown";
  CSVFileName.ReplaceAll("_SystematicBreakdown_StatisticBreakdown", "_UncertaintyBreakdown");
  std::ofstream outFile(CSVFileName + ".csv");
  calib.makeCSV(outFile);
  outFile.close();
  
}

/*
void MergeSystematicCategories(string BTagger, string InputCSVFileName, string MeasType) {

  BTagCalibration inputcalib(BTagger, InputCSVFileName);

  BTagCalibration outputcalib(BTagger);

  nBinsCampaignForFit = 6;

  MinPtCampaign = xPt[0] - exPt[0];
  MaxPtCampaign = xPt[nBinsCampaignForFit-1] + exPt[nBinsCampaignForFit-1];

  float DiscMin = 0., DiscMax = 1.;
  if (BTagger=="JP")   { DiscMin =  0.; DiscMax = 5.; }
  if (BTagger=="cMVA") { DiscMin = -1.; DiscMax = 1.; }

  const int nUncorrelatedSources = 5;
 string UncorrelatedSourceName[nUncorrelatedSources] = {"statistic", "backg", "met", "hpp", "topmass"};
  bool OnlyFromTTbar[nUncorrelatedSources] = {false, true, true, true, true};
  
  for (int op = 0; op<3; op++) {

    BTagEntry::OperatingPoint WP;
    if (op==0) WP = BTagEntry::OP_LOOSE;
    if (op==1) WP = BTagEntry::OP_MEDIUM;
    if (op==2) WP = BTagEntry::OP_TIGHT;

    BTagCalibrationReader *CentralReader = new BTagCalibrationReader(&inputcalib, WP, MeasType, "central");

    for (int ThisFlavour = 4; ThisFlavour<=5; ThisFlavour++) {
      
      BTagEntry::JetFlavor jfl;
      if (ThisFlavour==4) jfl = BTagEntry::FLAV_C;
      if (ThisFlavour==5) jfl = BTagEntry::FLAV_B;
  
      TF1 CentralFunction = CentralReader->Func(jfl, 0., (MinPtCampaign+MaxPtCampaign)/2.);
      TString CentralFormula = CentralFunction.GetExpFormula();
      
      BTagEntry::Parameters params(WP, MeasType, "central", jfl, -2.4, 2.4, MinPtCampaign, MaxPtCampaign, DiscMin, DiscMax);
      
      const TF1 SFFun("SFFun", CentralFormula, MinPtCampaign, MaxPtCampaign);
      BTagEntry e1(&SFFun, params);  
    
      outputcalib.addEntry(e1);
      
      for (int bpt = 0; bpt<nBinsCampaignForFit; bpt++) {

	// total uncertainty
	BTagCalibrationReader *UncertaintyReader = new BTagCalibrationReader(&inputcalib, WP, MeasType, "up");

	double CentralSF = CentralReader->eval(jfl, 0., xPt[bpt]);
	double UpSF = UncertaintyReader->eval(jfl, 0., xPt[bpt]);
	double TotalUncertainty = UpSF - CentralSF;

	// up
	BTagEntry::Parameters paramsUp(WP, MeasType, "up", jfl, -2.4, 2.4, xPt[bpt]-exPt[bpt], xPt[bpt]+exPt[bpt], DiscMin, DiscMax);
	
	TString UpFormula = CentralFormula;
	UpFormula += "+"; UpFormula += TotalUncertainty;
	const TF1 SFUpFun("SFUpFun", UpFormula, xPt[bpt]-exPt[bpt], xPt[bpt]+exPt[bpt]);
	BTagEntry eup(&SFUpFun, paramsUp);  
	
	outputcalib.addEntry(eup);

	// down
	BTagEntry::Parameters paramsDown(WP, MeasType, "down", jfl, -2.4, 2.4, xPt[bpt]-exPt[bpt], xPt[bpt]+exPt[bpt], DiscMin, DiscMax);
	
	TString DownFormula = CentralFormula;
	DownFormula += "-"; DownFormula += TotalUncertainty;
	const TF1 SFDownFun("SFDownFun", DownFormula, xPt[bpt]-exPt[bpt], xPt[bpt]+exPt[bpt]);
	BTagEntry edown(&SFDownFun, paramsDown);  
	
	outputcalib.addEntry(edown);

	// uncorrelated uncertainty
	double UncorrelatedUncertainty = 0;
	for (int us = 0; us<nUncorrelatedSources; us++) {
	 
	  if (OnlyFromTTbar[us] && MeasType=="mujets") continue;
	  
	  string ThisSource = "up_"; ThisSource += UncorrelatedSourceName[us];
	  BTagCalibrationReader *ThisReader = new BTagCalibrationReader(&inputcalib, WP, MeasType, ThisSource);
	  
	  double ThisSF = ThisReader->eval(jfl, 0., xPt[bpt]);
	  double ThisUncertainty = ThisSF - CentralSF;
	  UncorrelatedUncertainty += ThisUncertainty*ThisUncertainty;
	  
	}
	UncorrelatedUncertainty = sqrt(UncorrelatedUncertainty);

	// up uncorrelated
	BTagEntry::Parameters paramsUpUnc(WP, MeasType, "up_uncorrelated", jfl, -2.4, 2.4, xPt[bpt]-exPt[bpt], xPt[bpt]+exPt[bpt], DiscMin, DiscMax);
	
	TString UpUncFormula = CentralFormula;
	UpUncFormula += "+"; UpUncFormula += UncorrelatedUncertainty;
	const TF1 SFUpUncFun("SFUpUncFun", UpUncFormula, xPt[bpt]-exPt[bpt], xPt[bpt]+exPt[bpt]);
	BTagEntry eupunc(&SFUpUncFun, paramsUpUnc);  
	
	outputcalib.addEntry(eupunc);

	// down uncorrelated
	BTagEntry::Parameters paramsDownUnc(WP, MeasType, "down_uncorrelated", jfl, -2.4, 2.4, xPt[bpt]-exPt[bpt], xPt[bpt]+exPt[bpt], DiscMin, DiscMax);
	
	TString DownUncFormula = CentralFormula;
	DownUncFormula += "-"; DownUncFormula += UncorrelatedUncertainty;
	const TF1 SFDownUncFun("SFDownUncFun", DownUncFormula, xPt[bpt]-exPt[bpt], xPt[bpt]+exPt[bpt]);
	BTagEntry edownunc(&SFDownUncFun, paramsDownUnc);  
	
	outputcalib.addEntry(edownunc);

	// correlated uncertainty
	double CorrelatedUncertainty = sqrt(TotalUncertainty*TotalUncertainty - UncorrelatedUncertainty*UncorrelatedUncertainty);

	// up correlated
	BTagEntry::Parameters paramsUpCorr(WP, MeasType, "up_correlated", jfl, -2.4, 2.4, xPt[bpt]-exPt[bpt], xPt[bpt]+exPt[bpt], DiscMin, DiscMax);
	
	TString UpCorrFormula = CentralFormula;
	UpCorrFormula += "+"; UpCorrFormula += CorrelatedUncertainty;
	const TF1 SFUpCorrFun("SFUpCorrFun", UpCorrFormula, xPt[bpt]-exPt[bpt], xPt[bpt]+exPt[bpt]);
	BTagEntry eupcorr(&SFUpCorrFun, paramsUpCorr);  
	
	outputcalib.addEntry(eupcorr);

	// down correlated
	BTagEntry::Parameters paramsDownCorr(WP, MeasType, "down_correlated", jfl, -2.4, 2.4, xPt[bpt]-exPt[bpt], xPt[bpt]+exPt[bpt], DiscMin, DiscMax);
	
	TString DownCorrFormula = CentralFormula;
	DownCorrFormula += "-"; DownCorrFormula += CorrelatedUncertainty;
	const TF1 SFDownCorrFun("SFDownCorrFun", DownCorrFormula, xPt[bpt]-exPt[bpt], xPt[bpt]+exPt[bpt]);
	BTagEntry edowncorr(&SFDownCorrFun, paramsDownCorr);  
	
	outputcalib.addEntry(edowncorr);
	
      }
    
    }

  }
  
  TString CSVFileName = BTagger;
  CSVFileName += "_";
  CSVFileName += CampaignNameString;
  CSVFileName += "_";
  CSVFileName += MeasType;
  CSVFileName += "_CategoriesBreakdown";
  std::ofstream outFile(CSVFileName + ".csv");
  outputcalib.makeCSV(outFile);
  outFile.close();
  
}
*/
