#!/usr/bin/env python
import os
import sys
import ROOT
import math
import copy
import optparse
from array import *
import CondTools.BTau.dataLoader as dataLoader
from collections import defaultdict
from collections import OrderedDict

if __name__ == '__main__':

    # Input parameters
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)   
    parser.add_option('--campaign'       , dest='campaign'       , help='Measurement campaign'           , default='') 
    parser.add_option('--algorithm'      , dest='algorithm'      , help='Algorithms to be analysed'      , default='all')
    parser.add_option('--combination'    , dest='combination'    , help='Combination to be performed'    , default='all')
    parser.add_option('--workingpoint'   , dest='workingpoint'   , help='Working points to be analysed'  , default='all')
    parser.add_option('--vetomethod'     , dest='vetomethod'     , help='Measurement methods to veto'    , default='NONE')
    parser.add_option('--maskmethod'     , dest='maskmethod'     , help='Measurement methods to mask'    , default='NONE')
    parser.add_option('--plotoff'        , dest='plotoff'        , help='Don\'t make plots'              , default=False, action='store_true')
    parser.add_option('--store'          , dest='store'          , help='Store csv files'                , default=False, action='store_true')
    parser.add_option('--storebybins'    , dest='storebybins'    , help='Store csv files by bins'        , default=False, action='store_true')
    parser.add_option('--breaksyst'      , dest='breaksyst'      , help='Store systematics breakdown'    , default=False, action='store_true')
    parser.add_option('--yearcorr'       , dest='yearcorr'       , help='Store year correlations'        , default=False, action='store_true')
    parser.add_option('--cjetsoff'       , dest='cjetsoff'       , help='Don\'t store SFs for c jets'    , default=False, action='store_true')
    parser.add_option('--forceptfit'     , dest='forceptfit'     , help='Force the pt-dependence fit'    , default=False, action='store_true')
    parser.add_option('--plotdir'        , dest='plotdir'        , help='Output directory for plots'     , default='./Plots')
    parser.add_option('--csvfiledir'     , dest='csvfiledir'     , help='Output directory for csv files' , default='./CSVFiles')
    (opt, args) = parser.parse_args()

    # Read campaign info
    if os.path.exists('./Measurements/'+opt.campaign+'/CampaignInfo.py'):
        handle = open('./Measurements/'+opt.campaign+'/CampaignInfo.py','r')
        exec(handle)
        handle.close() 

    else:
        print 'Campaign', opt.campaign, 'not found. Please, specify a valid campaign name'
        exit()

    # Some settings

    #ROOT.gROOT.ProcessLine('gErrorIgnoreLevel = 1001;')

    if opt.storebybins or opt.breaksyst or opt.yearcorr: opt.store = True

    if opt.store:
        opt.storebyfunction = not opt.storebybins
        jetFlavoursToBeStored = [ ROOT.BTagEntry.FLAV_B ]
        if not opt.cjetsoff: jetFlavoursToBeStored.append(ROOT.BTagEntry.FLAV_C)

    else:
        opt.storebyfunction = False

    opt.doptfit = not opt.plotoff or opt.storebyfunction or opt.forceptfit

    # Get list of combinations, algorithms, working points, and measurement methods
    if opt.algorithm!='all':
       for algo in algorithms.keys():
           if algo.lower() not in opt.algorithm.lower(): del algorithms[algo]

    if opt.combination!='all':
       for comb in copy.deepcopy(combinations):
           if comb.lower() not in opt.combination.lower(): combinations.remove(comb)

    if opt.workingpoint!='all':
       for wp in workingPoints.keys():
           if wp.lower() not in opt.workingpoint.lower(): del workingPoints[wp]

    if opt.vetomethod!='None':
        for method in measurements.keys():
            if method.lower() in opt.vetomethod.lower() and method not in vetoedMethods: vetoedMethods.append(method)

    if opt.maskmethod!='None':
        for method in measurements.keys():
            if method.lower() in opt.maskmethod.lower() and method not in maskedMethods: maskedMethods.append(method) 

    # Loop on algorithms, combinations, and working points
    for algo in algorithms:

        if opt.store:
            
            wpFlag = '' if opt.workingpoint=='all' else '-'.join([ ''.join([ x for x in y if x.isupper() ]) for y in workingPoints ])
            csvFileNameList = [ algo+wpFlag, opt.campaign ]
            if opt.combination!='all': csvFileNameList.insert(1, '-'.join([ x for x in combinations ]))
            if len(maskedMethods)>0: csvFileNameList.insert(csvFileNameList.index(opt.campaign), 'masked'+'-'.join([ measurements[x]['plotname'] for x in set(vetoedMethods+maskedMethods) ]))

            if opt.storebybins: csvFileNameList.append('Binned')
            if opt.breaksyst: csvFileNameList.append('BategoryBreakdown')
            if opt.yearcorr: csvFileNameList.append('YearCorrelation')

            os.system('mkdir -p '+opt.csvfiledir)
            csvFile = ROOT.BTagCalibration(opt.csvfiledir+'/'+'_'.join(csvFileNameList)+'.csv')

        for comb in combinations:
            for wp in workingPoints:

                measuredScaleFactors = OrderedDict() 
                graphScaleFactors, graphScaleFactorsTotal = OrderedDict(), OrderedDict()

                # loop on the measurement methods ...
                for method in measurements:
                    if comb in measurements[method]['data'] and method not in vetoedMethods:

                        measurementFileName = '/'.join([ 'Measurements', opt.campaign, algo+'_'+method+'.csv' ])
                        if not os.path.exists(measurementFileName): 
                            measurementFileName = measurementFileName.replace(algo, algo+workingPoints[wp])
                        if not os.path.exists(measurementFileName): 
                            print 'Measurement file for', method, algo, wp, 'not found'
                            exit()

                        loaders = dataLoader.get_data(measurementFileName)

                        # ... to read the results of scale factor measurements, and fill ...
                        for data in loaders:
                            for e in data.entries:
                                if str(e.params.operatingPoint)==workingPoints[wp] and float(e.formula)>0. and e.params.ptMin<maxPtCampaign:

                                    ptbin = 'Pt-'+str(int(e.params.ptMin))+'to'+str(int(e.params.ptMax))
                                    if ptbin not in measuredScaleFactors: measuredScaleFactors[ptbin] = OrderedDict()
                                    if method not in measuredScaleFactors[ptbin]: measuredScaleFactors[ptbin][method] = {}
                                    measuredScaleFactors[ptbin][method][e.params.sysType] = float(e.formula)

                        # ... the systematics to be used to construct the covariance matrix and ...
                        for ptbin in measuredScaleFactors.keys():
                            if method in  measuredScaleFactors[ptbin]:

                                measuredScaleFactors[ptbin][method]['systematics'] = {}
                                for syst in measuredScaleFactors[ptbin][method].keys():
                                    if syst!='central':
                                    
                                        systName = 'total' if '_' not in syst else syst.split('_')[1]
                                        if systName not in measuredScaleFactors[ptbin][method]:

                                            upSF   = measuredScaleFactors[ptbin][method]['up']   if systName=='total' else measuredScaleFactors[ptbin][method]['up_'+systName]
                                            downSF = measuredScaleFactors[ptbin][method]['down'] if systName=='total' else measuredScaleFactors[ptbin][method]['down_'+systName]
                                            upSyst   = upSF   - measuredScaleFactors[ptbin][method]['central']
                                            downSyst = downSF - measuredScaleFactors[ptbin][method]['central']
                                            measuredScaleFactors[ptbin][method]['systematics'][systName] = math.copysign((abs(upSyst)+abs(downSyst))/2., upSyst)

                                            if systName!='total':

                                                if systName not in systematicPtCorrelated and systName not in systematicPtUncorrelated: 
                                                    print 'Error:', systName, 'not assigned as pt-correlated nor as pt-uncorrelated'
                                                    exit()

                                                if systName in systematicPtCorrelated and systName in systematicPtUncorrelated:
                                                    print 'Error:', systName, 'assigned both as pt-correlated and pt-uncorrelated'
                                                    exit()

                                                if systName not in systematicYearCorrelated and systName not in systematicYearUncorrelated:
                                                    print 'Error:', systName, 'not assigned as year-correlated nor as year-uncorrelated'
                                                    exit()

                                                if systName in systematicYearCorrelated and systName in systematicYearUncorrelated:
                                                    print 'Error:', systName, 'assigned both as year-correlated and year-uncorrelated'
                                                    exit()
   
                                                if systName not in type1Systematics and systName not in type2Systematics and systName not in type3Systematics:
                                                    print 'Error:', systName, 'not assigned to any breakdown category'
                                                    exit()

                                                if (systName in type1Systematics and systName in type2Systematics) or (systName in type1Systematics and systName in type3Systematics) or (systName in type2Systematics and systName in type3Systematics):
                                                    print 'Error:', systName, 'assigned to more than one breakdown category'
                                                    exit()                                          

                        # ... the graph to be used for the final plots
                        graphMetod, graphMetodTotal = ROOT.TGraphErrors(), ROOT.TGraphErrors() 

                        ibin = 0
                        for ptbin in measuredScaleFactors:
                            if method in  measuredScaleFactors[ptbin]:

                                minPt, maxPt = float(ptbin.split('-')[1].split('to')[0]), float(ptbin.split('to')[1]) 
                                midPt = (maxPt+minPt)/2. + measurements[method]['shift']
                      
                                graphMetod.SetPoint(ibin, midPt, measuredScaleFactors[ptbin][method]['central'])
                                graphMetod.SetPointError(ibin, (maxPt-minPt)/2., measuredScaleFactors[ptbin][method]['systematics']['statistic'])
                                graphMetodTotal.SetPoint(ibin, midPt, measuredScaleFactors[ptbin][method]['central'])
                                graphMetodTotal.SetPointError(ibin, (maxPt-minPt)/2., measuredScaleFactors[ptbin][method]['systematics']['total'])

                                ibin += 1

                        graphMetod.SetFillColor(measurements[method]['color']);    graphMetodTotal.SetFillColor(measurements[method]['color']);
                        graphMetod.SetMarkerStyle(measurements[method]['marker']); graphMetodTotal.SetMarkerStyle(measurements[method]['marker']);
                        graphMetod.SetMarkerColor(measurements[method]['color']);  graphMetodTotal.SetMarkerColor(measurements[method]['color']);
                        graphMetod.SetMarkerSize(measurements[method]['size']);    graphMetodTotal.SetMarkerSize(measurements[method]['size']);
                        graphMetod.SetLineStyle(1);                                graphMetodTotal.SetLineStyle(1);
                        graphMetod.SetLineColor(measurements[method]['color']);    graphMetodTotal.SetLineColor(measurements[method]['color']);
                        graphMetod.SetLineWidth(measurements[method]['width']);    graphMetodTotal.SetLineWidth(1);

                        graphScaleFactors[method] = graphMetod
                        graphScaleFactorsTotal[method] = graphMetodTotal

                        # Finally, apply masks for this method
                        for ptbin in measuredScaleFactors.keys():
                            if method in measuredScaleFactors[ptbin]:
                                if method in maskedMethods or (method in maskedMeasurements and ptbin in maskedMeasurements[method]):
                                    del measuredScaleFactors[ptbin][method]
                            if len(measuredScaleFactors[ptbin].keys())==0: del measuredScaleFactors[ptbin]

                # Build the matrices for the fit
                nMeasurements, nCombinedScaleFactors = 0, len(measuredScaleFactors.keys())
                for ptbin in measuredScaleFactors: nMeasurements += len(measuredScaleFactors[ptbin].keys())

                matrixU = ROOT.TMatrixD(nMeasurements, nCombinedScaleFactors)
                matrixU *= 0.

                scaleFactorVector = ROOT.TMatrixD(nMeasurements, 1)
                scaleFactorVector *= 0.

                covarianceMatrix = ROOT.TMatrixD(nMeasurements, nMeasurements)
                covarianceMatrix *= 0.

                breakdownCovarianceMatrices = OrderedDict()

                if opt.breaksyst:
                    breakdownCovarianceMatrices['type1'] = ROOT.TMatrixD(nMeasurements, nMeasurements); breakdownCovarianceMatrices['type1'] *= 0.
                    breakdownCovarianceMatrices['type3'] = ROOT.TMatrixD(nMeasurements, nMeasurements); breakdownCovarianceMatrices['type3'] *= 0.

                if opt.yearcorr:
                    breakdownCovarianceMatrices['correlated']   = ROOT.TMatrixD(nMeasurements, nMeasurements); breakdownCovarianceMatrices['correlated']   *= 0.
                    breakdownCovarianceMatrices['uncorrelated'] = ROOT.TMatrixD(nMeasurements, nMeasurements); breakdownCovarianceMatrices['uncorrelated'] *= 0.

                row = 0
                for iptbin, row_ptbin in enumerate(measuredScaleFactors):
                    for row_method in measuredScaleFactors[row_ptbin]:  

                        matrixU[row][iptbin] = 1.
                        scaleFactorVector[row][0] = measuredScaleFactors[row_ptbin][row_method]['central']

                        column = 0
                        for column_ptbin in measuredScaleFactors:
                            for column_method in measuredScaleFactors[column_ptbin]:

                                for syst in measuredScaleFactors[row_ptbin][row_method]['systematics']:
                                    if syst in measuredScaleFactors[column_ptbin][column_method]['systematics'] and syst!='total':

                                        matrixElement = measuredScaleFactors[row_ptbin][row_method]['systematics'][syst]*measuredScaleFactors[column_ptbin][column_method]['systematics'][syst]

                                        if row_ptbin!=column_ptbin:
                                            if syst in systematicPtUncorrelated: continue
                                            elif syst in ptCorrelationCoefficients: matrixElement *= ptCorrelationCoefficients[syst]
                                        elif row_method!=column_method and syst in statisticalCorrelationCoefficients:
                                            if row_method+'-'+column_method in statisticalCorrelationCoefficients[syst]:
                                                if row_ptbin in statisticalCorrelationCoefficients[syst][row_method+'-'+column_method]:
                                                    matrixElement *= statisticalCorrelationCoefficients[syst][row_method+'-'+column_method][row_ptbin]

                                        covarianceMatrix[row][column] += matrixElement

                                        if opt.breaksyst:
                                            if syst in type1Systematics: breakdownCovarianceMatrices['type1'][row][column] += matrixElement
                                            elif syst in type3Systematics: breakdownCovarianceMatrices['type3'][row][column] += matrixElement
                                            elif syst in type2Systematics: 
                                                if syst not in breakdownCovarianceMatrices:
                                                    breakdownCovarianceMatrices[syst] = ROOT.TMatrixD(nMeasurements, nMeasurements)
                                                    breakdownCovarianceMatrices[syst] *= 0.
                                                breakdownCovarianceMatrices[syst][row][column] += matrixElement

                                        if opt.yearcorr:
                                            if syst in systematicYearCorrelated: breakdownCovarianceMatrices['correlated'][row][column] += matrixElement
                                            if syst in systematicYearUncorrelated: breakdownCovarianceMatrices['uncorrelated'][row][column] += matrixElement

                                column += 1

                        row += 1

                if covarianceMatrix.Determinant()==0.:
                    print 'Covariance matrix not invertible for combination', comb, 'algorithm', algo, 'working point', wp
                    continue

                # Make the fit
                transposedMatrixU = ROOT.TMatrixD(nCombinedScaleFactors, nMeasurements)
                transposedMatrixU.Transpose(matrixU)
 
                inverseCovarianceMatrix = ROOT.TMatrixD(covarianceMatrix)
                inverseCovarianceMatrix.Invert()

                auxiliaryMatrix1 = ROOT.TMatrixD(nCombinedScaleFactors, nMeasurements)
                auxiliaryMatrix1.Mult(transposedMatrixU, inverseCovarianceMatrix) 

                auxiliaryMatrix2 = ROOT.TMatrixD(nCombinedScaleFactors, nCombinedScaleFactors)
                auxiliaryMatrix2.Mult(auxiliaryMatrix1, matrixU)

                if auxiliaryMatrix2.Determinant()==0.:
                   print 'Auxiliary matrix 2 not invertible for combination', comb, 'algorithm', algo, 'working point', wp
                   continue

                auxiliaryMatrix2.Invert()

                coefficientsMatrix = ROOT.TMatrixD(nCombinedScaleFactors, nMeasurements)
                coefficientsMatrix.Mult(auxiliaryMatrix2, auxiliaryMatrix1)

                # Get the scale factors and their uncertainties
                combinedScaleFactorVector = ROOT.TMatrixD(nCombinedScaleFactors, 1)
                combinedScaleFactorVector.Mult(coefficientsMatrix, scaleFactorVector)

                transposedCoefficientsMatrix = ROOT.TMatrixD(nMeasurements, nCombinedScaleFactors)
                transposedCoefficientsMatrix.Transpose(coefficientsMatrix)

                auxiliaryMatrix3 = ROOT.TMatrixD(nCombinedScaleFactors, nMeasurements)
                auxiliaryMatrix3.Mult(coefficientsMatrix, covarianceMatrix)

                scaleFactorUncertaintyMatrix = ROOT.TMatrixD(nCombinedScaleFactors, nCombinedScaleFactors)
                scaleFactorUncertaintyMatrix.Mult(auxiliaryMatrix3, transposedCoefficientsMatrix)
             
                combinedScaleFactorUncertaintyVector = ROOT.TMatrixD(nCombinedScaleFactors, 1)
                for iptbin in range(nCombinedScaleFactors): combinedScaleFactorUncertaintyVector[iptbin][0] = math.sqrt(scaleFactorUncertaintyMatrix[iptbin][iptbin])

                combinedScaleFactorUncertaintyBreakdownVectors = OrderedDict()


                for syst in breakdownCovarianceMatrices:

                    breakdownAuxiliaryMatrix3 = ROOT.TMatrixD(nCombinedScaleFactors, nMeasurements)
                    breakdownAuxiliaryMatrix3.Mult(coefficientsMatrix, breakdownCovarianceMatrices[syst])

                    scaleFactorUncertaintyBreakdownMatrix = ROOT.TMatrixD(nCombinedScaleFactors, nCombinedScaleFactors)
                    scaleFactorUncertaintyBreakdownMatrix.Mult(breakdownAuxiliaryMatrix3, transposedCoefficientsMatrix)

                    combinedScaleFactorUncertaintyBreakdownVector = ROOT.TMatrixD(nCombinedScaleFactors, 1)
                    for iptbin in range(nCombinedScaleFactors): combinedScaleFactorUncertaintyBreakdownVector[iptbin][0] = math.sqrt(scaleFactorUncertaintyBreakdownMatrix[iptbin][iptbin])

                    combinedScaleFactorUncertaintyBreakdownVectors[syst] = combinedScaleFactorUncertaintyBreakdownVector

                # Compute the fit chi2 
                normalizedChi2 = 0.

                for meas0 in range(nMeasurements):
                    for meas1 in range(nMeasurements): 
                        for ptbin0 in range(nCombinedScaleFactors):
                            for ptbin1 in range(nCombinedScaleFactors):
                                normalizedChi2 += matrixU[meas0][ptbin0]*(scaleFactorVector[meas0][0]-combinedScaleFactorVector[ptbin0][0])*inverseCovarianceMatrix[meas0][meas1]*matrixU[meas1][ptbin1]*(scaleFactorVector[meas1][0]-combinedScaleFactorVector[ptbin1][0])

                if nMeasurements>nCombinedScaleFactors: normalizedChi2 /= (nMeasurements - nCombinedScaleFactors)

                # Special error treatments, in case of mis-agreement between the scale factor measurements
                if sampleDependence>0.:

                    for iptbin in range(nCombinedScaleFactors):

                        sampleDependenceUncertaintySquared = pow(sampleDependence*combinedScaleFactorVector[iptbin][0],2)
                        combinedScaleFactorUncertaintyVector[iptbin][0] = math.sqrt(pow(combinedScaleFactorUncertaintyVector[iptbin][0],2)+sampleDependenceUncertaintySquared)
                        for syst in [ 'type3', 'correlated' ]:
                            if syst in combinedScaleFactorUncertaintyBreakdownVectors:
                                combinedScaleFactorUncertaintyBreakdownVectors[syst][iptbin][0] = math.sqrt(pow(combinedScaleFactorUncertaintyBreakdownVectors[syst][iptbin][0],2)+sampleDependenceUncertaintySquared)

                elif normalizedChi2>normalizedChi2Tollerance:

                    for iptbin in range(nCombinedScaleFactors):

                        chi2InflationFactor = math.sqrt(normalizedChi2)
                        combinedScaleFactorUncertaintyVector[iptbin][0] *= chi2InflationFactor
                        for syst in [ 'type1', 'uncorrelated' ]:
                            if syst in combinedScaleFactorUncertaintyBreakdownVectors:
                                chi2InflationSquared = pow(combinedScaleFactorUncertaintyVector[iptbin][0],2)*(1.-1./pow(chi2InflationFactor,2))
                                combinedScaleFactorUncertaintyBreakdownVectors[syst][iptbin][0] = math.sqrt(pow(combinedScaleFactorUncertaintyBreakdownVectors[syst][iptbin][0],2)+chi2InflationSquared)

                # Fit pt-dependence of combined scale factors
                if opt.doptfit:

                    graphCombinedScaleFactors = ROOT.TGraphErrors()
                    minCombinedPt, maxCombinedPt = 999999., -1.

                    for iptbin, ptbin in enumerate(measuredScaleFactors):

                        minPt, maxPt = float(ptbin.split('-')[1].split('to')[0]), float(ptbin.split('to')[1])
                        midPt = (maxPt+minPt)/2.

                        graphCombinedScaleFactors.SetPoint(iptbin, midPt, combinedScaleFactorVector[iptbin][0])
                        graphCombinedScaleFactors.SetPointError(iptbin, (maxPt-minPt)/2., combinedScaleFactorUncertaintyVector[iptbin][0])

                        minCombinedPt = min(minCombinedPt, minPt)
                        maxCombinedPt = max(maxCombinedPt, maxPt)

                    fittingFunction = ROOT.TF1('fittingFunction', str(fittingFunctions[comb][algo][wp].GetExpFormula().ReplaceAll('p','')), minCombinedPt, maxCombinedPt)
                    for par in range(fittingFunction.GetNpar()): fittingFunction.SetParameter(par, fittingFunctions[comb][algo][wp].GetParameter(par)) 
                    fittingFunction.SetLineColor(ROOT.kBlack)
                    fittingFunction.SetLineWidth(2)
                    fittingFunction.SetLineStyle(1)

                    graphCombinedScaleFactors.Fit('fittingFunction', '0',    '', minCombinedPt, maxCombinedPt)
                    graphCombinedScaleFactors.Fit('fittingFunction', 'rve0', '', minCombinedPt, maxCombinedPt)

                # Print results
                print '\nFit performed for combination', comb, 'algorithm', algo, 'working point', wp, 'with normalized Chi2 =', normalizedChi2, '\n'
                for iptbin, ptbin in enumerate(measuredScaleFactors):
                    print '    Combined scale factor for pt bin', ptbin, ':', round(combinedScaleFactorVector[iptbin][0],3), '+-', round(combinedScaleFactorUncertaintyVector[iptbin][0],3)
                print '\n'

                if opt.doptfit:
                    print 'Pt-dependence function:', str(fittingFunction.GetExpFormula('p')), '\n' 
                    for iptbin, ptbin in enumerate(measuredScaleFactors):
                        midPt = (float(ptbin.split('-')[1].split('to')[0])+float(ptbin.split('to')[1]))/2.
                        print '    Fitted scale factor for pt bin', ptbin, ':', round(fittingFunction.Eval(midPt),3), 'difference with combined scale factor:', round(fittingFunction.Eval(midPt)-combinedScaleFactorVector[iptbin][0],3)
                print '\n'

                # Store results for csv files 
                if opt.store:

                    for csvFlavour in jetFlavoursToBeStored:

                        systematicsToBeStored = [ 'up', 'down' ]

                        if opt.storebyfunction:

                            centralScaleFactor = str(fittingFunction.GetExpFormula('p'))
                            params = ROOT.BTagEntry.Parameters(csvWorkingPoints[wp], comb, 'central', csvFlavour, -maxEtaCampaign, maxEtaCampaign, 
                                                               minCombinedPt, maxCombinedPt, algorithms[algo][0], algorithms[algo][1])
                            SFFun = ROOT.TF1 ('SFFun', centralScaleFactor, minCombinedPt, maxCombinedPt)
                            entry = ROOT.BTagEntry(SFFun, params)
                            csvFile.addEntry(entry)

                        else: systematicsToBeStored.insert(0, 'central')

                        for syst in combinedScaleFactorUncertaintyBreakdownVectors:
                            systematicsToBeStored.append('up_'+syst)
                            systematicsToBeStored.append('down_'+syst)

                        for iptbin, ptbin in enumerate(measuredScaleFactors):

                            minPt, maxPt = float(ptbin.split('-')[1].split('to')[0]), float(ptbin.split('to')[1])
                            midPt = (maxPt+minPt)/2.

                            for syst in systematicsToBeStored:

                                if syst=='central': 

                                    centralScaleFactor = str(combinedScaleFactorVector[iptbin][0])
                                    scaleFactorSystematic = ''

                                else:

                                    if '_' not in syst: scaleFactorUncertainty = combinedScaleFactorUncertaintyVector[iptbin][0]
                                    else: scaleFactorUncertainty = combinedScaleFactorUncertaintyBreakdownVectors[syst.split('_')[1]][iptbin][0]
                                    
                                    if opt.storebyfunction:        
                                        fittedScaleFactor = fittingFunction.Integral(midPt-1.,midPt+1.)/2.
                                        scaleFactorUncertainty *= fittedScaleFactor/combinedScaleFactorVector[iptbin][0]

                                    if csvFlavour==ROOT.BTagEntry.FLAV_C:
                                        scaleFactorUncertainty *= cJetsInflationFactor[wp] 

                                    if syst.split('_')[0]=='up': scaleFactorSystematic = '+'+str(scaleFactorUncertainty) if scaleFactorUncertainty>=0. else str(scaleFactorUncertainty)
                                    elif syst.split('_')[0]=='down': scaleFactorSystematic = '-'+str(scaleFactorUncertainty) if scaleFactorUncertainty>=0. else '+'+str(abs(scaleFactorUncertainty))

                                params = ROOT.BTagEntry.Parameters(csvWorkingPoints[wp], comb, syst.replace('type1','statistic'), csvFlavour, 
                                                                   -maxEtaCampaign, maxEtaCampaign, minPt, maxPt, algorithms[algo][0], algorithms[algo][1])
                                SFFun = ROOT.TF1 ('SFFun', centralScaleFactor+scaleFactorSystematic, minPt, maxPt)
                                entry = ROOT.BTagEntry(SFFun, params)
                                csvFile.addEntry(entry)

                # Plot the results of the scale factor combination
                if not opt.plotoff:
     
                    c1 = ROOT.TCanvas('c1','plots',200,0,700,700)
                    c1.SetFillColor(10)
                    c1.SetFillStyle(4000)
                    c1.SetBorderSize(2)
     
                    pad1 = ROOT.TPad('pad1','This is pad1',0.02,0.52,0.98,0.98,21)
                    pad2 = ROOT.TPad('pad2','This is pad2',0.02,0.03,0.98,0.49,21)
          
                    # Run2015B Setting
                    ROOT.gStyle.SetOptFit(0)
                    ROOT.gStyle.SetOptStat(0)
                    ROOT.gStyle.SetOptTitle(0)
                    c1.Range(0,0,1,1)
                    c1.SetFillColor(10)
                    c1.SetBorderMode(0)
                    c1.SetBorderSize(2)
                    c1.SetTickx(1)
                    c1.SetTicky(1)
                    c1.SetLeftMargin(0.16)
                    c1.SetRightMargin(0.02)
                    c1.SetTopMargin(0.05)
                    c1.SetBottomMargin(0.13)
                    c1.SetFrameFillColor(0)
                    c1.SetFrameFillStyle(0)
                    c1.SetFrameBorderMode(0)
       
                    pad1.SetFillColor(0)
                    pad1.SetBorderMode(0)
                    pad1.SetBorderSize(2)
                    #pad1.SetLogy()
                    pad1.SetLogx()
                    pad1.SetTickx(1)
                    pad1.SetTicky(1)
                    pad1.SetLeftMargin(0.16)
                    pad1.SetRightMargin(0.02)
                    pad1.SetTopMargin(0.065)
                    pad1.SetBottomMargin(0.13)
                    pad1.SetFrameFillStyle(0)
                    pad1.SetFrameBorderMode(0)
                    pad1.SetFrameFillStyle(0)
                    pad1.SetFrameBorderMode(0)
                    pad1.Draw()
       
                    pad2.SetFillColor(0)
                    pad2.SetBorderMode(0)
                    pad2.SetBorderSize(2)
                    #pad2.SetGridy()
                    pad2.SetLogx()
                    pad2.SetTickx(1)
                    pad2.SetTicky(1)
                    pad2.SetLeftMargin(0.16)
                    pad2.SetRightMargin(0.02)
                    #pad2.SetTopMargin(0.05)
                    #pad2.SetBottomMargin(0.31)
                    pad2.SetTopMargin(0.065)
                    pad2.SetBottomMargin(0.13)
                    pad2.SetFrameFillStyle(0)
                    pad2.SetFrameBorderMode(0)
                    pad2.SetFrameFillStyle(0)
                    pad2.SetFrameBorderMode(0)
                    pad2.Draw()
                    # End Run2015B Setting 
       
                    # Plot the measurements in the top pad
                    pad1.cd()

                    histo = ROOT.TH2F('histo','',58,minCombinedPt,maxCombinedPt,100,1.-widthYAxis,1.+widthYAxis)
                    histo.SetLabelSize(0.05, 'XYZ')
                    histo.SetTitleSize(0.06, 'XYZ') 
                    histo.SetLabelFont(42, 'XYZ') 
                    histo.SetTitleFont(42, 'XYZ')
                    #histo.GetXaxis().SetTitle('Jet p_{T} [GeV]')
                    histo.GetXaxis().SetTitle('p_{T} [GeV]')
                    #histo.GetYaxis().SetTitle('Data/Simulation SF_{b}')
                    histo.GetYaxis().SetTitle('SF_{b}')
                    #histo.SetTitleOffset(1.1,'X') # Ideal for .png
                    histo.SetTitleOffset(0.95,'X')
                    histo.SetTitleOffset(0.8,'Y')
                    histo.SetTickLength(0.06,'X')
                    histo.SetNdivisions(509, 'XYZ')
                    histo.GetXaxis().SetMoreLogLabels()
                    histo.GetXaxis().SetNoExponent()
                    histo.Draw('')
                   
                    for method in graphScaleFactors: graphScaleFactors[method].Draw('P')
                    for method in graphScaleFactorsTotal: graphScaleFactorsTotal[method].Draw('P')                  
 
                    # Run2015B Style
                    tex = ROOT.TLatex(0.2,0.88,'CMS') 
                    tex.SetNDC() 
                    tex.SetTextAlign(13)
                    tex.SetTextFont(61)
                    tex.SetTextSize(0.07475)
                    tex.SetLineWidth(2) 
                    tex.Draw()                                                                                       
                    tex1 = ROOT.TLatex(0.2,0.79,'Preliminary') 
                    tex1.SetNDC()
                    tex1.SetTextAlign(13)
                    tex1.SetTextFont(52)
                    tex1.SetTextSize(0.05681)
                    tex1.SetLineWidth(2)   
                    tex1.Draw()   
                    text1 = ROOT.TLatex(0.98,0.95125, campaignLuminosity + centerOfMassEnergy) 
                    text1.SetNDC()                                              
                    text1.SetTextAlign(31)                          
                    text1.SetTextFont(42)    
                    text1.SetTextSize(0.04875)   
                    text1.SetLineWidth(2)    
                    text1.Draw() 
                    # End Run2015B Style

                    # Print the combination result in the top pad
                    graphCombinedScaleFactors.SetFillStyle(3005)
                    graphCombinedScaleFactors.SetFillColor(ROOT.kGray+3)
                    graphCombinedScaleFactors.Draw('e2')

                    # Add legend
                    plotHeader = algo + ' ' + ''.join([ x for x in wp if x.isupper() ])
                    leg1 = ROOT.TLegend(0.48,0.64,0.70,0.89)
                    leg1.SetBorderSize(0)
                    leg1.SetFillColor(ROOT.kWhite)
                    leg1.SetTextFont(62)
                    leg1.SetHeader(plotHeader) 
                    leg1.SetNColumns((len(graphScaleFactors))/3+1)

                    for method in graphScaleFactors:
	                leg1.AddEntry(graphScaleFactors[method], measurements[method]['legname'], 'PL')
	            leg1.AddEntry(graphCombinedScaleFactors,'weighted average','PF')

                    leg1.SetY1(0.89 - 0.25*leg1.GetNRows()/4)
                    if leg1.GetNColumns()>1: leg1.SetX2(0.92)
                    if leg1.GetNColumns()>2: 
                        leg1.SetX1(0.36) 
                        leg1.SetX2(0.95) 
                    leg1.Draw()
       
                    # Superimpose the single measurements on the combination and legend   
                    for method in graphScaleFactors: graphScaleFactors[method].Draw('P')
                    for method in graphScaleFactorsTotal: graphScaleFactorsTotal[method].Draw('P')

                    # Now plotting the combination in the bottom pad
                    pad2.cd()
                    histo.Draw('')
                    graphCombinedScaleFactors.Draw('e2')
       
                    # Add the functions values to the bottom pad
                    graphFittedScaleFactors = ROOT.TGraphErrors()

                    for iptbin, ptbin in enumerate(measuredScaleFactors):

                        minPt, maxPt = float(ptbin.split('-')[1].split('to')[0]), float(ptbin.split('to')[1])
                        midPt = (maxPt+minPt)/2.

                        fittedScaleFactor = fittingFunction.Integral(midPt-1.,midPt+1.)/2.
                        fittedScaleFactorError = combinedScaleFactorUncertaintyVector[iptbin][0]*fittedScaleFactor/combinedScaleFactorVector[iptbin][0]

                        graphFittedScaleFactors.SetPoint(iptbin, midPt, fittedScaleFactor)
                        graphFittedScaleFactors.SetPointError(iptbin, (maxPt-minPt)/2., fittedScaleFactorError)

                    graphFittedScaleFactors.SetMarkerColor(ROOT.kRed)
                    graphFittedScaleFactors.SetLineColor(ROOT.kRed)
                    graphFittedScaleFactors.SetLineStyle(1)
                    graphFittedScaleFactors.SetMarkerStyle(24)
                    graphFittedScaleFactors.SetMarkerSize(0.001)
                    graphFittedScaleFactors.SetLineWidth(2)
                    graphFittedScaleFactors.Draw('P') 
                    fittingFunction.SetMarkerColor(1)
                    fittingFunction.Draw('same') 
        
                    # Add legend
                    leg2 = ROOT.TLegend(0.48,0.64,0.70,0.89)
                    leg2.SetBorderSize(0)
                    leg2.SetFillColor(ROOT.kWhite)
                    leg2.SetTextFont(62)
                    leg2.SetTextSize(0.05)   
                    leg2.SetHeader(plotHeader)
                    leg2.AddEntry(graphCombinedScaleFactors,'weighted average','PF')
                    leg2.AddEntry(fittingFunction,'fit','L')
                    leg2.AddEntry(graphFittedScaleFactors,'fit #pm (stat #oplus syst)','LE')
                    leg2.Draw()
       
                    # Run2015B Style
                    tex2 = ROOT.TLatex(0.2,0.88,'CMS') 
                    tex2.SetNDC() 
                    tex2.SetTextAlign(13)
                    tex2.SetTextFont(61)
                    tex2.SetTextSize(0.07475)
                    tex2.SetLineWidth(2) 
                    tex2.Draw()                                                                                       
                    tex3 = ROOT.TLatex(0.2,0.79,'Preliminary') 
                    tex3.SetNDC()
                    tex3.SetTextAlign(13)
                    tex3.SetTextFont(52)
                    tex3.SetTextSize(0.05681)
                    tex3.SetLineWidth(2)   
                    tex3.Draw()      
       
                    text2 = ROOT.TLatex(0.98,0.95125, campaignLuminosity + centerOfMassEnergy) 
                    text2.SetNDC()                                              
                    text2.SetTextAlign(31)                          
                    text2.SetTextFont(42)    
                    text2.SetTextSize(0.04875)   
                    text2.SetLineWidth(2)    
                    text2.Draw()  
                    # End Run2015B Style 

                    # Save plots
                    plotDirectory = '/'.join([ opt.plotdir, opt.campaign, '' ])
                    os.system('mkdir -p '+plotDirectory+'; cp '+opt.plotdir+'/index.php '+plotDirectory)

                    plotTitleList = [ 'SFb', opt.campaign, plotHeader.replace(' ','') ]
                    for method in graphScaleFactors: plotTitleList.append(('masked' if method in maskedMethods else '')+measurements[method]['plotname'])

                    for fileExtension in [ '.png', '.pdf', '.root', '.C' ]:
                        c1.Print(plotDirectory+'_'.join(plotTitleList)+fileExtension)
                
        # Store the results of the scale factor combinations for this algorithm
        if opt.store:
            with open(opt.csvfiledir+'/'+'_'.join(csvFileNameList)+'.csv', 'w') as f:
                f.write(csvFile.makeCSV())

