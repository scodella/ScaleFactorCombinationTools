from collections import OrderedDict
import math
import ROOT

# General data info
campaignLuminosity = '19.5 fb^{-1}'
centerOfMassEnergy = ' (13 TeV, 2016preVFP)'

minPtCampaign  =   20.
maxPtCampaign  = 1000.
maxEtaCampaign = 2.4
widthYAxis     = 0.4

# Taggers
algorithms = OrderedDict()
algorithms['DeepCSV'] = [ 0., 1. ]
algorithms['DeepJet'] = [ 0., 1. ]

workingPoints = OrderedDict()
workingPoints['Loose']  = '0'
workingPoints['Medium'] = '1'
workingPoints['Tight']  = '2'

csvWorkingPoints = { 'Loose'  : ROOT.BTagEntry.OperatingPoint.OP_LOOSE,
                     'Medium' : ROOT.BTagEntry.OperatingPoint.OP_MEDIUM,
                     'Tight'  : ROOT.BTagEntry.OperatingPoint.OP_TIGHT,
                    } 

# Some general uncertainty options
sampleDependence         =   0. 
normalizedChi2Tollerance = 999.
cJetsInflationFactor     = { 'Loose' : 2.5, 'Medium' : 3.0, 'Tight' : 3.5 }

# Measurements
combinations = [ 'comb', 'mujets' ]

measurements =  OrderedDict()
measurements['PtRel']   = { 'plotname' : 'PT',  'data' : [ 'mujets', 'comb' ], 'ptbins' : 8, 'color' : 600, 'marker' : 22, 'size' : 1.3, 'shift' :  0., 'width' : 3 }
measurements['System8'] = { 'plotname' : 'S8',  'data' : [ 'mujets', 'comb' ], 'ptbins' : 5, 'color' : 416, 'marker' : 29, 'size' : 1.3, 'shift' :  1., 'width' : 3 }
measurements['LT']      = { 'plotname' : 'LT',  'data' : [ 'mujets', 'comb' ], 'ptbins' : 9, 'color' : 632, 'marker' : 20, 'size' : 1.3, 'shift' : -1., 'width' : 2 }
#measurements['TwoTag']  = { 'plotname' : 'TT',  'data' : [ 'ttbar',  'comb' ], 'ptbins' : 6, 'color' : 800, 'marker' : 21, 'size' : 1.0, 'shift' :  2., 'width' : 2 }
measurements['Kin']     = { 'plotname' : 'Kin', 'data' : [ 'ttbar',  'comb' ], 'ptbins' : 7, 'color' : 1,   'marker' : 33, 'size' : 1.0, 'shift' : -2., 'width' : 2 }
measurements['TnP']     = { 'plotname' : 'TnP', 'data' : [ 'ttbar',  'comb' ], 'ptbins' : 6, 'color' : 30,  'marker' : 23, 'size' : 1.0, 'shift' :  3., 'width' : 2 }

for method in measurements:
    measurements[method]['csvname'] = method.lower()
    measurements[method]['legname'] = method.replace('System8','System-8')

vetoedMethods = []
maskedMethods = []
maskedMeasurements = { 'Kin' : [ 'Pt-300to600' ] }

# Function for pt-dependence fit
fittingFunctions = {}
for comb in combinations:
    fittingFunctions[comb] = {}
    for algo in algorithms:
        fittingFunctions[comb][algo] = {}
        for wp in workingPoints:
            if algo=='DeepCSV' and (wp=='Tight' or (wp=='Medium' and comb=='comb')):
                fittingFunctions[comb][algo][wp] = ROOT.TF1('fittingFunction', '[0]+[1]*log(x+19)*log(x+18)*(3-[2]*log(x+18))', minPtCampaign, maxPtCampaign)
            else:
                fittingFunctions[comb][algo][wp] = ROOT.TF1('fittingFunction', '[0]*(1.+[1]*x)/(1.+[2]*x)', minPtCampaign, maxPtCampaign)

# Systematic breakdown categories
type1Systematics = [ 'statistic' ]
type2Systematics = [ 'pileup', 'bfragmentation', 'cfragmentation', 'jes', 'jer', 'topmass', 'mtop', 'hdamp', 'qcdscale', 'isrDef', 'fsrDef', 'isrdef', 'fsrdef', 'tuneCP5', 'tunecp5' ]  
type3Systematics = [ 'mupt', 'gluonsplitting', 'jetaway', 'mudr', 'cb', 'ttbarmodelling', 'l2c', 'ptrel', 'ipbias', 'ltothers', 'tt', 'kt', 'tct', 'ntrk', 'njet', 'jeta', 'dmux', 'ksl', 'ntrkgen', 'btempcorr', 'ltempcorr', 'flavFrac', 'bkg', 'scale1', 'ps', 'met', 'hpp', 'sel', 'trig', 'tWth', 'csv', 'cjets', 'mistag', 'nonttXsec', 'bdecay', 'mescales', 'psscale', 'tune' ]

# Systematic correlations
systematicPtCorrelated = [ 'pileup', 'mupt', 'bfragmentation', 'cfragmentation', 'mudr', 'cb', 'jes', 'ltothers', 'jer', 'flavFrac', 'scale1', 'ps', 'topmass', 'mtop', 'hdamp', 'qcdscale', 'isrDef', 'fsrDef', 'isrdef', 'fsrdef', 'tuneCP5', 'tunecp5', 'bdecay', 'mescales', 'psscale', 'tune' ]
systematicPtUncorrelated = [ 'statistic', 'gluonsplitting', 'jetaway', 'ttbarmodelling', 'l2c', 'ptrel', 'ipbias', 'tt', 'kt', 'tct', 'ntrk', 'njet', 'jeta', 'dmux', 'ksl', 'ntrkgen', 'btempcorr', 'ltempcorr', 'bkg', 'met', 'hpp', 'sel', 'trig', 'tWth', 'csv', 'cjets', 'mistag', 'nonttXsec' ]
ptCorrelationCoefficients = { 'ltothers' : 0.5 }

systematicYearCorrelated = [ 'pileup', 'bfragmentation', 'cfragmentation', 'topmass', 'mtop', 'hdamp', 'qcdscale', 'isrDef', 'fsrDef', 'isrdef', 'fsrdef', 'tuneCP5', 'tunecp5' ]
systematicYearCorrelated.extend(type3Systematics)
systematicYearUncorrelated = [ 'jes', 'jer' ]
systematicYearUncorrelated.extend(type1Systematics)

# Statistical correlations
statisticalCorrelationCoefficients = {}
statisticalCorrelationCoefficients['statistic'] = {}
statisticalCorrelationCoefficients['statistic']['LT-PtRel'] = {}
statisticalCorrelationCoefficients['statistic']['LT-PtRel']['Pt-20to30']    = 0.09
statisticalCorrelationCoefficients['statistic']['LT-PtRel']['Pt-30to50']    = 0.096491
statisticalCorrelationCoefficients['statistic']['LT-PtRel']['Pt-50to70']    = 0.128465
statisticalCorrelationCoefficients['statistic']['LT-PtRel']['Pt-70to100']   = 0.162212
statisticalCorrelationCoefficients['statistic']['LT-PtRel']['Pt-100to140']  = 0.182612
statisticalCorrelationCoefficients['statistic']['LT-PtRel']['Pt-140to200']  = 0.203323
statisticalCorrelationCoefficients['statistic']['LT-PtRel']['Pt-200to300']  = 0.235041
statisticalCorrelationCoefficients['statistic']['LT-PtRel']['Pt-300to600']  = 0.239796
statisticalCorrelationCoefficients['statistic']['LT-PtRel']['Pt-600to1000'] = 0.24
statisticalCorrelationCoefficients['statistic']['LT-System8'] = {}
statisticalCorrelationCoefficients['statistic']['LT-System8']['Pt-20to30']    = 0.20
statisticalCorrelationCoefficients['statistic']['LT-System8']['Pt-30to50']    = 0.206224
statisticalCorrelationCoefficients['statistic']['LT-System8']['Pt-50to70']    = 0.230399
statisticalCorrelationCoefficients['statistic']['LT-System8']['Pt-70to100']   = 0.310995
statisticalCorrelationCoefficients['statistic']['LT-System8']['Pt-100to140']  = 0.331386
statisticalCorrelationCoefficients['statistic']['LT-System8']['Pt-140to200']  = 0.356070
statisticalCorrelationCoefficients['statistic']['LT-System8']['Pt-200to300']  = 0.403861
statisticalCorrelationCoefficients['statistic']['LT-System8']['Pt-300to600']  = 0.432426
statisticalCorrelationCoefficients['statistic']['LT-System8']['Pt-600to1000'] = 0.43
statisticalCorrelationCoefficients['statistic']['System8-PtRel'] = {}
statisticalCorrelationCoefficients['statistic']['System8-PtRel']['Pt-20to30']    = 0.46
statisticalCorrelationCoefficients['statistic']['System8-PtRel']['Pt-30to50']    = 0.467892
statisticalCorrelationCoefficients['statistic']['System8-PtRel']['Pt-50to70']    = 0.488014
statisticalCorrelationCoefficients['statistic']['System8-PtRel']['Pt-70to100']   = 0.413078
statisticalCorrelationCoefficients['statistic']['System8-PtRel']['Pt-100to140']  = 0.434579
statisticalCorrelationCoefficients['statistic']['System8-PtRel']['Pt-140to200']  = 0.44793
statisticalCorrelationCoefficients['statistic']['System8-PtRel']['Pt-200to300']  = 0.401653
statisticalCorrelationCoefficients['statistic']['System8-PtRel']['Pt-300to600']  = 0.413572
statisticalCorrelationCoefficients['statistic']['System8-PtRel']['Pt-600to1000'] = 0.41

for methodPair in statisticalCorrelationCoefficients['statistic'].keys():
    for ptbin in statisticalCorrelationCoefficients['statistic'][methodPair].keys():
        statisticalCorrelationCoefficients['statistic'][methodPair][ptbin] = math.sqrt(statisticalCorrelationCoefficients['statistic'][methodPair][ptbin])
    invertedMethodPair = '-'.join([ methodPair.split('-')[1],  methodPair.split('-')[0] ])
    statisticalCorrelationCoefficients['statistic'][invertedMethodPair] = statisticalCorrelationCoefficients['statistic'][methodPair]

# TTbar global measurement
chi2Strategy = 'Overall'

# This is for avereging the results on ttbar spectrum
ttbarPt =  [ 0., 85., 50., 36., 26., 14., 2.1, 0. ]
ttbarSpectrumScale = 1.


