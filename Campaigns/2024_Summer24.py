from collections import OrderedDict
import math
#import ROOT

# General data info
campaignLuminosity = '109 fb^{-1}'
centerOfMassEnergy = ' (13.6 TeV, 2024_Summer24)'

minPtCampaign  =   20.
maxPtCampaign  = 1000.
maxEtaCampaign = 2.5
widthYAxis     = 0.4

csvFileNameFlag = '130XSummer24SF'

# Taggers
algorithms = OrderedDict()
algorithms['UParTAK4']                   = [ 0., 1. ]

workingPoints = OrderedDict()
workingPoints['Loose']            =   'L'
workingPoints['Medium']           =   'M'
workingPoints['Tight']            =   'T'
workingPoints['eXtraTight']       =  'XT'
workingPoints['eXtraeXtraTight']  = 'XXT'

if opt.store:

    csvWorkingPoints = { 'Loose'            : ROOT.BTagEntry.OperatingPoint.OP_LOOSE,
                         'Medium'           : ROOT.BTagEntry.OperatingPoint.OP_MEDIUM,
                         'Tight'            : ROOT.BTagEntry.OperatingPoint.OP_TIGHT,
                         'eXtraTight'       : ROOT.BTagEntry.OperatingPoint.OP_EXTRATIGHT,
                         'eXtraeXtraTight'  : ROOT.BTagEntry.OperatingPoint.OP_EXTRAEXTRATIGHT,
                        } 

# Some general uncertainty options
sampleDependence         =   0. 
normalizedChi2Tollerance = 999.
cJetsInflationFactor     = { 'Loose' : 2.5, 'Medium' : 3.0, 'Tight' : 3.5, 'eXtraTight' : 3.5, 'eXtraeXtraTight' : 3.5 }

# Measurements
combinations = [ 'comb', 'mujets' ]

measurements =  OrderedDict()
measurements['ptrel']  = { 'plotname' : 'PtRel',  'data' : [ 'mujets', 'comb' ], 'ptbins' : 8, 'color' : 600, 'marker' : 22, 'size' : 1.3, 'shift' :  0., 'width' : 3 }
measurements['sys8']   = { 'plotname' : 'Syst8',  'data' : [ 'mujets', 'comb' ], 'ptbins' : 5, 'color' : 416, 'marker' : 29, 'size' : 1.3, 'shift' :  1., 'width' : 3 }
measurements['ltsv']   = { 'plotname' : 'LTSV',   'data' : [ 'mujets', 'comb' ], 'ptbins' : 9, 'color' : 632, 'marker' : 20, 'size' : 1.3, 'shift' : -1., 'width' : 2 }
measurements['kinfit'] = { 'plotname' : 'KinFit', 'data' : [ 'ttbar',  'comb' ], 'ptbins' : 7, 'color' : 1,   'marker' : 33, 'size' : 1.0, 'shift' : -2., 'width' : 2 }
measurements['tnp']    = { 'plotname' : 'TnP',    'data' : [ 'ttbar',  'comb' ], 'ptbins' : 6, 'color' : 30,  'marker' : 23, 'size' : 1.0, 'shift' :  2., 'width' : 2 }

for method in measurements:
    measurements[method]['version'] = 'last'
    measurements[method]['legname'] = measurements[method]['plotname'].replace('Syst8','System-8')

vetoedMethods = []
maskedMethods = []
maskedMeasurements = { }

# Function for pt-dependence fit
fittingFunctions, ptBinErrorScales = {}, {}
for comb in combinations:
    fittingFunctions[comb] = {}
    ptBinErrorScales[comb] = {}
    for algo in algorithms:
        fittingFunctions[comb][algo] = {}
        ptBinErrorScales[comb][algo] = {}
        for wp in workingPoints:
            ptBinErrorScales[comb][algo][wp] = {}
            if comb=='mujets':
              if wp=='eXtraeXtraTight' or wp=='eXtraTight':
                  fittingFunctions[comb][algo][wp] = ROOT.TF1('fittingFunction', '[0]+[1]*log(x-[3])+[2]*log(x-[3])*log(x-[3])', minPtCampaign, maxPtCampaign)
              elif wp=='Loose':
                  fittingFunctions[comb][algo][wp] = ROOT.TF1('fittingFunction', '[0]+[1]*log(x)+[2]*log(x)*log(x)', minPtCampaign, maxPtCampaign)
              else:
                  fittingFunctions[comb][algo][wp] = ROOT.TF1('fittingFunction', '[0]*(1.+[1]*x)/(1.+[2]*x)', minPtCampaign, maxPtCampaign)
                  #if wp=='Loose':
                  #    fittingFunctions[comb][algo][wp].SetParameters(0.94573766673456239, 0.010463597845190494, 0.010043776616054479)
            else:
                fittingFunctions[comb][algo][wp] = ROOT.TF1('fittingFunction', '[0]+[1]*log(x)+[2]*log(x)*log(x)', minPtCampaign, maxPtCampaign)
                #fittingFunctions[comb][algo][wp] = ROOT.TF1('fittingFunction', '[0]+[1]*log(x-[3])+[2]*log(x-[3])*log(x-[3])', minPtCampaign, maxPtCampaign)
                #fittingFunctions[comb][algo][wp] = ROOT.TF1('fittingFunction', '[0]*(1.+[1]*x)/(1.+[2]*x)', minPtCampaign, maxPtCampaign)

# Systematic breakdown categories
type1Systematics = [ 'statistic' ]
type2Systematics = [ 'pileup', 'bfragmentation', 'cfragmentation', 'jes', 'jer', 'topmass', 'mtop', 'hdamp', 'pdf', 'qcdscale', 'isrDef', 'fsrDef', 'isrdef', 'fsrdef', 'colorreconnection', 'tuneCP5', 'tunecp5', 'muf', 'mur', 'pdfas' ]  
type3Systematics = [ 'mupt', 'gluonsplitting', 'jetaway', 'mudr', 'cb', 'ttbarmodelling', 'l2c', 'ptrel', 'ipbias', 'ltothers', 'tt', 'kt', 'tct', 'ntrk', 'njet', 'jeta', 'dmux', 'ksl', 'ntrkgen', 'btempcorr', 'ltempcorr', 'flavFrac', 'bkg', 'scale1', 'ps', 'met', 'hpp', 'sel', 'trig', 'tWth', 'csv', 'cjets', 'mistag', 'nonttXsec', 'bdecays', 'mescales', 'psscale', 'tune' ]

# Systematic correlations
systematicPtCorrelated = [ 'pileup', 'bfragmentation', 'cfragmentation', 'cb', 'ltothers', 'jer', 'flavFrac', 'scale1', 'ps', 'topmass', 'mtop', 'hdamp', 'pdf',  'qcdscale', 'isrDef', 'fsrDef', 'isrdef', 'fsrdef', 'colorreconnection', 'tuneCP5', 'tunecp5', 'muf', 'mur', 'pdfas', 'bdecays', 'mescales', 'psscale', 'tune', 'jes' ]
systematicPtUncorrelated = [ 'statistic', 'gluonsplitting', 'jetaway', 'ttbarmodelling', 'l2c', 'ptrel', 'ipbias', 'tt', 'kt', 'tct', 'ntrk', 'njet', 'jeta', 'dmux', 'ksl', 'ntrkgen', 'btempcorr', 'ltempcorr', 'bkg', 'met', 'hpp', 'sel', 'trig', 'tWth', 'csv', 'cjets', 'mistag', 'nonttXsec', 'mupt', 'mudr' ]
ptCorrelationCoefficients = { 'ltothers' : 0.5 }

systematicYearCorrelated = [ 'pileup', 'bfragmentation', 'cfragmentation', 'topmass', 'mtop', 'hdamp', 'pdf', 'qcdscale', 'isrDef', 'fsrDef', 'isrdef', 'fsrdef', 'colorreconnection', 'tuneCP5', 'tunecp5', 'muf', 'mur', 'pdfas' ]
systematicYearCorrelated.extend(type3Systematics)
systematicYearUncorrelated = [ 'jes', 'jer' ]
systematicYearUncorrelated.extend(type1Systematics)

# Statistical correlations
statisticalCorrelationCoefficients = {}
statisticalCorrelationCoefficients['statistic'] = {}
statisticalCorrelationCoefficients['statistic']['ltsv-ptrel'] = {}
statisticalCorrelationCoefficients['statistic']['ltsv-ptrel']['Pt-20to30']    = 0.09
statisticalCorrelationCoefficients['statistic']['ltsv-ptrel']['Pt-30to50']    = 0.096491
statisticalCorrelationCoefficients['statistic']['ltsv-ptrel']['Pt-50to70']    = 0.128465
statisticalCorrelationCoefficients['statistic']['ltsv-ptrel']['Pt-70to100']   = 0.162212
statisticalCorrelationCoefficients['statistic']['ltsv-ptrel']['Pt-100to140']  = 0.182612
statisticalCorrelationCoefficients['statistic']['ltsv-ptrel']['Pt-140to200']  = 0.203323
statisticalCorrelationCoefficients['statistic']['ltsv-ptrel']['Pt-200to300']  = 0.235041
statisticalCorrelationCoefficients['statistic']['ltsv-ptrel']['Pt-300to600']  = 0.239796
statisticalCorrelationCoefficients['statistic']['ltsv-ptrel']['Pt-600to1000'] = 0.24
statisticalCorrelationCoefficients['statistic']['ltsv-sys8'] = {}
statisticalCorrelationCoefficients['statistic']['ltsv-sys8']['Pt-20to30']    = 0.20
statisticalCorrelationCoefficients['statistic']['ltsv-sys8']['Pt-30to50']    = 0.206224
statisticalCorrelationCoefficients['statistic']['ltsv-sys8']['Pt-50to70']    = 0.230399
statisticalCorrelationCoefficients['statistic']['ltsv-sys8']['Pt-70to100']   = 0.310995
statisticalCorrelationCoefficients['statistic']['ltsv-sys8']['Pt-100to140']  = 0.331386
statisticalCorrelationCoefficients['statistic']['ltsv-sys8']['Pt-140to200']  = 0.356070
statisticalCorrelationCoefficients['statistic']['ltsv-sys8']['Pt-200to300']  = 0.403861
statisticalCorrelationCoefficients['statistic']['ltsv-sys8']['Pt-300to600']  = 0.432426
statisticalCorrelationCoefficients['statistic']['ltsv-sys8']['Pt-600to1000'] = 0.43
statisticalCorrelationCoefficients['statistic']['sys8-ptrel'] = {}
statisticalCorrelationCoefficients['statistic']['sys8-ptrel']['Pt-20to30']    = 0.46
statisticalCorrelationCoefficients['statistic']['sys8-ptrel']['Pt-30to50']    = 0.467892
statisticalCorrelationCoefficients['statistic']['sys8-ptrel']['Pt-50to70']    = 0.488014
statisticalCorrelationCoefficients['statistic']['sys8-ptrel']['Pt-70to100']   = 0.413078
statisticalCorrelationCoefficients['statistic']['sys8-ptrel']['Pt-100to140']  = 0.434579
statisticalCorrelationCoefficients['statistic']['sys8-ptrel']['Pt-140to200']  = 0.44793
statisticalCorrelationCoefficients['statistic']['sys8-ptrel']['Pt-200to300']  = 0.401653
statisticalCorrelationCoefficients['statistic']['sys8-ptrel']['Pt-300to600']  = 0.413572
statisticalCorrelationCoefficients['statistic']['sys8-ptrel']['Pt-600to1000'] = 0.41

for methodPair in list(statisticalCorrelationCoefficients['statistic'].keys()):
    for ptbin in list(statisticalCorrelationCoefficients['statistic'][methodPair].keys()):
        statisticalCorrelationCoefficients['statistic'][methodPair][ptbin] = math.sqrt(statisticalCorrelationCoefficients['statistic'][methodPair][ptbin])
    invertedMethodPair = '-'.join([ methodPair.split('-')[1],  methodPair.split('-')[0] ])
    statisticalCorrelationCoefficients['statistic'][invertedMethodPair] = statisticalCorrelationCoefficients['statistic'][methodPair]

# TTbar global measurement
chi2Strategy = 'Overall'

# This is for avereging the results on ttbar spectrum
ttbarPt =  [ 0., 85., 50., 36., 26., 14., 2.1, 0. ]
ttbarSpectrumScale = 1.


