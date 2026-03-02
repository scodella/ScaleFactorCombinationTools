from collections import OrderedDict
import math
#import ROOT

# General data info
campaignLuminosity = '17.1 fb^{-1}'
centerOfMassEnergy = ' (13.6 TeV, 2023_Summer23)'

minPtCampaign  =   20.
maxPtCampaign  = 1000.
maxEtaCampaign = 2.5
widthYAxis     = 0.4

csvFileNameFlag = '130XSummer23SF'

# Taggers
algorithms = OrderedDict()
algorithms['deepJet']                   = [ 0., 1. ]
algorithms['particleNet']               = [ 0., 1. ]
algorithms['robustParticleTransformer'] = [ 0., 1. ]

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
            if comb=='mujets' and (algo=='deepJet' or algo=='robustParticleTransformer'):
                if ((algo=='deepJet' and wp=='eXtraTight') or algo=='robustParticleTransformer') and 'ltsv' not in opt.vetomethod and 'ltsv' not in opt.maskmethod:
                    fittingFunctions[comb][algo][wp] = ROOT.TF1('fittingFunction', '[0]+[1]*log(x-[3])+[2]*log(x-[3])*log(x-[3])', minPtCampaign, maxPtCampaign)
                else:
                    fittingFunctions[comb][algo][wp] = ROOT.TF1('fittingFunction', '[0]*(1.+[1]*x)/(1.+[2]*x)', minPtCampaign, maxPtCampaign)
            elif comb=='comb' and ('kinfit' in opt.vetomethod or 'kinfit' in opt.maskmethod):
                fittingFunctions[comb][algo][wp] = ROOT.TF1('fittingFunction', '[0]+[1]*log(x)+[2]*log(x)*log(x)', minPtCampaign, maxPtCampaign)
                if algo=='deepJet':
                    if wp=='Loose': fittingFunctions[comb][algo][wp].SetParameters(0.354212,0.237673,-0.0218838,0.)
                    if wp=='Medium': fittingFunctions[comb][algo][wp].SetParameters(1.47533,-0.24666,0.0296768,0.)
                    if wp=='Tight': fittingFunctions[comb][algo][wp].SetParameters(0.519429,0.156053,-0.0129339,0.)
                    if wp=='eXtraTight': fittingFunctions[comb][algo][wp].SetParameters(0.339425,0.217108,-0.0177989,0.)
                    if wp=='eXtraeXtraTight': 
                        #fittingFunctions[comb][algo][wp] = ROOT.TF1('fittingFunction', '[0]+[1]*log(x-[3])+[2]*log(x-[3])*log(x-[3])', minPtCampaign, maxPtCampaign)
                        fittingFunctions[comb][algo][wp].SetParameters(-0.471009,0.520799,-0.0465554,0.)
                elif algo=='robustParticleTransformer':
                    if wp=='Loose': fittingFunctions[comb][algo][wp].SetParameters(-0.12055,0.423318,-0.0399768,0.)
                    if wp=='Medium': fittingFunctions[comb][algo][wp].SetParameters(0.990195,0.0139309,-0.00290727,0.)
                    if wp=='Tight': fittingFunctions[comb][algo][wp].SetParameters(0.951052,0.0139379,-0.00183577,0.)
                    if wp=='eXtraTight': fittingFunctions[comb][algo][wp].SetParameters(-1.07075,0.785628,-0.0743669,0.)  
                    if wp=='eXtraeXtraTight': fittingFunctions[comb][algo][wp].SetParameters(-0.741095,0.626532,-0.0567302,0.)
            elif comb=='comb' and algo=='particleNet':
                fittingFunctions[comb][algo][wp] = ROOT.TF1('fittingFunction', '[0]+[1]*log(x-[3])+[2]*log(x-[3])*log(x-[3])', minPtCampaign, maxPtCampaign)
                if wp=='Loose':
                    fittingFunctions[comb][algo][wp] = ROOT.TF1('fittingFunction', '[0]+[1]*log(x)+[2]*log(x)*log(x)', minPtCampaign, maxPtCampaign)
                    fittingFunctions[comb][algo][wp].SetParameters(2.03869,-0.436488,0.0447336,0.)
                    #fittingFunctions[comb][algo][wp].SetParameters(4.04675,-1.21536,0.120034,0.)
                elif wp=='Medium':
                    fittingFunctions[comb][algo][wp] = ROOT.TF1('fittingFunction', '[0]+[1]*log(x)+[2]*log(x)*log(x)', minPtCampaign, maxPtCampaign)
                    fittingFunctions[comb][algo][wp].SetParameters(2.03869,-0.436488,0.0447336,0.)
                elif wp=='Tight':
                    fittingFunctions[comb][algo][wp].SetParameters(0.805525,0.0295216,-2.32352e-06,0.)
                elif wp=='eXtraTight':
                     fittingFunctions[comb][algo][wp].SetParameters(0.930226,-0.0325424,0.00663287,0.)
                elif wp=='eXtraeXtraTight':
                     fittingFunctions[comb][algo][wp].SetParameters(-70.0757,21.7726,-1.6688,-535.352)
            elif comb=='comb' and algo=='deepJet':
                fittingFunctions[comb][algo][wp] = ROOT.TF1('fittingFunction', '[0]+[1]*log(x-[3])+[2]*log(x-[3])*log(x-[3])', minPtCampaign, maxPtCampaign)
                if wp=='Loose':
                    fittingFunctions[comb][algo][wp] = ROOT.TF1('fittingFunction', '[0]+[1]*log(x)+[2]*log(x)*log(x)', minPtCampaign, maxPtCampaign)
                    fittingFunctions[comb][algo][wp].SetParameters(2.15295,-0.389697,0.0306347,0.)
                elif wp=='Medium':
                    fittingFunctions[comb][algo][wp].SetParameters(2.02773,-0.427995,0.0428375,0.)
                elif wp=='Tight':
                    fittingFunctions[comb][algo][wp] = ROOT.TF1('fittingFunction', '[0]+[1]*log(x)+[2]*log(x)*log(x)', minPtCampaign, maxPtCampaign)
                    fittingFunctions[comb][algo][wp].SetParameters(1.44515,-0.203004,0.020718,0.)
                elif wp=='eXtraTight':
                     fittingFunctions[comb][algo][wp].SetParameters(-0.0724806,0.465473,-0.0523755,0.)
                elif wp=='eXtraeXtraTight':
                     fittingFunctions[comb][algo][wp].SetParameters(-73.7967,21.728,-1.57902,-764.933)
            elif comb=='comb' and algo=='robustParticleTransformer':
                fittingFunctions[comb][algo][wp] = ROOT.TF1('fittingFunction', '[0]+[1]*log(x-[3])+[2]*log(x-[3])*log(x-[3])', minPtCampaign, maxPtCampaign)
                if wp=='Loose':
                    fittingFunctions[comb][algo][wp] = ROOT.TF1('fittingFunction', '[0]+[1]*log(x)+[2]*log(x)*log(x)', minPtCampaign, maxPtCampaign)
                    fittingFunctions[comb][algo][wp].SetParameters(2.15295,-0.389697,0.0306345,0.)
                elif wp=='Medium':
                    fittingFunctions[comb][algo][wp] = ROOT.TF1('fittingFunction', '[0]+[1]*log(x)+[2]*log(x)*log(x)', minPtCampaign, maxPtCampaign)
                    fittingFunctions[comb][algo][wp].SetParameters(1.62248,-0.252638,0.0243624,0.)
                elif wp=='Tight':
                    fittingFunctions[comb][algo][wp].SetParameters(1.44515,-0.203004,0.020718,0.)
                elif wp=='eXtraTight':
                    fittingFunctions[comb][algo][wp] = ROOT.TF1('fittingFunction', '[0]+[1]*log(x)+[2]*log(x)*log(x)', minPtCampaign, maxPtCampaign)
                    fittingFunctions[comb][algo][wp].SetParameters(-0.0724806,0.465473,-0.0523755,0.)
                elif wp=='eXtraeXtraTight':
                    fittingFunctions[comb][algo][wp] = ROOT.TF1('fittingFunction', '[0]+[1]*log(x)+[2]*log(x)*log(x)', minPtCampaign, maxPtCampaign)
                    fittingFunctions[comb][algo][wp].SetParameters(-73.7967,21.728,-1.57902,-764.933)

if 'ltsv' not in opt.vetomethod and 'ltsv' not in opt.maskmethod:
    ptBinErrorScales['mujets']['robustParticleTransformer']['Medium']['Pt-20to30'] = 1.5

if 'kinfit' in opt.vetomethod or 'kinfit' in opt.maskmethod:
    pass

else:
    ptBinErrorScales['comb']['particleNet']['Loose']['Pt-100to140'] = 100.
    ptBinErrorScales['comb']['particleNet']['Medium']['Pt-140to200'] = 100.
    ptBinErrorScales['comb']['particleNet']['Tight']['Pt-50to70'] = 100.
    ptBinErrorScales['comb']['particleNet']['Tight']['Pt-300to600'] = 0.1
    ptBinErrorScales['comb']['particleNet']['eXtraTight']['Pt-50to70'] = 100.
    ptBinErrorScales['comb']['particleNet']['eXtraTight']['Pt-300to600'] = 0.1
    ptBinErrorScales['comb']['particleNet']['eXtraeXtraTight']['Pt-50to70'] = 100.
    for pTbin in 'Pt-70to100', 'Pt-100to140', 'Pt-140to200', 'Pt-200to300':
        ptBinErrorScales['comb']['particleNet']['Tight'][pTbin] = 5.
        ptBinErrorScales['comb']['particleNet']['eXtraTight'][pTbin] = 5.
        ptBinErrorScales['comb']['particleNet']['eXtraeXtraTight'][pTbin] = 5.
    ptBinErrorScales['comb']['particleNet']['eXtraeXtraTight']['Pt-300to600'] = 0.1

    ptBinErrorScales['comb']['deepJet']['Loose']['Pt-20to30'] = 0.5
    ptBinErrorScales['comb']['deepJet']['Loose']['Pt-30to50'] = 100.
    ptBinErrorScales['comb']['deepJet']['Loose']['Pt-100to140'] = 100.
    ptBinErrorScales['comb']['deepJet']['Loose']['Pt-300to600'] = 100.
    ptBinErrorScales['comb']['deepJet']['Medium']['Pt-300to600'] = 3.
    for pTbin in 'Pt-30to50', 'Pt-50to70' 'Pt-70to100', 'Pt-100to140', 'Pt-140to200', 'Pt-200to300':
        ptBinErrorScales['comb']['deepJet']['Tight'][pTbin] = 5.
        ptBinErrorScales['comb']['deepJet']['eXtraTight'][pTbin] = 5.
    ptBinErrorScales['comb']['deepJet']['eXtraeXtraTight']['Pt-50to70'] = 100.
    ptBinErrorScales['comb']['deepJet']['eXtraeXtraTight']['Pt-300to600'] = 0.1
    ptBinErrorScales['comb']['deepJet']['eXtraeXtraTight']['Pt-600to1000'] = 0.1

    ptBinErrorScales['comb']['robustParticleTransformer']['Loose']['Pt-20to30'] = 0.5
    ptBinErrorScales['comb']['robustParticleTransformer']['Loose']['Pt-30to50'] = 100.
    ptBinErrorScales['comb']['robustParticleTransformer']['Loose']['Pt-70to100'] = 100.
    ptBinErrorScales['comb']['robustParticleTransformer']['Loose']['Pt-100to140'] = 100.
    ptBinErrorScales['comb']['robustParticleTransformer']['Loose']['Pt-300to600'] = 100.
    ptBinErrorScales['comb']['robustParticleTransformer']['Medium']['Pt-70to100'] = 100.
    ptBinErrorScales['comb']['robustParticleTransformer']['Medium']['Pt-140to200'] = 100.
    ptBinErrorScales['comb']['robustParticleTransformer']['Tight']['Pt-30to50'] = 100.
    ptBinErrorScales['comb']['robustParticleTransformer']['Tight']['Pt-140to200'] = 100.
    ptBinErrorScales['comb']['robustParticleTransformer']['eXtraTight']['Pt-50to70'] = 100.
    ptBinErrorScales['comb']['robustParticleTransformer']['eXtraTight']['Pt-140to200'] = 100.
    ptBinErrorScales['comb']['robustParticleTransformer']['eXtraeXtraTight']['Pt-50to70'] = 100.
    ptBinErrorScales['comb']['robustParticleTransformer']['eXtraeXtraTight']['Pt-140to200'] = 100.

# Systematic breakdown categories
type1Systematics = [ 'statistic' ]
type2Systematics = [ 'pileup', 'bfragmentation', 'cfragmentation', 'jes', 'jer', 'topmass', 'mtop', 'hdamp', 'pdf', 'qcdscale', 'isrDef', 'fsrDef', 'isrdef', 'fsrdef', 'colorreconnection', 'tuneCP5', 'tunecp5' ]  
type3Systematics = [ 'mupt', 'gluonsplitting', 'jetaway', 'mudr', 'cb', 'ttbarmodelling', 'l2c', 'ptrel', 'ipbias', 'ltothers', 'tt', 'kt', 'tct', 'ntrk', 'njet', 'jeta', 'dmux', 'ksl', 'ntrkgen', 'btempcorr', 'ltempcorr', 'flavFrac', 'bkg', 'scale1', 'ps', 'met', 'hpp', 'sel', 'trig', 'tWth', 'csv', 'cjets', 'mistag', 'nonttXsec', 'bdecay', 'mescales', 'psscale', 'tune' ]

# Systematic correlations
systematicPtCorrelated = [ 'pileup', 'bfragmentation', 'cfragmentation', 'cb', 'ltothers', 'jer', 'flavFrac', 'scale1', 'ps', 'topmass', 'mtop', 'hdamp', 'pdf',  'qcdscale', 'isrDef', 'fsrDef', 'isrdef', 'fsrdef', 'colorreconnection', 'tuneCP5', 'tunecp5', 'bdecay', 'mescales', 'psscale', 'tune', 'jes' ]
systematicPtUncorrelated = [ 'statistic', 'gluonsplitting', 'jetaway', 'ttbarmodelling', 'l2c', 'ptrel', 'ipbias', 'tt', 'kt', 'tct', 'ntrk', 'njet', 'jeta', 'dmux', 'ksl', 'ntrkgen', 'btempcorr', 'ltempcorr', 'bkg', 'met', 'hpp', 'sel', 'trig', 'tWth', 'csv', 'cjets', 'mistag', 'nonttXsec', 'mupt', 'mudr' ]
ptCorrelationCoefficients = { 'ltothers' : 0.5 }

systematicYearCorrelated = [ 'pileup', 'bfragmentation', 'cfragmentation', 'topmass', 'mtop', 'hdamp', 'pdf', 'qcdscale', 'isrDef', 'fsrDef', 'isrdef', 'fsrdef', 'colorreconnection', 'tuneCP5', 'tunecp5' ]
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


