#!/usr/bin/env python
import os
import sys
#!/usr/bin/env python
import ROOT
import math
import optparse
from array import *
import CondTools.BTau.dataLoader as dataLoader
from collections import defaultdict

if __name__ == '__main__':

    # Input parameters
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    
    parser.add_option('--inputfile'   , dest='inputfile'   , help='Name of the input csv file '           , default='')
    parser.add_option('--outputfile'  , dest='outputfile'  , help='Name of the output csv file'           , default='default')
    parser.add_option('--csvpath  '   , dest='csvpath'     , help='Path where csv files are stored '      , default='./CSVFiles/')
    parser.add_option('--yearcorroff' , dest='yearcorroff' , help='Turn off year correlations'            , default=False, action='store_true')
    parser.add_option('--splittype2'  , dest='splittype2'  , help='Split type2 uncertainties'             , default=False, action='store_true')
    parser.add_option('--custom'      , dest='custom'      , help='Custom list of wanted uncertainties'   , default=None)
    (opt, args) = parser.parse_args()

    uncorrelatedList = [ 'statistic' ] # To be completed

    type2List = [ "_pileup", "_jes", "_jer" ] # To be checked
    
    splitList = opt.custom.split(',') if (opt.custom!=None) else [ ] 
        
    outputFlag = ''
    if not opt.yearcorroff:
        outputFlag += '_years'
    if opt.splittype2:
        outputFlag += '_type2'
    if outputFlag=='':
        outputFlag = '_basic'
    for custom in splitList:
        outputFlag += '_' + custom
        
    if opt.splittype2:
        splitList.extend(type2List)
    
    outputfilename = opt.csvpath + '/' + opt.outputfile.replace('default', opt.inputfile.replace('.csv', '')) + outputFlag + '.csv'

    loaders = dataLoader.get_data(opt.csvpath + '/' + opt.inputfile)

    calib = ROOT.BTagCalibration(outputfilename.split('/')[-1].split('_')[0])

    mergedSystematics = [ ] 
    mergedUncertainty = { }

    if outputFlag!='_basic':

        mergedSystematics = [ 'correlated', 'uncorrelated' ] if ('years' in outputFlag) else [ 'unsplit' ] 

        # Get the structure of the CSV file
        for data in loaders:
            for e in data.entries:

                if e.params.sysType=='up':

                    for mergedSystematic in mergedSystematics:

                        paramList = [ e.params.operatingPoint, e.params.measurementType, mergedSystematic, e.params.jetFlavor, e.params.etaMin, e.params.etaMax, e.params.ptMin, e.params.ptMax, e.params.discrMin, e.params.discrMax ]

                        newMergedUncertainty = 0.

                        for param in reversed(paramList):
                            tmpMergedUncertainty = { } 
                            tmpMergedUncertainty[param] = newMergedUncertainty
                            newMergedUncertainty = tmpMergedUncertainty

                        newMergedUncertaintyList = [ newMergedUncertainty ] 
                        auxMergedUncertaintyList = [ mergedUncertainty ] 
            
                        for par1 in range(len(paramList)):
                            if paramList[par1] in auxMergedUncertaintyList[par1]:
                                auxMergedUncertaintyList.append(auxMergedUncertaintyList[par1][paramList[par1]])
                                newMergedUncertaintyList.append(newMergedUncertaintyList[par1][paramList[par1]])
                            else:
                                auxMergedUncertaintyList[par1][paramList[par1]] = newMergedUncertaintyList[par1][paramList[par1]]
                                break

        # Merge uncertainties
        for data in loaders:
            for e in data.entries:

                # Systematics are symmetric 
                if 'up_' in e.params.sysType and e.params.sysType.split('_')[1] not in splitList:

                    mergedSystematic = 'unsplit'
                    if 'years' in outputFlag:
                        if e.params.sysType.split('_')[1] in uncorrelatedList:
                            mergedSystematic = 'uncorrelated'
                        else:
                            mergedSystematic = 'correlated'

                    systematicValue = e.formula
                    if '+' in systematicValue: # This is ok for b-tag SF combination
                        systematicValue = float(systematicValue.split('+')[-1])

                    mergedUncertainty[e.params.operatingPoint][e.params.measurementType][mergedSystematic][e.params.jetFlavor][e.params.etaMin][e.params.etaMax][e.params.ptMin][e.params.ptMax][e.params.discrMin][e.params.discrMax] += systematicValue*systematicValue

    # Now we can write the output: first the cetral values, total uncertainties, and splitted systematics ...
    for data in loaders:
        for e in data.entries:
     
            if e.params.sysType=='central' or e.params.sysType=='up' or e.params.sysType=='down' or e.params.sysType.split('_')[1] in splitList:

                params = ROOT.BTagEntry.Parameters( e.params.operatingPoint, 
                                                    e.params.measurementType, 
                                                    e.params.sysType,
                                                    e.params.jetFlavor,
                                                    e.params.etaMin,
                                                    e.params.etaMax, 
                                                    e.params.ptMin, 
                                                    e.params.ptMax,
                                                    e.params.discrMin,
                                                    e.params.discrMax )
                
                SFFun = ROOT.TF1 ('SFFun', e.formula, e.params.ptMin, e.params.ptMax)
                entry = ROOT.BTagEntry(SFFun, params)
                calib.addEntry(entry)
                
    # ... then the merged systematics
    for mergedSystematic in mergedSystematics:
        for data in loaders:
            for e in data.entries:
                if e.params.sysType=='central': 

                    auxMergedUncertainty = mergedUncertainty[e.params.operatingPoint][e.params.measurementType][mergedSystematic][e.params.jetFlavor]
                    for etamin in auxMergedUncertainty:
                        if etamin>=e.params.etaMin:
                            for etamax in auxMergedUncertainty[etamin]:
                                if etamax<=e.params.etaMax:
                                    for ptmin in auxMergedUncertainty[etamin][etamax]:
                                        if ptmin>=e.params.ptMin:
                                            for ptmax in auxMergedUncertainty[etamin][etamax][ptmin]:
                                                if ptmax<=e.params.ptMax:

                                                    for variation in [ 'up_', 'down_' ]:

                                                        params = ROOT.BTagEntry.Parameters( e.params.operatingPoint, 
                                                                                            e.params.measurementType, 
                                                                                            variation + mergedSystematic,
                                                                                            e.params.jetFlavor,
                                                                                            etamin,
                                                                                            etamax, 
                                                                                            ptmin, 
                                                                                            ptmax,
                                                                                            e.params.discrMin,
                                                                                            e.params.discrMax )

                                                        systematicValue = str(math.sqrt(auxMergedUncertainty[etamin][etamax][ptmin][ptmax][e.params.discrMin][e.params.discrMax]))
                                                        sign = '+' if (variation=='up_') else '-'
                                                        SFFun = ROOT.TF1 ('SFFun', e.formula + sign + systematicValue, ptmin, ptmax)
                                                        entry = ROOT.BTagEntry(SFFun, params)
                                                        calib.addEntry(entry)

    with open(outputfilename, 'w') as f:
        f.write(calib.makeCSV())

                
