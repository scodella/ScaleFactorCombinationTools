#!/usr/bin/env python
import os
import sys
import ROOT
import math
import optparse
from array import *
import CondTools.BTau.dataLoader as dataLoader

# From https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/btv/btagSFProducer.py
supported_btagSF = {
    'csvv2' : {
        '2016' : {
            'inputFileName' : "btagSF_CSVv2_ichep2016.csv",
            'measurement_types' : {
                0 : "comb",  # b
                1 : "comb",  # c
                2 : "incl"   # light
            },
            'supported_wp' : [ "L", "M", "T", "shape_corr" ]
        },
        '2017' : {
            'inputFileName' : "CSVv2_94XSF_V2_B_F.csv",
            'measurement_types' : {
                0 : "comb",  # b
                1 : "comb",  # c
                2 : "incl"   # light
            },
            'supported_wp' : [ "L", "M", "T", "shape_corr"]
        }
    },
    'deepcsv' : {
        'Legacy2016' : {
            'inputFileName' : "DeepCSV_2016LegacySF_V1.csv",
            'measurement_types' : {
                0 : "comb-mujets",  # b
                1 : "comb-mujets",  # c
                2 : "incl"   # light
            },
            'supported_wp' : [ "L", "M", "T", "shape_corr"]
        },
        '2017' : {
            'inputFileName' : "DeepCSV_94XSF_V4_B_F.csv",
            'measurement_types' : {
                0 : "comb-mujets",  # b
                1 : "comb-mujets",  # c
                2 : "incl"   # light
            },
            'supported_wp' : [ "L", "M", "T", "shape_corr"]
        },
        '2018' : {
            'inputFileName' : "DeepCSV_102XSF_V1.csv",
            'measurement_types' : {
                0 : "comb-mujets",  # b
                1 : "comb-mujets",  # c
                2 : "incl"   # light
            },
            'supported_wp' : [ "L", "M", "T", "shape_corr"] 
        }    
    },
    'deepjet' : {
        'Legacy2016' : {
            'inputFileName' : "DeepJet_2016LegacySF_V1.csv",
            'measurement_types' : {
                0 : "comb-mujets",  # b
                1 : "comb-mujets",  # c
                2 : "incl"   # light
            },
            'supported_wp' : [ "L", "M", "T", "shape_corr"]
        },
        '2017' : {
            'inputFileName' : "DeepFlavour_94XSF_V3_B_F.csv",
            'measurement_types' : {
                0 : "comb-mujets",  # b
                1 : "comb-mujets",  # c
                2 : "incl"   # light
            },
            'supported_wp' : [ "L", "M", "T", "shape_corr"]
        },
        '2018' : {
            'inputFileName' : "DeepJet_102XSF_V1.csv",
            'measurement_types' : {
                0 : "comb-mujets",  # b
                1 : "comb-mujets",  # c
                2 : "incl"   # light
            },
            'supported_wp' : [ "L", "M", "T", "shape_corr"]
        }    
    },
    'cmva' : {
        '2016' : {
            'inputFileName' : "btagSF_cMVAv2_ichep2016.csv",
            'measurement_types' : {
                0 : "ttbar", # b
                1 : "ttbar", # c
                2 : "incl"   # light
            },
            'supported_wp' : [ "L", "M", "T", "shape_corr" ]
        }
    }
}

def plotScaleFactors(plottitle, scaleFactors, plotformat):

    ROOT.gROOT.SetBatch(ROOT.kTRUE)

    yoff = 0.85 if ('_T_' in plottitle) else 0.35
    leg = ROOT.TLegend(0.18, yoff, 0.50, yoff-0.2)
    leg.SetFillColor(ROOT.kWhite) 
    leg.SetBorderSize(0)
    leg.SetTextColor(1) 
    leg.SetTextSize(0.04)
    leg.SetTextFont(62) 
    
    #histoiter = 1 
    minX, maxX, minY, maxY = 999., -1., 999., -1.
       
    for function in scaleFactors:
 
        if 'central' in function.GetName():
           leg.AddEntry(function, function.GetName().replace('_central', ''), "l")

        minX = min(minX, function.GetXmin())
        maxX = max(maxX, function.GetXmax())
        minY = min(minY, function.GetMinimum())
        maxY = max(maxY, function.GetMaximum())

    minY = max(minY, 0.7)
    maxY = min(maxY, 1.5)

    histo = ROOT.TH1F('histo', '', 200, minX, maxX)

    histo.SetMaximum(maxY+0.1)
    histo.SetMinimum(minY-0.1)
    
    histo.SetXTitle("p_{T} [GeV]")
    histo.SetYTitle("SF_{b}")
        
    histo.SetLabelSize(0.05, "XYZ")
    histo.SetTitleSize(0.06, "XYZ") 
    histo.SetLabelFont(42, "XYZ") 
    histo.SetTitleFont(42, "XYZ")
    histo.SetTitleOffset(0.95,"X")
    histo.SetTitleOffset(0.8,"Y")
    histo.SetTickLength(0.06,"X")
    histo.SetNdivisions(509, "XYZ")
    histo.GetXaxis().SetMoreLogLabels()
    histo.GetXaxis().SetNoExponent()

    ROOT.gStyle.SetOptFit(ROOT.kFALSE)
    ROOT.gStyle.SetOptStat(ROOT.kFALSE)
    ROOT.gStyle.SetOptTitle(ROOT.kFALSE)
    
    CC = ROOT.TCanvas("CC", "", 1200, 800)

    CC.SetFillColor(10)
    CC.SetFillStyle(4000)
    CC.SetBorderSize(2)
    CC.Divide(1, 1)
    CC.Range(0,0,1,1)
    CC.SetFillColor(10)
    CC.SetBorderMode(0)
    CC.SetBorderSize(2)
    CC.SetTickx(1)
    CC.SetTicky(1)
    CC.SetLeftMargin(0.16)
    CC.SetRightMargin(0.02)
    CC.SetTopMargin(0.05)
    CC.SetBottomMargin(0.13)
    CC.SetFrameFillColor(0)
    CC.SetFrameFillStyle(0)
    CC.SetFrameBorderMode(0)
    
    pad1 = CC.GetPad(1)

    pad1.SetFillColor(0)
    pad1.SetBorderMode(0)
    pad1.SetBorderSize(2)
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
    
    pad1.cd()
    
    pad1.SetLogx()
    
    histo.DrawCopy()

    for function in scaleFactors:
        function.DrawCopy("same")

    leg.Draw()

    for ext in plotformat.split('-'):
        CC.Print("./Plots/" + plottitle + '.' + ext)

if __name__ == '__main__':

    # Input parameters
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    
    parser.add_option('--taggers'     , dest='taggers'     , help='Tagger(s) to be compared'            , default='deepcsv')
    parser.add_option('--years'       , dest='years'       , help='Year(s) to be compared'              , default='2016-2017-2018')
    parser.add_option('--inputpath'   , dest='inputpath'   , help='Path where csv files are stored '    , default='./CSVFiles')
    parser.add_option('--inputfiles'  , dest='inputfiles'  , help='Additional csv files to be compared' , default='')
    parser.add_option('--customfiles' , dest='customfiles' , help='Only use custom csv files'           , default=False, action='store_true')
    parser.add_option('--wps'         , dest='wps'         , help='Working point(s) to be compared'     , default='M')
    parser.add_option('--meastypes'   , dest='meastypes'   , help='Measurement type(s) to be compared'  , default='default')
    parser.add_option('--flavour'     , dest='flavour'     , help='Flavour to be studied'               , default='b')
    parser.add_option('--plotformat'  , dest='plotformat'  , help='Formats of the plot, e.g.: png-pdf'  , default='png')
    (opt, args) = parser.parse_args()
    
    years = opt.years.split('-')

    flavour = { 'b' : 0, 'c' : 1, 'l' : 2, '0' : 0, '1' : 1, '2' : 2 }.get(opt.flavour, None)
    if flavour==None:
        print 'Wrong choice of jet flavour:', opt.flavour, '-> exiting'
        exit()

    if opt.customfiles:
        supported_btagSF.clear()
                     
    if opt.inputfiles!='':

        for csvfile in opt.inputfiles.split(','):

            if '/' not in csvfile:
                csvfile = opt.inputpath + '/' + csvfile

            csvyear = csvfile.split('/')[-1].replace('.csv', '')
            years.append(csvyear)

            tagger = csvyear.split('_')[0].lower()

            custom_btagSF = { } 
            custom_btagSF['inputFileName'] = csvfile
            custom_btagSF['measurement_types'] = { 0 : "", 1 : "", 2 : ""}
            custom_btagSF['supported_wp'] = [ ]

            for data in dataLoader.get_data(csvfile):
                for e in data.entries:

                    wp = { 0 : "L", 1 : "M", 2 : "T" }.get(e.params.operatingPoint, None)

                    if wp not in custom_btagSF['supported_wp']:
                        custom_btagSF['supported_wp'].append(wp)

                    if e.params.measurementType not in custom_btagSF['measurement_types'][e.params.jetFlavor]:
                        custom_btagSF['measurement_types'][e.params.jetFlavor] += "-" + e.params.measurementType

            if tagger not in supported_btagSF:
                supported_btagSF[tagger] = { }

            supported_btagSF[tagger][csvyear] = custom_btagSF

    scaleFactors = [ ]

    color = 1

    for taggername in opt.taggers.split('-'):
        tagger = taggername.lower()
        if tagger in supported_btagSF:

            for campaign in supported_btagSF[tagger]:

                useThisCampaign = False
                for year in years:
                    if year in campaign:
                        useThisCampaign = True
                
                if useThisCampaign:

                    inputFileName = supported_btagSF[tagger][campaign]['inputFileName']
                    if '/' not in inputFileName:
                        inputFileName = os.path.join(opt.inputpath, inputFileName)
                    loaders = dataLoader.get_data(inputFileName)

                    for wp in supported_btagSF[tagger][campaign]['supported_wp']:
                        if wp in opt.wps:

                            wp_btv = { "l" : 0, "m" : 1, "t" : 2 }.get(wp.lower(), None)
                            if wp_btv==None:
                                print 'Working point', wp, 'not supported -> skipping'
                                continue

                            meastypes = opt.meastypes if (opt.meastypes!='default') else supported_btagSF[tagger][campaign]['measurement_types'][flavour]
                            for meastype in meastypes.split('-'):
                                if meastype!='' and meastype in supported_btagSF[tagger][campaign]['measurement_types'][flavour]:
                                    
                                    minPt, maxPt = 999999., -1.
                                    title = tagger + '_' + campaign + '_' + wp + '_' + meastype
                                    function = { 'central' : ' 0. ', 'up' : ' 0. ', 'down' : ' 0. ' }

                                    for data in loaders:
                                        for e in data.entries:

                                            if e.params.measurementType==meastype and e.params.operatingPoint==wp_btv and e.params.jetFlavor==flavour:
                                                if e.params.sysType in function:

                                                    formula = '(x>='+str(e.params.ptMin)+' && x<'+str(e.params.ptMax)+') ? '+e.formula
                                                    function[e.params.sysType] = function[e.params.sysType].replace(' 0. ',formula+' : 0. ')
                                                    minPt = min(minPt, e.params.ptMin)
                                                    maxPt = max(maxPt, e.params.ptMax)
                                                
                                    for syst in function:
                                        systFunction = ROOT.TF1(title+'_'+syst, function[syst], minPt, maxPt)
                                        width = 3 if (syst=='central') else 1
                                        style = 1 if (syst=='central') else 2
                                        systFunction.SetLineWidth(width)
                                        systFunction.SetLineStyle(style)
                                        systFunction.SetLineColor(color)
                                        scaleFactors.append(systFunction)

                                    color += 1

    if len(scaleFactors)==0:
        print 'Exiting with no scale factors to plot'
        exit()

    plottitle = opt.taggers + '_' + opt.years + '_' + opt.wps + '_' + opt.flavour
    if opt.meastypes!='default':
        plottitle += '_' + opt.meastypes 
    plotScaleFactors(plottitle, scaleFactors, opt.plotformat)
