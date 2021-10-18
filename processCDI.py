import ROOT
import ctypes
import yaml
import numpy as np
from argparse import ArgumentParser
from pprint import pprint
from itertools import product

def getArgumentParser():
    parser = ArgumentParser()
    parser.add_argument('--cdi_file', default='/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/xAODBTaggingEfficiency/13TeV/2020-21-13TeV-MC16-CDI-2021-04-16_v1.root', help='path to CDI file')
    parser.add_argument('--binning', default='data/binning.yml')
    return parser


def processCDI(args, tagger='DL1r', jet='AntiKt4EMPFlowJets_BTagging201903', wp='FixedCutBEff_70'):
    # with uproot.open(cdi_file) as f:
    #     tags = ['B', 'C', 'T', 'Light']
    #     for tag in tags:
    #         path = f'{algo}/{jets}/{wp}/{tag}/default_Eff'
    #         pprint(f[path])

    # setup xAOD
    ROOT.xAOD.Init().ignore()
    ROOT.gROOT.SetBatch(ROOT.kTRUE)

    # settings
    flavour = 'B'
    sample_id = 410470

    # set up b-tagging efficiency tool
    tool_name = 'ftag_btaggingefficiencytool'
    btagtool = ROOT.BTaggingEfficiencyTool(tool_name)
    btagtool.setProperty('ScaleFactorFileName', args.cdi_file).ignore()
    btagtool.setProperty("TaggerName", tagger).ignore()
    btagtool.setProperty("OperatingPoint", wp).ignore()
    btagtool.setProperty("JetAuthor", jet).ignore()
    btagtool.setProperty("MinPt", 20000).ignore()
    btagtool.setProperty('Efficiency'+flavour+'Calibrations', sample_id).ignore()
    btagtool.initialize().ignore()

    # set up dummy jet variables and target efficiency float variable
    dummy_jet_vars = ROOT.Analysis.CalibrationDataVariables()
    dummy_jet_type = 5  # b-jet
    dummy_eff = ctypes.c_float(0.)

    # define binning
    with open(args.binning, "r") as f:
        binning = yaml.safe_load(f)
        binning_pt = np.array(binning['PFlow_DL1r_B']['pt'])
        binning_eta = np.array(binning['PFlow_DL1r_B']['abseta'])

    # get bin centers
    binning_pt = (binning_pt[0:-1] + binning_pt[1:]) * 0.5
    binning_eta = (binning_eta[0:-1] + binning_eta[1:]) * 0.5

    # get MC efficiency
    output = ""
    for pt, eta in product(binning_pt, binning_eta):
        dummy_jet_vars.jetPt = pt * 1000.  # MeV
        dummy_jet_vars.jetEta = eta
        btagtool.getMCEfficiency(dummy_jet_type, dummy_jet_vars, dummy_eff).ignore()
        eff = dummy_eff.value
        output += '{0},{1},{2},{3},{4},{5},{6}\n'.format(
            tagger,
            jet,
            wp,
            flavour,
            pt,
            eta,
            eff
        )
    return output


def main():
    args = getArgumentParser().parse_args()
    output = ''
    output += processCDI(args, wp='FixedCutBEff_60')
    output += processCDI(args, wp='FixedCutBEff_70')
    output += processCDI(args, wp='FixedCutBEff_77')
    output += processCDI(args, wp='FixedCutBEff_85')

    with open('data/cdi_dump.csv', 'w') as f:
        f.write(output)


if __name__ == '__main__':
    main()
