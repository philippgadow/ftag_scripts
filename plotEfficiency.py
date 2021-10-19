import h5py
import os
import yaml
import logging
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from glob import glob
from argparse import ArgumentParser
from pprint import pprint
from tqdm import tqdm


def getArgumentParser():
    parser = ArgumentParser()
    parser.add_argument('--ntuple_path', default='/nfs/dust/atlas/user/pgadow/ftag/data/ntuple_links/', help='path to ntuples')
    parser.add_argument('--ntuple_file_pattern', default='user.alfroch.410470.btagTraining.e6337_s3126_r10201_p3985.EMPFlow.2021-09-07-T122808-R14883_output.h5/*.h5', help='ntuples pattern')
    parser.add_argument('--n_jets', default=5e5)
    parser.add_argument('--binning', default=os.path.join('data', 'binning.yml'))
    parser.add_argument('--compare_cdi', default=os.path.join('data', 'cdi_dump.csv'))
    return parser


def get_jets(filename, n_jets=-1, jet_flavour=-1):
    """Helper function to extract jet and track information from a h5 ntuple.

    :param filename: path to the h5 ntuple
    :param n_jets: number of jets to be extracted from the ntuple
    :returns: (jets), where jets is a numpy array of jets.
    """
    logging.info("Opening file " + filename)
    data_set = h5py.File(filename, "r")
    jets = data_set["jets"]
    if jet_flavour > 0:
        jets = jets[jets["HadronConeExclTruthLabelID"] == int(jet_flavour)]
    logging.info(f"Total number of jets in file: {jets.size}")
    if n_jets > 0:
        jets = jets[:n_jets]
    return jets


def GetScore(pb, pc, pu, ptau=None, fc=0.018, ftau=None):
    pb = pb.astype("float64")
    pc = pc.astype("float64")
    pu = pu.astype("float64")
    add_small = 1e-10
    if ptau is not None:
        if ftau is None:
            flight = 1 - fc
            ftau = flight
        else:
            flight = 1 - fc - ftau
        ptau = ptau.astype("float64")
        return np.log(
            (pb + add_small)
            / (flight * pu + ftau * ptau + fc * pc + add_small)
        )
    return np.log((pb + add_small) / ((1.0 - fc) * pu + fc * pc + add_small))


def computeEfficiency(jets, cut):
    pb  = jets["DL1r_pb"]
    pc = jets["DL1r_pc"]
    pu = jets["DL1r_pu"]
    d_b = GetScore(pb, pc, pu)
    try:
        return float((jets[d_b > cut]).size) / float(jets.size)
    except ZeroDivisionError:
        return -1


def plotEfficiency(eff_dict, eff_wp, label, variable, bins, cdi_data=None):
    fig, ax = plt.subplots()
    x = [float(k) for k in eff_dict.keys()]
    y = [float (v) for v in eff_dict.values()]
    plt.plot(x, y, 'o')
    if cdi_data is not None:
        x_cdi = cdi_data[variable].values
        y_cdi = cdi_data['eff'].values
        plt.plot(x_cdi, y_cdi, 'o')

    plt.xlabel(variable)
    plt.ylabel('Efficiency')
    plt.title(f"{label}, {eff_wp}% eff. working point")
    ax.set_ylim([0.5, 0.9])
    fig.savefig(f"eff_{variable}_{int(eff_wp)}.png")


def plotEfficiency2D(eff, eff_wp, label, var_x, var_y, binning_x, binning_y, ratio=False):
    fig, ax = plt.subplots()
    if ratio:
        cmap = sns.color_palette("vlag", as_cmap=True)
        sns.heatmap(eff, linewidths=.5, cmap=cmap, vmin=.5, vmax=1.5, annot=True, fmt=".2f", )
    else:
        cmap = sns.color_palette("viridis", as_cmap=True)
        sns.heatmap(eff, linewidths=.5, cmap=cmap, vmin=0., vmax=1., annot=True, fmt=".2f", )
    plt.xlabel(var_x)
    plt.ylabel(var_y)
    fig.savefig(f'eff_{eff_wp}_{label}.png')


def main():
    args = getArgumentParser().parse_args()
    # list of h5 ntuples
    ntuples = glob(os.path.join(args.ntuple_path, args.ntuple_file_pattern))
    jets = None

    logging.info("Processing ntuples...")
    pbar = tqdm(total=args.n_jets)
    for i, filename in enumerate(ntuples):
        if jets is None:
            # iteration over first ntuple
            jets = get_jets(filename, jet_flavour=5)
            pbar.update(jets.size)
        else:
            # iterations over second and following ntuples
            add_jets = get_jets(filename, jet_flavour=5)
            pbar.update(add_jets.size)
            jets = np.concatenate([jets, add_jets])
            del add_jets
        if args.n_jets < jets.size: break
    pbar.close()

    # compute efficiencies
    cut_values = {
        '60': 4.565,
        '70': 3.245,
        '77': 2.195,
        '85': 0.665,
    }

    # define binning
    with open(args.binning, "r") as f:
        binning = yaml.safe_load(f)
        binning_pt = binning['PFlow_DL1r_B']['pt']
        binning_eta = binning['PFlow_DL1r_B']['abseta']

    # get reference from CDI file
    columns = ['tagger', 'jet', 'wp', 'flavour', 'pt', 'eta', 'eff']
    data_cdi = pd.read_csv(
        args.compare_cdi,
        names=columns)
    data_cdi = data_cdi[(data_cdi['tagger'] == 'DL1r') &\
                        (data_cdi['jet'] == 'AntiKt4EMPFlowJets_BTagging201903') &\
                        (data_cdi['flavour'] == 'B')]

    pt = jets["pt_btagJes"] * 0.001  # GeV
    abs_eta = jets["absEta_btagJes"]

    for eff_wp, cut in cut_values.items():
        efficiency = computeEfficiency(jets, cut)
        print(eff_wp, efficiency)

        # # efficiency vs pt
        # eff_dict_pt = {}
        # for pt_low, pt_high in zip(binning_pt[0:-1], binning_pt[1:]):
        #     selected_jets = jets[(pt > float(pt_low)) & (pt <= float(pt_high))]
        #     efficiency = computeEfficiency(selected_jets, cut)
        #     pt_center = f"{(pt_low + pt_high) * 0.5:0.2f}"
        #     eff_dict_pt[pt_center] = efficiency
        #     print(pt_low, pt_high, pt_center, eff_wp, efficiency)
        # selected_data_cdi = data_cdi[(data_cdi['wp'] == f'FixedCutBEff_{int(eff_wp)}') & (data_cdi['eta'] == 0.30)][['pt', 'eff']]
        # plotEfficiency(eff_dict_pt, eff_wp, 'DL1r b-jets (ttbar)', 'pt', binning_pt, selected_data_cdi)

        # # efficiency vs eta
        # eff_dict_eta = {}
        # for eta_low, eta_high in zip(binning_eta[0:-1], binning_eta[1:]):
        #     selected_jets = jets[(abs_eta > float(eta_low)) & (abs_eta <= float(eta_high))]
        #     efficiency = computeEfficiency(selected_jets, cut)
        #     eta_center = f"{(eta_low + eta_high) * 0.5:0.2f}"
        #     eff_dict_eta[eta_center] = efficiency
        #     print(eta_low, eta_high, eta_center, eff_wp, efficiency)
        # plotEfficiency(eff_dict_eta, eff_wp, 'DL1r b-jets (ttbar)', 'eta', binning_eta)

        # 2d efficiency map
        eff_dict = {}
        for pt_low, pt_high in zip(binning_pt[0:-1], binning_pt[1:]):
            pt_center = f"{(pt_low + pt_high) * 0.5:0.2f}"
            eff_dict[pt_center] = {}
            for eta_low, eta_high in zip(binning_eta[0:-1], binning_eta[1:]):
                selected_jets = jets[(abs_eta > float(eta_low)) & (abs_eta <= float(eta_high)) & (pt > float(pt_low)) & (pt <= float(pt_high))]
                efficiency = computeEfficiency(selected_jets, cut)
                eta_center = f"{(eta_low + eta_high) * 0.5:0.2f}"
                eff_dict[pt_center][eta_center] = efficiency
        
        selected_data_eff = pd.DataFrame.from_dict(eff_dict)
        selected_data_cdi = (data_cdi[(data_cdi['wp'] == f'FixedCutBEff_{int(eff_wp)}')][['pt', 'eta', 'eff']]).pivot(index='eta', columns='pt', values='eff')

        print(selected_data_eff)
        print(selected_data_cdi)

        selected_data_ratio = selected_data_eff / selected_data_cdi

        print(selected_data_ratio)

        plotEfficiency2D(selected_data_eff, eff_wp, 'dl1r_b_ttbar_tdd', 'pt', 'eta', binning_pt, binning_eta)
        plotEfficiency2D(selected_data_cdi, eff_wp, 'dl1r_b_ttbar_cdi', 'pt', 'eta', binning_pt, binning_eta)

        plotEfficiency2D(selected_data_ratio, eff_wp, 'dl1r_b_ttbar_ratio', 'pt', 'eta', binning_pt, binning_eta, ratio=True)


if __name__ == '__main__':
    main()
