from __future__ import division, print_function, unicode_literals
import numpy as np
import sys
import ROOT
from collections import OrderedDict
import pandas as pd
import os
import matplotlib.pyplot as plt
import pdb

typeDict = {"Char_t": 'S16', 'Double_t': 'f8', 'Int_t': 'i8'}

def load_tuple_to_array(filename, tuplename, directory=None):
    """
    load_tuple_to_array(filename, tuplename, directory=None)
    Return root tuple as numpy structured array
    filename: string
    tuplename: string
    directory: string (name of directory that holds the tuple)
    """
    f = ROOT.TFile(filename)
    if directory is not None:
        f.cd(directory)

    tree = ROOT.gDirectory.Get(tuplename)
    names = [x.GetName() for x in tree.GetListOfBranches()]
    types = [tree.GetLeaf(n).GetTypeName() for n in names]

    dt = np.dtype([(n, typeDict[t]) for (n,t) in zip(names, types)])

    n_entries = x.GetEntries()
    lists = {n: [] for n in names}

    # Populate the lists
    for ii in xrange(0, n_entries):
        if (ii//1000000)*1000000 == ii:
            print("> Read " + str(ii) + " of " + str(n_entries) + " values.")
        tree.GetEntry(ii)
        for (n, l) in lists.items():
            if dt[n] == np.dtype("S16"):
                l.append(tree.GetLeaf(n).GetValueString())
            else:
                l.append(tree.GetLeaf(n).GetValue())

    # list as array
    arrlist = zip(*[lists[n] for n in names])
    print("> Building array of " + str(n_entries) + " elements")
    return np.array(arrlist, dtype=dt)


def get_damage_summary(rootfile):
    """Get a dictionary that summarises damage for the given root file

    dict = get_damage_summary(rootfile)

    Return the dictionary:("Complexity", "Source", "IndirectEvents", "SSB/DSB",
        "IndirectDSB/TotalDSB", "Indirect/Total", "FragmentHisto", "BaseDamage",
        "StrandDamage", "IndirectHits", "DirectHits")
    """
    complexity = OrderedDict([("None", 0), ("SSB", 0), ("SSB+", 0),
                              ("2SSB", 0), ("DSB", 0), ("DSB+", 0),
                              ("DSB++", 0)])
    source = OrderedDict([("None", 0), ("SSBi", 0), ("SSBd", 0), ("SSBm", 0),
                          ("DSBi", 0), ("DSBd", 0), ("DSBm", 0)])
    indirect_events = OrderedDict([("base_oh", 0), ("base_eaq", 0),
                                   ("base_h", 0), ("strand_oh", 0),
                                   ("strand_eaq", 0), ("strand_h", 0),])
    output = OrderedDict([("Complexity", complexity),
                          ("Source", source),
                          ("IndirectEvents", indirect_events),
                          ("SSB/DSB", 0),
                          ("IndirectDSB/TotalDSB", 0),
                          ("Indirect/Total", 0),
                          ("FragmentHisto", None),
                          ("BaseDamage", 0),
                          ("StrandDamage", 0),
                          ("IndirectHits", 0),
                          ("DirectHits", 0)])

    # read the tuple, fill the output dictionary
    # could probably do this more elegantly
    root_tuple = load_tuple_to_array(rootfile, "damage", "tuples")
    if "TypeClassification" in root_tuple.dtype.fields:
        for val in root_tuple["TypeClassification"]:
            for key in output["Complexity"].keys():
                if val == key:
                    output["Complexity"][key] += 1
    if "SourceClassification" in root_tuple.dtype.fields:
        for val in root_tuple["SourceClassification"]:
            for key in output["Source"].keys():
                if val == key:
                    output["Source"][key] += 1

    if "OHBaseHits" in root_tuple.dtype.fields:
        for val in root_tuple["OHBaseHits"]:
            output["IndirectEvents"]["base_oh"] += val
        for val in root_tuple["OHStrandHits"]:
            output["IndirectEvents"]["strand_oh"] += val
        for val in root_tuple["EaqBaseHits"]:
            output["IndirectEvents"]["base_eaq"] += val
        for val in root_tuple["EaqStrandHits"]:
            output["IndirectEvents"]["strand_eaq"] += val
        for val in root_tuple["HBaseHits"]:
            output["IndirectEvents"]["base_h"] += val
        for val in root_tuple["HStrandHits"]:
            output["IndirectEvents"]["strand_h"] += val

    if "IndirectHits" in root_tuple.dtype.fields:
        for val in root_tuple["IndirectHits"]:
            output["IndirectHits"] += val

    if "DirectHits" in root_tuple.dtype.fields:
        for val in root_tuple["DirectHits"]:
            output["DirectHits"] += val

    if "BaseDamage" in root_tuple.dtype.fields:
        for val in root_tuple["BaseDamage"]:
            output["BaseDamage"] += val

    if "StrandDamage" in root_tuple.dtype.fields:
        for val in root_tuple["StrandDamage"]:
            output["StrandDamage"] += val

    cl = output["Complexity"]
    sr = output["Source"]
    ndsb = float(cl["DSB"] + cl["DSB+"] + cl["DSB++"])
    nssb = float(cl["SSB"] + cl["2SSB"] + cl["SSB+"])
    ndsb_indirect = float(sr["DSBi"] + sr["DSBm"])
    nssb_indirect = float(sr["SSBi"] + sr["SSBm"])
    output["SSB/DSB"] = nssb/ndsb if ndsb > 0 else -999
    output["IndirectDSB/TotalDSB"] = ndsb_indirect/ndsb if ndsb > 0 else -999
    output["Indirect/Total"] = (ndsb_indirect + nssb_indirect)/(nssb + ndsb)
    fragment_bins, fragments = np.histogram(root_tuple["FragmentLength"],
                                            bins=216, range=[0, 215])
    output["FragmentHisto"] = (fragment_bins, fragments)
    return output

def print_summary():
    files = [("e-", 300, "6A", 17.5, 17.5, "direct_electrons/cyl_e_300ev.root"),  # NOQA
             ("e-", 1000, "6A", 17.5, 17.5, "direct_electrons/cyl_e_1000ev.root"),  # NOQA
            #  ("e-", 3000, "6A", 17.5, 17.5, "direct_electrons/cyl_e_3000ev.root"),  # NOQA
             ("e-", 4500, "6A", 17.5, 17.5, "direct_electrons/cyl_e_4500ev.root"),  # NOQA
             ("e-", 300, "5A", 17.5, 17.5, "direct_electrons/cyl_e_300ev_n10M_r5A.root"),  # NOQA
             ("e-", 300, "6A", 17.5, 17.5, "direct_electrons/cyl_e_300ev_n10M_r6A.root"),  # NOQA
             ("e-", 300, "7A", 17.5, 17.5, "direct_electrons/cyl_e_300ev_n10M_r7A.root"),  # NOQA
             ("proton", 4500, "6A", 17.5, 17.5, "direct_protons/cyl_p_4500ev.root"),  # NOQA
             ("proton", 10000, "6A", 17.5, 17.5, "direct_protons/cyl_p_10kev.root"),  # NOQA
             ("proton", 15000, "6A", 17.5, 17.5, "direct_protons/cyl_p_15kev.root"),  # NOQA
             ("proton", 20000, "6A", 17.5, 17.5, "direct_protons/cyl_p_20kev.root"),  # NOQA
             ("proton", 4500, "4A", 17.5, 17.5, "direct_radius/cyl4A.root"),
             ("proton", 4500, "5A", 17.5, 17.5, "direct_radius/cyl5A.root"),
             ("proton", 4500, "6A", 17.5, 17.5, "direct_radius/cyl6A.root"),
             ("proton", 4500, "7A", 17.5, 17.5, "direct_radius/cyl7A.root"),
             ("proton", 4500, "8A", 17.5, 17.5, "direct_radius/cyl8A.root"),
             ("proton", 4500, "9A", 17.5, 17.5, "direct_radius/cyl9A.root"),
             ("proton", 4500, "6A", 8.22, 8.22, "direct_energy/cyl8eV.root"),
             ("proton", 4500, "6A", 11.0, 11.0, "direct_energy/cyl11eV.root"),
             ("proton", 4500, "6A", 14.0, 14.0, "direct_energy/cyl14eV.root"),
             ("proton", 4500, "6A", 17.5, 17.5, "direct_energy/cyl17eV.root"),
             ("proton", 4500, "6A", 20.0, 20.0, "direct_energy/cyl20eV.root"),
             ("proton", 4500, "6A", 24.0, 24.0, "direct_energy/cyl24eV.root"),
             ("proton", 4500, "6A", 5.00, 38.0, "direct_energy/cylPARTRAC.root")]  # NOQA

    analysed = []
    classifications = ["None", "SSB", "SSB+", "2SSB", "DSB", "DSB+", "DSB++"]
    for f in files:
        classified = {key: 0 for key in classifications}
        damage = load_tuple_to_array(f[5], "damage", "tuples")
        for val in damage["Classification"]:
            for key in classified.keys():
                if val == key:
                    classified[key] += 1
        analysed.append( (f, classified) )

    print("Primary  En (eV)  Rad  E1 (eV)  E2 (eV)  None   SSB  SSB+  2SSB   DSB  DSB+  DSB++  SSB/DSB")  # NOQA
    print("-------------------------------------------------------------------------------------------")  # NOQA
    for line in analysed:
        f = line[0]
        cl = line[1]
        ndsb = float(cl["DSB"] + cl["DSB+"] + cl["DSB++"])
        nssb = float(cl["SSB"] + cl["2SSB"] + cl["SSB+"])
        if ndsb == 0:
            ndsb = -1
        ratio = nssb/ndsb
        print("{0[0]:7s}  {0[1]:7d}  {0[2]:3s}  {0[3]:7.2f}  {0[4]:7.2f}  {1[None]:4d}  {1[SSB]:4d}  {1[SSB+]:4d}  {1[2SSB]:4d}  {1[DSB]:4d}  {1[DSB+]:4d}   {1[DSB++]:4d}  {2:7.4f}".format(f, cl, ratio))  # NOQA
    return analysed


def get_parameters(indices, directory):
    """Get input parameters for a given list of indices

    out = get_parameters(indices, directory)

    out: pandas dataframe containing desired indices
    dir: directory to look for the params.txt file
    """
    filename = os.path.join(directory, "params.txt")
    assert os.path.exists(filename),\
        "{} does not contain a params.txt file".format(directory)
    tab = pd.read_table(filename, delimiter=" ")
    return tab[tab.FILENUM.isin(indices)]


def do_damage_distance_graph(indices, directory):
    """Make a graph of Indirect Damage
    """
    assert os.path.exists("figs"), "Please make a directory called figs to" +\
        " save the figure"
    params = get_parameters(indices, directory)
    vals = [(n, d,
             get_damage_summary(os.path.join(directory,"{}.root".format(n))))
             for (n, d) in zip(params.FILENUM, params.KILLDIST)]
    x = [v[1] for v in vals]
    dsb = [v[2]["Source"]["DSBi"] +
           v[2]["Source"]["DSBm"] for v in vals]
    ssb = [v[2]["Source"]["SSBi"] +
           v[2]["Source"]["SSBm"] for v in vals]
    fig = plt.figure(figsize=[3, 2])
    ax = fig.add_subplot(111)
    ax.set_xlabel("Distance (nm)")
    ax.set_ylabel("Indirect Breaks")
    # ax.plot(x, ssb, "r-", label="SSBs")
    ax.plot(x, dsb, "b--", label="DSBs")
    ax.legend(frameon=False, loc="upper left")
    fig.savefig("figs/{}-damage-distance.pdf".format(directory),
                bbox_inches="tight")
    return None



def print_usage():
    print("""python analysis.py [-h] directory --method indices
    -h print this help message

Example usage:
    python analysis.py cylinders --indirect-range 1000-1030 1111
        Run the indirect-range method, for root files in the cyliners directory
        with file numbers in the range [1000, 1030) and also number 1111
        Files do not need to exist

Warning:
    Figures are created in a directory called "./figs". This script will
    overwrite existing figures that have the same filename as that being
    saved.

Available Methods:
    --indirect-range: Make a figure
""")
    return None

if __name__ == "__main__":
    if "-h" in sys.argv:
        print_usage()
        sys.exit()
    directory = sys.argv[1]
    if directory[0:2] == "--":
        print("Invalid directory name")
        print_usage()
        sys.exit()
    if "--indirect-range" in sys.argv:
        start_idx = sys.argv.index("--indirect-range") + 1
        ranges = []
        for val in sys.argv[start_idx:]:
            if "-" not in val:
                ranges.append(int(val))
            else:
                [mn, mx] = val.split("-")
                for ii in range(int(mn), int(mx)):
                    ranges.append(ii)
        do_damage_distance_graph(ranges, directory)
