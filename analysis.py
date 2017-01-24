from __future__ import division, print_function, unicode_literals
import numpy as np
import sys
import ROOT
from collections import OrderedDict
import pandas as pd
import os
import matplotlib.pyplot as plt
import pdb
import matplotlib as mpl

mpl.rc("axes", titlesize=9, labelsize=8)
mpl.rc("legend", fontsize=8)


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


def do_damage_distance_graph(indices, directory, outfile):
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
    ax.set_ylabel("Relative Indirect Breaks")
    ssb = np.asarray(ssb)
    dsb = np.asarray(dsb)
    ax.plot(x, ssb/max(ssb), "r-", label="SSBs")
    ax.plot(x, dsb/max(dsb), "b--", label="DSBs")
    ax.legend(frameon=False, loc="upper left")
    fig.savefig(outfile, bbox_inches="tight")
    with open(outfile + ".log", 'w') as f:
        f.write("Max DSB: {}; Max SSB: {}".format(max(dsb), max(ssb)))
    return None


def make_damage_table(indices, directory):
    """Make a graph of Indirect Damage
    """
    assert os.path.exists("figs"), "Please make a directory called figs to" +\
        " save the figure"
    params = get_parameters(indices, directory)
    # vals = [(n, ein, primary, events*ein, d, elow, ehigh, physics,
    #          get_damage_summary(os.path.join(directory,"{}.root".format(n))))
    #          for (n, ein, primary, d, events, elow, ehigh, physics) in
    #          zip(params.FILENUM, params.EN_EV, params.PRIMARY,
    #              params.DMG_DIST, params.EVENTS,
    #              params.DMG_EN_LOW, params.DMG_EN_HIGH, params.PHYSICS)]
    dmgs = [get_damage_summary(os.path.join(directory,"{}.root".format(n)))
            for n in params.FILENUM]
    valdict = OrderedDict([("FILENUM", [n for n in params.FILENUM])])
    valdict["SSB/DSB"] =\
        [dmg["SSB/DSB"] for dmg in dmgs]
    valdict["IndirectDSB/TotalDSB"] =\
        [dmg["IndirectDSB/TotalDSB"] for dmg in dmgs]
    valdict["Indirect/Total"] =\
        [dmg["Indirect/Total"] for dmg in dmgs]
    # valdict["FragmentHisto"] =\
    #     [dmg["FragmentHisto"] for dmg in dmgs]
    valdict["BaseDamage"] =\
        [dmg["BaseDamage"] for dmg in dmgs]
    valdict["StrandDamage"] =\
        [dmg["StrandDamage"] for dmg in dmgs]
    valdict["IndirectHits"] =\
        [dmg["IndirectHits"] for dmg in dmgs]
    valdict["DirectHits"] =\
        [dmg["DirectHits"] for dmg in dmgs]
    for key in dmgs[0]["Source"].keys():
        valdict[key] = [dmg["Source"][key] for dmg in dmgs]
    for key in dmgs[0]["Complexity"].keys():
        valdict[key] = [dmg["Complexity"][key] for dmg in dmgs]
    for key in dmgs[0]["IndirectEvents"].keys():
        valdict[key] = [dmg["IndirectEvents"][key] for dmg in dmgs]
    df = pd.DataFrame(valdict)

    table = pd.merge(params, df, how="inner", on="FILENUM")
    return table


def save_damage_table(indices, directory, output):
    """Save a damage table
    """
    table = make_damage_table(indices, directory)
    table.to_csv(outfile)
    return None


def output_direct_damage(indices, directory, output):
    """Prepare formatted output of direct damage
    """
    table = make_damage_table(indices, directory)
    table = table[table["CHEM"] == False]
    phys4 = table[table["PHYSICS"] == 4]
    phys0 = table[table["PHYSICS"] == 0]
    columns = ["PRIMARY", "EN_EV", "EVENTS", "DMG_DIST", "DMG_EN_LOW",
               "DMG_EN_HIGH", "None", "SSB", "SSB+", "2SSB", "DSB", "DSB+",
               "DSB++", "SSB/DSB"]

    l1 = "Primary   E (eV) N       Rad  E1(eV) E2(eV) None   SSB   SSB+   2SSB    DSB   DSB+   DSB++   SSB/DSB"
    l2 = "----------------------------------------------------------------------------------------------------"
    for tab, txt in [[phys4, "Physics4"], [phys0, "Physics0"]]:
        print(txt)
        for rad in [3, 4, 5, 6, 7, 8, 9]:
            print(l1, "\n", l2)
            for elow in [12.6, 15.0, 17.5, 21.1, 30.0, 5.0]:
                for ein in [300, 500, 1000, 3000, 4500]:
                    row = tab[(tab["EN_EV"] == ein) &
                              (tab["DMG_DIST"] == rad) &
                              (tab["DMG_EN_LOW"] == elow)]
                    if len(row) > 0:
                        vals = [row.loc[row.index[0]][col] for col in columns]
                        valsfmt = [str(v) if type(v) != np.float64 else "{:>6.2f}".format(v) for v in vals]
                        valsfmt = ["{:>6}".format(v) if len(v) < 6 else v for v in valsfmt]
                        # pdb.set_trace()
                        try:
                            print(" ".join(valsfmt))
                        except:
                            pdb.set_trace()
                            print(row, file=sys.stderr)

    return None

def print_usage():
    print("""python analysis.py [-h] directory --method --indices indices
    -h print this help message

Example usage:
    python analysis.py dir --indirect-range out.pdf --indices 1000-1030 1111
        Run the indirect-range method, for root files in the cyliners directory
        with file numbers in the range [1000, 1030) and also number 1111
        Files do not need to exist

Warning:
    Created files will overwrite existing Files


Available Methods:
    --indirect-range fname: Make a figure of kill distance w. SSBs/DSBs
    --data-table fname: Make a data table of damage data
""")
    return None

if __name__ == "__main__":
    available_methods = ["--indirect-range", "--data-table", "--direct-damage"]
    if "-h" in sys.argv:
        print_usage()
        sys.exit()
    if "--indices" not in sys.argv:
        print("No indices specified")
        print_usage()
        sys.exit()
    directory = sys.argv[1]
    if directory[0:2] == "--":
        print("Invalid directory name")
        print_usage()
        sys.exit()
    method = sys.argv[2]
    if method not in available_methods:
        print("Invalid method specified")
        print_usage()
        sys.exit()
    start_idx = sys.argv.index("--indices") + 1
    ranges = []
    for val in sys.argv[start_idx:]:
        if "-" not in val:
            ranges.append(int(val))
        else:
            [mn, mx] = val.split("-")
            for ii in range(int(mn), int(mx)):
                ranges.append(ii)
    outfile = sys.argv[3]
    if outfile == "--indices":
        print("No destination file set")
        print_usage()
        sys.exit()
    if method == "--indirect-range":
        do_damage_distance_graph(ranges, directory, outfile)
        sys.exit()
    elif method == "--data-table":
        save_damage_table(ranges, directory, outfile)
        sys.exit()
    elif method == "--direct-damage":
        output_direct_damage(ranges, directory, outfile)
        sys.exit()
    else:
        print("Method not implemented")
        sys.exit()
