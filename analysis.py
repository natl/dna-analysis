import numpy as np
import ROOT

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


if __name__ == "__main__":
    print_summary()
