from __future__ import division, print_function, unicode_literals
import os, shutil, sys, pdb


def read_param_file(fname):
    """[dict1, dict2, ...] = read_param_file(filename)

    Reads a parameter file with filename, and returns a list of dictionaries
        for the entries.
    The returned dictionaries are guaranteed to contain the key ROOTFILE,
        which is the absolute path that this parameter set corresponds to.
    """
    with open(fname, "r") as f:
        lines = [line.split() for line in f.readlines()]

    dicts = []
    idx = 0
    while idx < len(lines):
        if (lines[idx][0] == lines[idx+1][0]) and\
                (lines[idx][1] != lines[idx+1][1]):
            print("Found a parameter set for {}".format(lines[idx][0]))
            abspath = os.path.split(os.path.abspath(fname))[0]
            r_idx = lines[idx][0]
            rootfile = os.path.join(abspath, "{}.root".format(r_idx))
            if os.path.exists(rootfile):
                dicts.append(
                    {k: v for (k, v) in zip(lines[idx][1:], lines[idx+1][1:])})
                dicts[-1]["ROOTFILE"] = rootfile
                dicts[-1]["FILENUM"] = lines[idx][0]
            else:
                print("But {} does not exist".format(rootfile))
            idx += 2
        else:
            idx += 1
    return dicts


def find_param_files(path):
    """find_param_files(path)

    search for parameter files in path, matching the regex build_.*/params.py
    """
    path = os.path.abspath(path)
    print("Searching in: {}".format(path))
    param_files = []
    files = os.listdir(path)
    for f in files:
        if f[0:6] == "build_":
            p = os.path.join(path, f, 'params.txt')
            if os.path.exists(p):
                print("Found {}".format(p))
                param_files.append(p)

    return param_files


def import_param_files(paramlist, destination_dir):
    """
    """
    if paramlist == []:
        print("Nothing to import")
        return None
    # Check Destination for local params.txt
    oldparams = []
    if os.path.exists(os.path.join(destination_dir, "params.txt")):
        print("Reading old parameter file")
        with open(os.path.join(destination_dir, "params.txt")) as f:
            keys = f.readline().split()
            for line in f:
	        oldparams.append({k: v for k, v in
			          zip(keys, line.split())})
        os.rename(os.path.join(destination_dir, "params.txt"),
                  os.path.join(destination_dir, "params.txt.old"))
    
    # make a set of all parameters
    all_params = []
    for item in paramlist + oldparams:
	all_params += [k for k in item.keys()] 
    all_params = set(all_params)

    # fill an array of filenames to copy, and parameter sets to write
    to_copy = []
    to_write = oldparams[:]
    existing = [item["FILENUM"] for item in oldparams]

    for paramset in paramlist:
        if paramset["FILENUM"] not in existing:
            print("Setting {} to copy".format(paramset["ROOTFILE"]))
            to_write.append(paramset)
            to_copy.append(paramset["ROOTFILE"])

    # flesh out paramset so that all files have shared parameters
    for paramset in to_write:
        for p in all_params:
            if p not in paramset.keys():
                paramset[p] = -999

    all_params = list(all_params)
    all_params.sort()
    all_params.remove("FILENUM")
    all_params = ["FILENUM"] + all_params
    to_write.sort(key=lambda x: x["FILENUM"])
    # write params.txt file
    with open(os.path.join(destination_dir, "params.txt"), 'w') as f:
        f.write(" ".join(all_params) + "\n")
        for out in to_write:
            print("Writing {}".format(out["FILENUM"]))
            f.write(" ".join([str(out[k]) for k in all_params]) + "\n")
    
    for f in to_copy:
        print("Copying {}".format(f))
        shutil.copy2(f, destination_dir)

    return None


def print_usage():
    print("""load_files.py source_directory destination_directory

    Search for simulation root files in a source directory and move them to
    a destination directory

    source_directory: path to be search for build_*/params.py
    destination_directory: path to copy root files and build parameter dict
""")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print_usage()
        sys.exit()
    param_files = find_param_files(sys.argv[1])
    params = []
    for param_file in param_files:
        params = params + read_param_file(param_file)
    if len(set([p["FILENUM"] for p in params])) == len(params):
        import_param_files(params, sys.argv[2])
        sys.exit()
    else:
        print("There are duplicated root files in the source directory")
        print("Exiting")
        sys.exit()
