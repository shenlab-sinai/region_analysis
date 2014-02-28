#! /usr/bin/env python

import os
import sys
import json
import string
from argparse import ArgumentParser
import regionanalysis.packageinfo
import regionanalysis.annotationdb

def read_gnlist(h_sp, installed_db, mode):
    """
    Read installed genomes list and make human readable format.
       Args:
         h_sp: vector of header.
         installed_db: list of dicts with genome information.
       Returns: (header split vector, hash table of installed genomes, 
                 vector of column widths)
    """
    v_cw = map(len, h_sp)  # column widths initialize to header widths.

    g_tbl = {}
    for rec in installed_db:
        for db_info in rec[u"databases"]:
            if db_info[u"database"] == u"ensembl":
                ensembl_ver = db_info[u"version"]
        r_sp = [rec[u"genome"], rec[u"assembly"], rec[u"species"],\
                ensembl_ver, rec[u"version"], rec[u"path"]]
        r_sp = [x.encode('ascii','ignore') for x in r_sp]
        r_cw = map(len, r_sp)
        v_cw = map(max, v_cw, r_cw)

        r_tbl = dict(zip(h_sp, r_sp))
        if(mode == "vector"):
            g_tbl[r_tbl["ID"]+r_tbl["RA_Ver"]] = r_sp
        elif(mode == "hash"):
            g_tbl[r_tbl["ID"]+r_tbl["RA_Ver"]] = r_tbl
        else:
            pass

    return (h_sp, g_tbl, v_cw)

def listgn(args):
    """
    List informations of all installed databases.
    """
    import math

    module_dir = os.path.dirname(os.path.realpath(regionanalysis.__file__))
    installed_db = regionanalysis.annotationdb.getAllInstalledDB(module_dir)
    h_sp = ["ID", "Assembly", "Species", "ENSEMBL_Ver", "RA_Ver", "Location"]
    (h_sp, g_tbl, v_cw) = read_gnlist(h_sp, installed_db, "vector")
    # Format column widths to beautify output.
    tab_u = 4
    v_cw = map(lambda x: int(math.ceil(float(x) / tab_u) * tab_u + 1), v_cw)
    # List genomes to screen.
    print "".join(map(lambda x, y: x.ljust(y), h_sp, v_cw))  # header.
    for k in sorted(g_tbl.viewkeys(), key=str.lower):
        print "".join(map(lambda x, y: x.ljust(y), g_tbl[k], v_cw))

    # print(installed_db) # for debuging

def getAns():
    """
    Get the answer "yes" or "No" from user's input.
    """
    ans = raw_input("Continue?(y/n): ")
    while True:
        if ans == 'y' or ans == 'Y' or ans == 'n' or ans == 'N':
            break
        else:
            ans = raw_input("The answer must be y/Y or n/N: ")
    return ans

def install(args):
    """
    Install databases from tar.gz file to database folder.
    """
    import tarfile
    import shutil
    pkg_file = args.pkg
    module_dir = os.path.dirname(os.path.realpath(regionanalysis.__file__))
    installed_db = regionanalysis.annotationdb.getAllInstalledDB(module_dir)
    yestoall = args.yes

    try:
        pkg_f = tarfile.open(pkg_file, "r:gz")
    except tarfile.ReadError:
        print "Read package file {0} error.".format(pkg_file),
        print "The downloaded file may be corrupted."
        sys.exit()

    print "Extracting information from package...\n",
    sys.stdout.flush()
    pkg_files = pkg_f.getnames()
    # extract package info from json.
    json_file = [x for x in pkg_files if os.path.splitext(x)[1] == ".json"][0]
    json_fp = pkg_f.extractfile(json_file)
    cur_genome_info = json.load(json_fp)
    json_fp.close()
    # check if the genome is in the the databases.
    for genome_info in installed_db:
        (location, genome_name) = os.path.split(genome_info["path"])
        if genome_name == pkg_files[0]:
            sys.stderr.write("%s, RAver %s already installed at %s!\n"%(genome_info["genome"], genome_info["version"], genome_info["path"]))
            if yestoall == True:
                sys.stderr.write("The installed database will be removed!\n")
                shutil.rmtree(genome_info["path"])
            else:
                sys.stderr.write("Would you want to continue the installation?\n")
                ans = getAns()
                if ans == 'y' or ans == 'Y':
                    sys.stderr.write("The installed database will be removed!\n")
                    shutil.rmtree(genome_info["path"])
                else:
                    sys.stderr.write("The installation is cancelled!")
                    sys.exit()
    install_path = regionanalysis.annotationdb.getInstallPath(module_dir)
    try:
        if not os.path.isdir(install_path):
            os.mkdir(install_path)
        sys.stdout.write("Installing %s, RAver %s in %s...\n"%(cur_genome_info["genome"], cur_genome_info["version"], install_path))
        pkg_f.extractall(install_path)
    except tarfile.ExtractError:
        print "Extract files from package error.", 
        print "The downloaded file may be corrupted."

def remove(args):
    import shutil
    cur_genome = args.gn
    module_dir = os.path.dirname(os.path.realpath(regionanalysis.__file__))
    installed_db = regionanalysis.annotationdb.getAllInstalledDB(module_dir)
    yestoall = args.yes
    for genome_info in installed_db:
        do_rm = False
        if (genome_info["genome"] == cur_genome):
            sys.stdout.write("%s, RAver %s in the databases will be removed.\n"%(cur_genome, genome_info["version"]))
            if yestoall:
                do_rm = True
            else:
                ans = getAns()
                if ans == 'y' or ans == 'Y':
                    do_rm = True
        if do_rm:
            sys.stdout.write("Removing genome...\n")
            sys.stdout.flush()
            shutil.rmtree(genome_info["path"])

def main():
    opt_parser = ArgumentParser(
        description="Manage annotation databases of region_analysis.", \
        prog="region_analysis_db.py")
    subparsers = opt_parser.add_subparsers(title="Subcommands",
                                           help="additional help")

    # list parser.
    parser_list = subparsers.add_parser("list", help="List genomes installed \
                                                      in database")
    parser_list.set_defaults(func=listgn)

    # install parser.
    parser_install = subparsers.add_parser("install",
                                           help="Install genome from tar.gz \
                                                  package file")
    parser_install.add_argument("pkg", help="Package file(.tar.gz) to install",
                                type=str)
    parser_install.set_defaults(func=install)
    parser_install.add_argument("-y", "--yes", help="Say yes to all prompted questions",
                                action="store_true")

    # remove parser.
    parser_remove = subparsers.add_parser("remove",
                                          help="Remove genome from database")
    parser_remove.add_argument("gn", help="Name of genome to be \
                                           removed(e.g. hg19)", type=str)
    parser_remove.set_defaults(func=remove)
    parser_remove.add_argument("-y", "--yes", help="Say yes to all prompted questions",
                        action="store_true")

    args = opt_parser.parse_args()
    args.func(args)


#-------------------------------------------------------------------------
if __name__ == '__main__':
    main()
#-------------------------------------------------------------------------
# EOF
