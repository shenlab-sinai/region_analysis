import os
import sys
import json
import glob

def loadJSON(json_file):
    fp = open(json_file)
    genome_info = json.load(fp)
    fp.close()
    genome_info[u"file"] = os.path.basename(json_file)
    genome_info[u"path"] = os.path.dirname(json_file)
    return genome_info

def expandOsPath(path):
    """
    To expand the path with shell variables.
    Arguments:
    - `path`: path string
    """
    return os.path.expanduser(os.path.expandvars(path))

def getAllPath(module_dir):
    search_path = []
    try:
        environ_path = expandOsPath(os.environ.get("RA_DB_PATH"))
    except:
        environ_path = None
    if environ_path != None:
        search_path.extend(glob.glob(os.path.join(environ_path,"*/")))
    home_dir = expandOsPath("~/.config/regionanalysis/")
    if os.path.isdir(home_dir):
        search_path.extend(glob.glob(os.path.join(home_dir, "*/")))
    search_path.extend(glob.glob(os.path.join(module_dir, "database/*/")))
    return search_path

def getInstallPath(module_dir):
    try:
        environ_path = expandOsPath(os.environ.get("RA_DB_PATH"))
    except:
        environ_path = None
    if environ_path != None:
        return environ_path
    try:
        home_dir = expandOsPath("~/.config/regionanalysis/")
        return home_dir
    except:
        return os.path.join(module_dir, "database/")

def getPathDB(currunt_path):
    jsons = expandOsPath(os.path.join(currunt_path, "*.json"))
    files = [x for x in glob.glob(jsons)]
    genome_infos = map(loadJSON, files)
    return genome_infos

def getAllInstalledDB(module_dir):
    search_path = getAllPath(module_dir)
    installed_db = []
    for query_path in search_path:
        installed_db.extend(getPathDB(query_path))
    return installed_db

def getAnnoABPath(module_dir, genome, anno_db, RA_ver=None):
    installed_db = getAllInstalledDB(module_dir)
    cur_genome = None
    anno_db_ver = None
    for genome_info in installed_db:
        if genome_info["genome"] != genome:
            continue
        if (RA_ver is not None) and (genome_info["version"] != RA_ver):
            continue
        for db_info in genome_info["databases"]:
            if db_info["database"] == anno_db:
                if cur_genome is None or cur_genome["version"] < genome_info["version"]:
                    cur_genome = genome_info
                    anno_db_ver = db_info["version"]
    if cur_genome is not None:
        sys.stdout.write("%s database version %s of %s, RAver %s will be used.\n"%(anno_db,anno_db_ver, genome, cur_genome["version"]))
        return cur_genome
    return None

