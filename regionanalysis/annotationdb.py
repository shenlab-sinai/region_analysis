import os
import sys
import json

def expandOsPath(path):
    """
    To expand the path with shell variables.
    Arguments:
    - `path`: path string
    """
    return os.path.expanduser(os.path.expandvars(path))

def getAnnoABPath(module_dir, genome, anno_db):
    search_path = [os.path.join(module_dir, "database/")]
    environ_path = expandOsPath(os.environ.get("RA_DB_PATH"))
    if environ_path != None:
        search_path.append(environ_path)
    home_dir = expandOsPath("~/.config/regionanalysis")
    if os.path.isdir(home_dir):
        search_path.append(home_dir)
    try:
        for query_dir in search_path:
            json_file = os.path.join(query_dir, genome+".json")
            if os.path.isfile(json_file):
                fp = open(json_file)
                genome_info = json.load(fp)
                fp.close()
                for db_info in genome_info["databases"]:
                    if db_info["database"] == anno_db:
                        sys.stderr.write("%s database version %s of %s will be used.\n"%(anno_db, db_info["version"], genome))
                        return (query_dir, genome_info)
        raise SystemExit
    except SystemExit:
        sys.stderr.write("%s not in the genome database!\n"%genome)
        exit()