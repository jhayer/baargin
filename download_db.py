#!/usr/bin/env python3

import requests
from git import Repo
import shutil
import re,os,sys
import logging
import argparse
import yaml
import tarfile
from yaml import load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

##########################
#        MAIN            #
##########################

if __name__ == '__main__':

    parser = argparse.ArgumentParser( description = __doc__ )

    #single character parameter
    parser.add_argument("-y", "--yaml", default="db.yaml", help="YAML file containing the databases to download.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing databases if any.")
    parser.add_argument("-o", "--output", default="db", help="Path where to save the downloaded databases.")
    parser.add_argument("-q", "--quiet", action="count", default=0, help="Decrease verbosity.")
    parser.add_argument("-v", "--verbose", action="count", default=2, help="Increase verbosity.")

    args = parser.parse_args()

    logging.basicConfig(format='%(asctime)s %(levelname)s %(module)s: %(message)s',
                        level = (4-args.verbose+args.quiet)*10,
                        datefmt="%H:%M:%S")

    # test output exists
    output_file = args.output
    if output_file and os.path.exists(output_file):
        if not args.overwrite :
            logging.info(output_file + " output already exists. Please use --overwrite parameter to skip this warning.")
            sys.exit()
        else:
            shutil.rmtree(output_file)
            os.mkdir(output_file)
    else:
        os.mkdir(output_file)

    # check yaml existence:
    yaml_file = args.yaml
    if yaml_file and os.path.exists(yaml_file) :
        with open(yaml_file, 'r') as file :
            yamlDict = yaml.safe_load(file)

            # loop over databases
            for db in yamlDict :
                logging.info("Downloading " + db + "database...")
                os.mkdir(args.output+"/"+db)

                # get URL
                url = yamlDict[db]
                db_basename = os.path.basename(url) # get file name
                path = args.output+"/"+db+"/"+db_basename


                # bitbucket case, use curl
                if "bitbucket.org" in url or "github.com" in url :
                    Repo.clone_from(url, path)
                    #cmd = "curl "+url+" --output "+path
                    #os.system(cmd)
                # use request approach
                else :
                    r = requests.get(url, allow_redirects=True)
                    # Save the content
                    open(path, 'wb').write(r.content)

                # check extension to uncompress if needed
                # bz2 extension
                if db_basename.endswith('.bz2'):
                    tar = tarfile.open(path, "r:bz2")
                    tar.extractall(args.output+"/"+db+"/")
                    tar.close()
    else :
        logging.error(str(yaml_file) + " input YAML is not a file (--yaml).\n")
        sys.exit()
