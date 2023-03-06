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


def downloadFileWithProgress(url, local_file, sha1=None, chunk_size=8192):
  import requests, sys
  """
  Function to download a file from the specified URL.
  Outputs a basic progress bar and download stats to the CLI.
  Optionally downloads to a specified download folder, and checks the SHA1 checksum of the file.
  Chunk size can also be specified, to control the download.
  Uses 'Requests' :: http://www.python-requests.org
  """
  def size_fmt(numBytes):
    for symbol in ['B','KB','MB','GB','TB','EB','ZB']:
      if numBytes < 1024.0:
        return "{0:3.1f} {1}".format(numBytes, symbol)
      else:
        numBytes /= 1024.0
    # Return Yottabytes if all else fails.
    return "{0:3.1f} {1}".format(numBytes, 'YB')

  r = requests.get(url, stream=True)
  r.raise_for_status()

  file_size = int(r.headers['Content-Length'])
  dl_size = 0

  print("Downloading: {0}; {1}".format(local_file.split('/')[-1], size_fmt(file_size)))
  with open(local_file, 'wb') as f:
    for chunk in r.iter_content(chunk_size):
      if chunk: # filter out keep-alive new chunks
        dl_size += len(chunk)
        f.write(chunk)
        f.flush()
        percentage = (dl_size * 100. / file_size)
        num_equals = int(round(percentage/4))

        sys.stdout.write("[{0:25}] {1:>10} [{2: 3.2f}%]\r".format('='*num_equals, size_fmt(dl_size), percentage))
        sys.stdout.flush()
  print
  if sha1 is not None:
    import hashlib
    print("Verifying download...")
    f = open(local_file, 'rb')
    download_checksum = hashlib.sha1(f.read()).hexdigest()
    f.close()
    if download_checksum == sha1:
      return local_file
    else:
      raise Exception("Checksum mismatch")
  else:
    return local_file

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
                logging.info("Downloading " + db + " database...")
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
                    downloadFileWithProgress(url, path)

                # check extension to uncompress if needed
                # bz2 extension
                if db_basename.endswith('.bz2'):
                    tar = tarfile.open(path, "r:bz2")
                    tar.extractall(args.output+"/"+db+"/")
                    tar.close()
                # gz extension
                elif db_basename.endswith('.gz') or db_basename.endswith('.tgz'):
                    tar = tarfile.open(path, "r:gz")
                    tar.extractall(args.output+"/"+db+"/")
                    tar.close()
                elif db_basename.endswith(".tar"):
                    tar = tarfile.open(fname, "r:")
                    tar.extractall(args.output+"/"+db+"/")
                    tar.close()
    else :
        logging.error(str(yaml_file) + " input YAML is not a file (--yaml).\n")
        sys.exit()
