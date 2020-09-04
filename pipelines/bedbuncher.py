#!/usr/bin/env python3
"""
bedset paths, statistics and PEP generating pipeline
"""

__author__ = ["Jose Verdezoto", "Michal Stolarczyk"]
__email__ = "jev4xy@virginia.edu"
__version__ = "0.0.2-dev"

import re
import os
import sys
import pandas as pd
import json
import tarfile

from argparse import ArgumentParser
from bbconf.const import *
from bbconf.exceptions import BedBaseConfError
from hashlib import md5

import bbconf
import pypiper
import yacman

parser = ArgumentParser(
    description="A pipeline to produce sets of bed files (bedsets)")

parser.add_argument("-q", "--query",
                    help="condition to restrict bedfiles table select with",
                    type=str)
parser.add_argument("-j", "--json",
                    help="the JSONB database column the provided query "
                         "should be applied to",
                    type=str, required=False, default=None)
parser.add_argument("-b", "--bedbase-config",
                    type=str, required=False, default=None,
                    help=f"path to the bedbase configuration file. "
                         f"Optional if {CFG_ENV_VARS} env vars are set.")
parser.add_argument("-n", "--bedset-name",
                    help="name assigned to the bedset to be created", type=str)

# add pypiper args to make pipeline looper compatible
parser = pypiper.add_pypiper_args(parser,
                                  groups=["pypiper"],
                                  required=["--query", "--bedset_name"])
args = parser.parse_args()

# Initialize bbc object
bbc = bbconf.BedBaseConf(filepath=bbconf.get_bedbase_cfg(args.bedbase_config))

# Create a folder to place pipeline logs
logs_dir = os.path.abspath(os.path.join(
    bbc[CFG_PATH_KEY][CFG_BEDBUNCHER_OUTPUT_KEY], "bedbuncher_pipeline_logs"))
if not os.path.exists(logs_dir):
    print("bedbuncher pipeline logs directory doesn't exist. Creating one...")
    os.makedirs(logs_dir)


def JSON_to_dict(file_name):
    """
    read json file to a dict

    :param str file_name: path to the file to read
    :return Mapping: read file
    """
    with open(file_name) as f_in:
        return json.load(f_in)


def flatten(tarinfo):
    """
    function to tar only the contents of a folder, excluding the dir hierarchy

    :param tarinfo:
    :return:
    """
    tarinfo.name = os.path.basename(tarinfo.name)
    return tarinfo


def get_bedset_digest(sr):
    """
    Get a unique id for a bedset, specific to its contents.

    :param list[psycopg2.extras.DictRow] sr: search results; output of
        BedBaseConf.select
    :return str: bedset digest
    """
    m = md5()
    m.update(";".join(sorted([f[JSON_MD5SUM_KEY] for f in sr])).encode('utf-8'))
    return m.hexdigest()


def main():
    pm = pypiper.PipelineManager(name="bedbuncher", outfolder=logs_dir, args=args)
    # Use bbconf method to look for files in the database
    search_results = bbc.select(
        condition=args.query, json=args.json, table_name=BED_TABLE)
    nhits = len(search_results)
    if nhits < 1:
        raise BedBaseConfError(f"No BED files match the query: {args.query}")
    pm.info(f"{nhits} BED files match the query")
    hit_ids = {i[JSON_NAME_KEY]: i[JSON_MD5SUM_KEY] for i in search_results}
    bedset_digest = get_bedset_digest(search_results)
    pm.info(f"bedset digest: {bedset_digest}")

    # check for presence of the output folder and create it if needed
    output_folder = os.path.abspath(os.path.join(
        bbc[CFG_PATH_KEY][CFG_BEDBUNCHER_OUTPUT_KEY], bedset_digest))
    if not os.path.exists(output_folder):
        pm.info("Output directory does not exist. Creating: {}".format(output_folder))
        os.makedirs(output_folder)

    # PRODUCE OUTPUT BEDSET PEP
    # Create PEP annotation and config files and TAR them along the queried
    # .bed.gz files produce basic csv annotation sheet based on IDs found in
    # the bbc search
    pm.info("Creating PEP annotation sheet and config.yaml for {}".
          format(args.bedset_name))
    pep_folder_path = os.path.join(output_folder, args.bedset_name + "_PEP")
    if not os.path.exists(pep_folder_path):
        os.makedirs(pep_folder_path)
    
    # define column names for annotation sheet based on LOLA requirements
    meta_list = [
        JSON_GENOME_KEY, JSON_PROTOCOL_KEY, JSON_CELL_TYPE_KEY,
        JSON_TISSUE_KEY, JSON_ANTIBODY_KEY, JSON_TREATMENT_KEY,
        JSON_DATA_SOURCE_KEY, JSON_DESCRIPTION_KEY
    ]
    bedset_pep_df = pd.DataFrame(
        columns=["sample_name", "output_file_path", "md5sum"] + meta_list)
    output_bed_path = "source1"
    for bedfiles in search_results:
        file_fmt = re.match('.*(.bed.*)$', bedfiles[BEDFILE_PATH_KEY]).group(1)
        pep_metadata = {"sample_name": bedfiles[JSON_NAME_KEY],
                        "output_file_path": output_bed_path,
                        "md5sum": bedfiles[JSON_MD5SUM_KEY],
                        "file_format": file_fmt}
        for key in meta_list:
            if key in bedfiles.keys():
                bed_file_meta = bedfiles[key]
                pep_metadata.update({key: bed_file_meta})
            else:
                pep_metadata.update({key: ""})
        bedset_pep_df = bedset_pep_df.append(pep_metadata, ignore_index=True)

    bedset_annotation_sheet = args.bedset_name + '_annotation_sheet.csv'
    bedset_pep_path = os.path.join(pep_folder_path, bedset_annotation_sheet)
    bedset_pep_df.to_csv(bedset_pep_path, index=False)

    # Create df with bedfiles metadata
    bedstats_df = pd.DataFrame(
        columns=[JSON_MD5SUM_KEY, JSON_NAME_KEY] + JSON_FLOAT_KEY_VALUES)

    # Access elements in search object produced by bbc.search (both in metadata and statistics sections keys)
    pm.info("Reading individual BED file statistics from the database")
    for bed_file in search_results:
        bid = bed_file[JSON_NAME_KEY]
        data = {JSON_MD5SUM_KEY: bed_file[JSON_MD5SUM_KEY], JSON_NAME_KEY: bid}
        pm.info(f"Processing: {bid}")
        for key in JSON_FLOAT_KEY_VALUES:
            try:
                bed_file_stat = bed_file[key]
            except KeyError:
                pm.info(f"'{key}' statistic not available for: {bid}")
            else:
                data.update({key: bed_file_stat})
        bedstats_df = bedstats_df.append(data, ignore_index=True)
    bedstats_df = bedstats_df.dropna(1)
    # Calculate bedset statistics
    pm.info("Calculating bedset statistics")
    bedfiles_means = bedstats_df.mean(axis=0)
    bedfiles_stdv = bedstats_df.std(axis=0)
    means_dictionary = dict(bedfiles_means)
    stdv_dictionary = dict(bedfiles_stdv)
    # Save bedstats_df matrix as csv file into the user defined output_folder
    bedfiles_stats_path = os.path.join(output_folder, args.bedset_name + '_bedstat.csv')
    pm.info(f"Saving bedfiles statistics to: {bedfiles_stats_path}")
    bedstats_df.to_csv(bedfiles_stats_path, index=False)

    # Save bedset_df  as csv file into the user defined output_folder
    means_df = pd.DataFrame(bedfiles_means, columns=["Mean"])
    stdv_df = pd.DataFrame(bedfiles_stdv, columns=["Standard Deviation"])
    bedset_df = pd.concat([means_df, stdv_df], axis=1)
    bedset_stats_path = os.path.join(output_folder, args.bedset_name + '_summaryStats.csv')
    pm.info(f"Saving bedset statistics to: {bedset_stats_path}")
    bedset_df.to_csv(bedset_stats_path)

    pm.info("Creating iGD database")
    # IGD DATABASE
    # Need a .txt file with the paths to the queried bed files as input to the igd create command
    txt_bed_path = os.path.join(output_folder, args.bedset_name + '.txt')
    txt_file = open(txt_bed_path, "a")
    for files in search_results:
        bedfile_path = files[BEDFILE_PATH_KEY]
        bedfile_target = os.readlink(bedfile_path) \
            if os.path.islink(bedfile_path) else bedfile_path
        txt_file.write("{}\n".format(bedfile_target))
    txt_file.close()
    pm.clean_add(txt_bed_path)
    
    # iGD database
    igd_folder_name = args.bedset_name + "_igd"
    igd_folder_path = os.path.join(output_folder, igd_folder_name)
    if not os.path.exists(igd_folder_path):
        os.makedirs(igd_folder_path)
    pm.clean_add(igd_folder_path)

    # Command templates for IGD database construction
    cmd = f"igd create {txt_bed_path} {igd_folder_path} {args.bedset_name} -f"
    pm.run(cmd, target=os.path.join(igd_folder_path + ".tar.gz"))

    # TAR the iGD database folder
    igd_tar_archive_path = os.path.abspath(os.path.join(igd_folder_path + '.tar.gz'))
    with tarfile.open(igd_tar_archive_path, mode="w:gz", dereference=True, debug=3) as igd_tar:
        pm.info(f"Creating iGD database TAR archive: {igd_tar_archive_path}")
        igd_tar.add(igd_folder_path, arcname="", recursive=True, filter=flatten)

    # plot
    rscript_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tools",
        "bedsetStat.R")
    assert os.path.exists(rscript_path), \
        FileNotFoundError(f"'{rscript_path}' script not found")
    json_file_path = os.path.join(output_folder, args.bedset_name + ".json")
    command = f"Rscript {rscript_path} --outputfolder={output_folder} " \
              f"--bedfilelist={txt_bed_path} --id={args.bedset_name} " \
              f"--json={json_file_path}"
    pm.run(cmd=command, target=json_file_path)

    # create yaml config file for newly produced bedset
    # create an empty file to write the cfg to
    cfg_path = os.path.join(pep_folder_path, args.bedset_name + "_cfg.yaml")
    open(cfg_path, 'a').close()
    config_attamp = yacman.YacAttMap(filepath=cfg_path)
    with config_attamp as y:
        y.pep_version = "2.0.0"
        y.sample_table = bedset_annotation_sheet
        y.looper = {}
        y.looper.output_dir = logs_dir 
        y.iGD_db = os.path.join(igd_folder_name, args.bedset_name + ".igd")
        y.iGD_index = os.path.join(igd_folder_name, args.bedset_name + "_index.tsv")
        y.sample_modifiers = {}
        y.sample_modifiers.append = {}
        y.sample_modifiers.append.output_file_path = "source1"
        y.sample_modifiers.derive = {}
        y.sample_modifiers.derive.attributes = ["output_file_path"] 
        y.sample_modifiers.derive.sources = {}
        y.sample_modifiers.derive.sources = {"source1": "{sample_name}{file_format}"}

    # Create a tar archive using bed files original paths and bedset PEP
    tar_archive_file = os.path.abspath(
        os.path.join(output_folder, args.bedset_name + '.tar'))
    tar_archive = tarfile.open(
        tar_archive_file, mode="w:", dereference=True, debug=3)
    pm.info(f"Creating TAR archive: {tar_archive_file}")
    for files in search_results:
        bedfile_path = files[BEDFILE_PATH_KEY]
        tar_archive.add(bedfile_path, arcname=os.path.basename(bedfile_path),
                        recursive=False, filter=None)
    tar_archive.add(pep_folder_path, arcname="", recursive=True, filter=flatten)
    tar_archive.close()

    # create a separate TAR.gz archive for the PEP annotation and config files
    pep_tar_archive_path = os.path.abspath(os.path.join(pep_folder_path + '.tar.gz'))
    with tarfile.open(pep_tar_archive_path, mode="w:gz", dereference=True, debug=3) \
            as pep_tar:
        pm.info(f"Creating PEP TAR archive: {pep_tar_archive_path}")
        pep_tar.add(pep_folder_path, arcname="", recursive=True, filter=flatten)
    pm.clean_add(pep_folder_path)

    # read JSON produced in bedsetStat.R (with plot paths)
    with open(json_file_path, 'r', encoding='utf-8') as f:
        bedset_summary_info = json.loads(f.read())
    # update read data
    bedset_summary_info.update(
        {JSON_NAME_KEY: args.bedset_name,
         JSON_BEDSET_MEANS_KEY: means_dictionary,
         JSON_BEDSET_SD_KEY: stdv_dictionary,
         JSON_BEDSET_TAR_PATH_KEY: [tar_archive_file],
         JSON_BEDSET_BEDFILES_GD_STATS_KEY: [bedfiles_stats_path],
         JSON_BEDSET_GD_STATS_KEY: [bedset_stats_path],
         JSON_BEDSET_IGD_DB_KEY: [igd_tar_archive_path],
         JSON_BEDSET_PEP_KEY: [pep_tar_archive_path],
         JSON_BEDSET_BED_IDS_KEY: [hit_ids],
         JSON_MD5SUM_KEY: [bedset_digest]})

    # Insert bedset information into BEDSET_INDEX
    try:
        bbc.insert_bedset_data(values=bedset_summary_info)
    pm.info(f"'{args.bedset_name}' summary info was successfully inserted into {BEDSET_TABLE}")
    pm.stop_pipeline()


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Pipeline aborted.")
        sys.exit(1)
