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
from elasticsearch.exceptions import ConflictError

import bbconf
import pypiper
import yacman

parser = ArgumentParser(description="A pipeline to produce sets of bed files (bedsets) from bedbase")

parser.add_argument("-q", "--JSON-query-path", help="path to JSON file containing the query", type=str) # path to JSON with query
parser.add_argument("-b", "--bedbase-config", type=str, required=False, default=None,
                    help="path to the bedbase configuration file")
parser.add_argument("-n", "--bedset-name", help="name assigned to queried bedset", type=str)

# add pypiper args to make pipeline looper compatible
parser = pypiper.add_pypiper_args(parser, groups=["pypiper"], required=["--JSON-query-path", "--bedset_name"])


args = parser.parse_args()

# Initialize bbc object
bbc = bbconf.BedBaseConf(filepath=bbconf.get_bedbase_cfg(args.bedbase_config))

# SET OUTPUT FOLDER
# Create a folder to place pipeline logs 
logs_name = "bedbuncher_pipeline_logs"
logs_dir = os.path.abspath(os.path.join(
    bbc[CFG_PATH_KEY][CFG_BEDBUNCHER_OUTPUT_KEY], logs_name))

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

    :param Mapping sr: search results. Output of BedBaseConf.search_bedfiles
    :return str: bedset digest
    """
    from hashlib import md5
    m = md5()
    m.update(";".join(sorted([bf[JSON_MD5SUM_KEY][0] for bf in sr])).encode('utf-8'))
    return m.hexdigest()


def main():
    pm = pypiper.PipelineManager(name="bedbuncher", outfolder=logs_dir, args=args)

    # Establish Elasticsearch connection and check status using bbconf
    bbc.establish_elasticsearch_connection()

    # Use bbconf method to look for files in the ES index
    q = JSON_to_dict(args.JSON_query_path)
    search_results = bbc.search_bedfiles(query=q)
    nhits = len(search_results)
    if nhits < 1:
        raise BedBaseConfError("No BED files match the query: {}".format(q))
    print("{} BED files match the query".format(nhits))
    hit_ids = {i[JSON_ID_KEY][0]: i[JSON_MD5SUM_KEY][0] for i in search_results}
    bedset_digest = get_bedset_digest(search_results)
    print("bedset digest: {}".format(bedset_digest))

    # check for presence of the output folder and create it if needed
    output_folder = os.path.abspath(os.path.join(
        bbc[CFG_PATH_KEY][CFG_BEDBUNCHER_OUTPUT_KEY], bedset_digest))
    if not os.path.exists(output_folder):
        print("Output directory does not exist. Creating: {}".format(output_folder))
        os.makedirs(output_folder)

    # PRODUCE OUTPUT BEDSET PEP
    # Create PEP annotation and config files and TAR them along the queried .bed.gz files
    # produce basic csv annotation sheet based on IDs found in the bbc search
    print("Creating PEP annotation sheet and config.yaml for {}".
          format(args.bedset_name))
    pep_folder_path = os.path.join(output_folder, args.bedset_name + "_PEP")
    if not os.path.exists(pep_folder_path):
        os.makedirs(pep_folder_path)
    
    # define column names for annotation sheet based on LOLA requirements
    meta_list = [
            JSON_GENOME_KEY, JSON_PROTOCOL_KEY, JSON_CELL_TYPE_KEY, JSON_TISSUE_KEY, JSON_ANTIBODY_KEY, 
            JSON_TREATMENT_KEY, JSON_DATA_SOURCE_KEY, JSON_DESCRIPTION_KEY]
    bedset_pep_df = pd.DataFrame(columns=["sample_name", "output_file_path", "md5sum"] + meta_list)
    output_bed_path = "source1"
    for bedfiles in search_results:
        pep_metadata = {"sample_name": bedfiles[JSON_ID_KEY][0],
                        "output_file_path": output_bed_path,
                        "md5sum": bedfiles[JSON_MD5SUM_KEY][0],
                        "file_format": re.match('.*(.bed.*)$', bedfiles[BEDFILE_PATH_KEY][0]).group(1)}
        for key in meta_list:
            if key in bedfiles.keys():
                bed_file_meta = bedfiles[key][0]
                pep_metadata.update({key: bed_file_meta})
            else:
                pep_metadata.update({key:""})
        bedset_pep_df = bedset_pep_df.append(pep_metadata, ignore_index=True)

    bedset_annotation_sheet = args.bedset_name + '_annotation_sheet.csv'
    bedset_pep_path = os.path.join(pep_folder_path, bedset_annotation_sheet)
    bedset_pep_df.to_csv(bedset_pep_path, index=False)

    # Create df with bedfiles metadata: gc_content, num_regions, mean_abs_tss_dist, genomic_partitions
    bedstats_df = pd.DataFrame(columns=[JSON_MD5SUM_KEY, JSON_ID_KEY] + JSON_NUMERIC_KEY_VALUES)

    # Access elements in search object produced by bbc.search (both in metadata and statistics sections keys)
    print("Reading individual BED file statistics from Elasticsearch")
    for bed_file in search_results:
        bid = bed_file[JSON_ID_KEY][0]
        data = {JSON_MD5SUM_KEY: bed_file[JSON_MD5SUM_KEY][0], JSON_ID_KEY: bid}
        print("Processing: {}".format(bid))
        for key in JSON_NUMERIC_KEY_VALUES:
            try:
                bed_file_stat = bed_file[key][0]
            except KeyError:
                print("'{}' statistic not available for: {}".format(key, bid))
            else:
                data.update({key: bed_file_stat})
        bedstats_df = bedstats_df.append(data, ignore_index=True)
    bedstats_df = bedstats_df.dropna(1)
    # Calculate bedset statistics
    print("Calculating bedset statistics")
    bedfiles_means = bedstats_df.mean(axis=0)
    bedfiles_stdv = bedstats_df.std(axis=0)
    means_dictionary = dict(bedfiles_means)
    stdv_dictionary = dict(bedfiles_stdv)
    # Save bedstats_df matrix as csv file into the user defined output_folder
    bedfiles_stats_path = os.path.join(output_folder, args.bedset_name + '_bedstat.csv')
    print("Saving bedfiles statistics to: {}".format(bedfiles_stats_path))
    bedstats_df.to_csv(bedfiles_stats_path, index=False)

    # Save bedset_df  as csv file into the user defined output_folder
    means_df = pd.DataFrame(bedfiles_means, columns=["Mean"])
    stdv_df = pd.DataFrame(bedfiles_stdv, columns=["Standard Deviation"])
    bedset_df = pd.concat([means_df, stdv_df], axis=1)
    bedset_stats_path = os.path.join(output_folder, args.bedset_name + '_summaryStats.csv')
    print("Saving bedset statistics to: {}".format(bedset_stats_path))
    bedset_df.to_csv(bedset_stats_path)

    print("Creating iGD database")
    # IGD DATABASE
    # Need a .txt file with the paths to the queried bed files as input to the igd create command
    txt_bed_path = os.path.join(output_folder, args.bedset_name + '.txt')
    txt_file = open(txt_bed_path, "a")
    for files in search_results:
        bedfile_path = files[BEDFILE_PATH_KEY][0]
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
    igd_template = "igd create {bed_source_path} {igd_folder_path} {database_name} -f"
    gzip_template = "gzip {dir}"
    cmd = igd_template.format(bed_source_path=txt_bed_path, igd_folder_path=igd_folder_path, database_name=args.bedset_name)
    pm.run(cmd, target=os.path.join(igd_folder_path + ".tar.gz"))

    # TAR the iGD database folder
    igd_tar_archive_path = os.path.abspath(os.path.join(igd_folder_path + '.tar.gz'))
    with tarfile.open(igd_tar_archive_path, mode="w:gz", dereference=True, debug=3) as igd_tar:
        print("Creating iGD database TAR archive: {}".format(os.path.basename(igd_tar_archive_path)))
        igd_tar.add(igd_folder_path, arcname="", recursive=True, filter=flatten)

    # plot
    rscript_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tools",
        "bedsetStat.R")
    assert os.path.exists(rscript_path), \
        FileNotFoundError("'{}' script not found".format(rscript_path))
    json_file_path = os.path.join(output_folder, args.bedset_name + ".json")
    cmd_vars = dict(rscript=rscript_path,
                    bedfilelist=txt_bed_path,
                    output_folder=output_folder,
                    digest=bedset_digest,
                    id=args.bedset_name,
                    json=json_file_path)
    command = "Rscript {rscript} --outputfolder={output_folder} " \
              "--bedfilelist={bedfilelist} --id={id} --json={json}".\
        format(**cmd_vars)
    pm.run(cmd=command, target=json_file_path)


    # create yaml config file for newly produced bedset
    # create an empty file to write the cfg to
    cfg_path = os.path.join(pep_folder_path, args.bedset_name + "_cfg.yaml")
    open(cfg_path, 'a').close()
    config_attamp = yacman.YacAttMap(filepath=cfg_path)
    with config_attamp as y:
        y.pep_version = {}
        y.pep_version = "2.0.0"
        y.sample_table = {}
        y.sample_table = bedset_annotation_sheet
        y.looper = {}
        y.looper.output_dir = logs_dir 
        y.iGD_db = {}
        y.iGD_db = os.path.join(igd_folder_name, args.bedset_name + ".igd")
        y.iGD_index = {}
        y.iGD_index = os.path.join(igd_folder_name, args.bedset_name + "_index.tsv")
        y.sample_modifiers = {}
        y.sample_modifiers.append = {}
        y.sample_modifiers.append.output_file_path = "source1"
        y.sample_modifiers.derive = {}
        y.sample_modifiers.derive.attributes = ["output_file_path"] 
        y.sample_modifiers.derive.sources = {}
        y.sample_modifiers.derive.sources = {"source1": "{sample_name}" + "{file_format}"}

    # Create a tar archive using bed files original paths and bedset PEP
    tar_archive_file = os.path.abspath(os.path.join(output_folder, args.bedset_name + '.tar'))
    tar_archive = tarfile.open(tar_archive_file, mode="w:", dereference=True, debug=3)
    print("Creating TAR archive: {}".format(tar_archive_file))
    for files in search_results:
        bedfile_path = files[BEDFILE_PATH_KEY][0]
        tar_archive.add(bedfile_path, arcname=os.path.basename(bedfile_path), recursive=False, filter=None)
    tar_archive.add(pep_folder_path, arcname="", recursive=True, filter=flatten)
    tar_archive.close()

    # create a separate TAR.gz archive for the PEP annotation and config files
    pep_tar_archive_path = os.path.abspath(os.path.join(pep_folder_path + '.tar.gz'))
    with tarfile.open(pep_tar_archive_path, mode="w:gz", dereference=True, debug=3) as pep_tar:
        print("Creating PEP TAR archive: {}".format(os.path.basename(pep_tar_archive_path)))
        pep_tar.add(pep_folder_path, arcname="", recursive=True, filter=flatten)
    pm.clean_add(pep_folder_path)

    # read JSON produced in bedsetStat.R (with plot paths)
    with open(json_file_path, 'r', encoding='utf-8') as f:
        bedset_summary_info = json.loads(f.read())
    # update read data
    bedset_summary_info.update(
        {JSON_ID_KEY: args.bedset_name,
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
        bbc.insert_bedsets_data(data=bedset_summary_info, doc_id=bedset_digest)
    except ConflictError:
        from warnings import warn
        warn("Document '{}' exists in Elasticsearch. Not inserting.".
             format(bedset_digest))
    else:
        print("'{}' summary info was successfully inserted into {}".
              format(args.bedset_name, BEDSET_INDEX))
    pm.stop_pipeline()


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Pipeline aborted.")
        sys.exit(1)
