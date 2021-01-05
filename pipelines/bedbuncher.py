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
from bbconf import BedBaseConf
from hashlib import md5
from psycopg2.extras import Json

import bbconf
import pypiper
import yacman

parser = ArgumentParser(
    description="A pipeline to produce sets of bed files (bedsets)")

parser.add_argument("-q", "--query",
                    help="condition string to restrict bedfiles table select with",
                    type=str)
parser.add_argument("-v", "--query-val",
                    help="condition values to populate condition with",
                    type=str)
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
bbc = BedBaseConf(config_path=args.bedbase_config, database_only=True)

# Create a folder to place pipeline logs
logs_dir = os.path.abspath(
    os.path.join(bbc.get_bedbuncher_output_path(), "bedbuncher_pipeline_logs"))
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
    m.update(";".join(sorted([f["md5sum"] for f in sr])).encode('utf-8'))
    return m.hexdigest()


def mk_file_type(pth, title, md5sum):
    """

    :param pth:
    :param title:
    :return:
    """
    rel_to = os.path.join(bbc.get_bedbuncher_output_path(), md5sum)
    return {"path": os.path.relpath(pth, rel_to), "title": title}


def main():
    pm = pypiper.PipelineManager(name="bedbuncher", outfolder=logs_dir, args=args)
    # Use bbconf method to look for files in the database
    search_results = bbc.bed.select(
        condition=args.query, condition_val=[args.query_val])
    nhits = len(search_results)
    if nhits < 2:
        raise BedBaseConfError(f"{nhits} BED files match the query: {args.query}")
    pm.info(f"{nhits} BED files match the query")
    hit_ids = [i["id"] for i in search_results]
    bedset_digest = get_bedset_digest(search_results)
    pm.info(f"bedset digest: {bedset_digest}")

    # check for presence of the output folder and create it if needed
    output_folder = os.path.abspath(
        os.path.join(bbc.get_bedbuncher_output_path(), bedset_digest))
    if not os.path.exists(output_folder):
        pm.info(f"Output directory does not exist. Creating: {output_folder}")
        os.makedirs(output_folder)

    # Create TrackHub directory
    pm.info(f"Creating TrackHub directory for {args.bedset_name}")
    hub_folder = os.path.join(output_folder, 'bedsetHub')
    if not os.path.exists(hub_folder):
        os.makedirs(hub_folder)
    # write hub.txt file
    hub_txt = {'hub': 'BEDBASE_'+args.bedset_name,
            'shortLabel': 'BEDBASE_'+args.bedset_name,
            'longLabel': args.bedset_name + ' signal tracks',
            'genomesFile': 'http://dev1.bedbase.org/api/bedset/' + bedset_digest + '/file/hub_file',
            'email': 'bx2ur@virginia.edu',
            'descriptionUrl': 'http://dev1.bedbase.org/'}
    f = open(os.path.join(hub_folder, "hub.txt"),"w")
    f.writelines('{}\t{}\n'.format(k,v) for k, v in hub_txt.items()) 
    f.close()
    # write genomes.txt and trackDb.txt
    genomes = {}
    for bedfiles in search_results:
        if bedfiles['other']["genome"] not in genomes.keys():
            genomes_txt = {'genome': bedfiles['other']["genome"],
                'trackDb': 'http://dev1.bedbase.org/api/bedset/' + bedset_digest + '/file/genomes_file?genome='+ bedfiles['other']["genome"]}
            if os.path.exists(os.path.join(hub_folder, "genomes.txt")):
                f = open(os.path.join(hub_folder, "genomes.txt"),"a")
                f.writelines('{}\t{}\n'.format(k,v) for k, v in genomes_txt.items())
                f.close()
            else:
                f = open(os.path.join(hub_folder, "genomes.txt"),"w")
                f.writelines('{}\t{}\n'.format(k,v) for k, v in genomes_txt.items())
                f.close()

            genome_folder = os.path.join(hub_folder,bedfiles['other']["genome"])
            if not os.path.exists(genome_folder):
                os.makedirs(genome_folder)
            
            trackDb_txt = {'track': bedfiles["name"],
                            'type': 'bigBed',
                            'bigDataUrl': 'http://dev1.bedbase.org/api/bed/' + bedfiles["md5sum"] + '/file/bigbed_file',
                            'shortLabel': bedfiles["name"],
                            'longLabel': bedfiles["other"]["description"]}
            if os.path.exists(os.path.join(genome_folder, "trackDb.txt")):
                f = open(os.path.join(genome_folder, "trackDb.txt"),"a")
                f.writelines('{}\t{}\n'.format(k,v) for k, v in trackDb_txt.items())
                f.close()
            else:
                f = open(os.path.join(genome_folder, "trackDb.txt"),"w")
                f.writelines('{}\t{}\n'.format(k,v) for k, v in trackDb_txt.items())
                f.close()
            genomes.update({bedfiles['other']["genome"] : os.path.join(genome_folder, "trackDb.txt")})

    print ("test print:\n",genomes)
    for genome, path in genomes.items():
        genomes[genome] = mk_file_type(
             path, genome+" trackDb.txt files for the BED set", bedset_digest) 

    # PRODUCE OUTPUT BEDSET PEP
    # Create PEP annotation and config files and TAR them along the queried
    # .bed.gz files produce basic csv annotation sheet based on IDs found in
    # the bbc search
    pm.info(f"Creating PEP annotation sheet and config.yaml for {args.bedset_name}")
    pep_folder_path = os.path.join(output_folder, args.bedset_name + "_PEP")
    if not os.path.exists(pep_folder_path):
        os.makedirs(pep_folder_path)
    
    # define column names for annotation sheet based on LOLA requirements
    meta_list = [
        "genome", "protocol", "cell_type",
        "tissue", "antibody", "treatment",
        "data_source", "description"
    ]
    bedset_pep_df = pd.DataFrame(
        columns=["sample_name", "output_file_path", "md5sum"] + meta_list)
    output_bed_path = "source1"
    for bedfiles in search_results:
        file_fmt = re.match('.*(.bed.*)$', bedfiles["bedfile"]["path"]).group(1)
        pep_metadata = {"sample_name": bedfiles["name"],
                        "output_file_path": output_bed_path,
                        "md5sum": bedfiles["md5sum"],
                        "file_format": file_fmt}
        for key in meta_list:
            if key in bedfiles['other'].keys():
                bed_file_meta = bedfiles['other'][key]
                pep_metadata.update({key: bed_file_meta})
            else:
                pep_metadata.update({key: ""})
        bedset_pep_df = bedset_pep_df.append(pep_metadata, ignore_index=True)

    bedset_annotation_sheet = args.bedset_name + '_annotation_sheet.csv'
    bedset_pep_path = os.path.join(pep_folder_path, bedset_annotation_sheet)
    bedset_pep_df.to_csv(bedset_pep_path, index=False)

    numeric_results = [k for k, v in bbc.bed.schema.items()
                       if v["type"] == "number"]
    # Create df with bedfiles metadata
    bedstats_df = pd.DataFrame(columns=["md5sum", "name"] + numeric_results)

    # Access elements in search object produced by bbc.search (both in metadata and statistics sections keys)
    pm.info("Reading individual BED file statistics from the database")
    for bed_file in search_results:
        bid = bed_file["name"]
        data = {"md5sum": bed_file["md5sum"], "name": bid}
        pm.info(f"Processing: {bid}")
        for key in numeric_results:
            try:
                bed_file_stat = bed_file[key]
            except KeyError:
                pm.info(f"'{key}' statistic not available for: {bid}")
            else:
                data.update({key: bed_file_stat})
        bedstats_df = bedstats_df.append(data, ignore_index=True)
    bedstats_df = bedstats_df.dropna(1)
    bedstats_df[numeric_results] = bedstats_df[numeric_results].apply(pd.to_numeric)
    # Calculate bedset statistics
    pm.info("Calculating bedset statistics")
    bedfiles_means = bedstats_df.mean(axis=0)
    bedfiles_stdv = bedstats_df.std(axis=0)
    print(f"bedstats_df: {bedstats_df}")
    print(f"bedfiles_stdv: {bedfiles_stdv}")
    means_dictionary = dict(bedfiles_means)
    stdv_dictionary = dict(bedfiles_stdv)
    print(f"stdv_dictionary: {stdv_dictionary}")
    # Save bedstats_df matrix as csv file into the user-defined output_folder
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
        bedfile_path = files["bedfile"]["path"]
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
    pm.run(cmd, target=os.path.join(igd_folder_path + ".tar.gz"), nofail=True)

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
        y.sample_modifiers.derive.sources = {"source1": "'{sample_name}{file_format}'"}

    # Create a tar archive using bed files original paths and bedset PEP
    tar_archive_file = os.path.abspath(
        os.path.join(output_folder, args.bedset_name + '.tar'))
    tar_archive = tarfile.open(
        tar_archive_file, mode="w:", dereference=True, debug=3)
    pm.info(f"Creating TAR archive: {tar_archive_file}")
    for files in search_results:
        bedfile_path = files["bedfile"]["path"]
        if not os.path.isabs(bedfile_path):
            bedfile_path = os.path.realpath(os.path.join(bbc.get_bedbuncher_output_path(), files["md5sum"], bedfile_path))
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

    for plot in bedset_summary_info["plots"]:
        plot_id = plot["name"]
        del plot["name"]
        bedset_summary_info.update({plot_id: plot})
    del bedset_summary_info["plots"]

    # TODO: source the key names from bbconf package?
    bedset_summary_info.update(
        {"name": args.bedset_name,
         "bedset_means": means_dictionary,
         "bedset_standard_deviation": stdv_dictionary,
         "bedset_tar_archive_path": mk_file_type(
             tar_archive_file, "TAR archive with BED files in this BED set", bedset_digest),
         "bedset_bedfiles_gd_stats": mk_file_type(
             bedfiles_stats_path, "Statistics of the BED files in this BED set", bedset_digest),
         "bedset_gd_stats": mk_file_type(
             bedset_stats_path, "Means and standard deviations of the BED files in this BED set", bedset_digest),
         "bedset_igd_database_path": mk_file_type(
             igd_tar_archive_path, "iGD database", bedset_digest),
         "bedset_pep": mk_file_type(
             pep_tar_archive_path, "PEP including BED files in this BED set", bedset_digest),
         "md5sum": bedset_digest, 
         "bedset_trackdb_path": genomes}, 
         )

    # select only first element of every list due to JSON produced by R putting
    # every value into a list
    data = {k.lower(): v[0] if (isinstance(v, list))
            else v for k, v in bedset_summary_info.items()}
    bedset_id = bbc.bedset.report(
        record_identifier=bedset_digest, values=data, return_id=True)
    for hit_id in hit_ids:
        bbc.report_relationship(bedset_id=bedset_id, bedfile_id=hit_id)
    pm.stop_pipeline()


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Pipeline aborted.")
        sys.exit(1)
