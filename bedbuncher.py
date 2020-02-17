#!/usr/bin/env python3
"""
bedset paths and info generating pipeline
"""

__author__ = ["Jose Verdezoto", "Michal Stolarczyk"]
__email__ = "jev4xy@virginia.edu"
__version__ = "0.0.1"

from argparse import ArgumentParser
import pypiper
import os
import sys
import pandas as pd
import bbconf
from bbconf.const import *
from bbconf.exceptions import BedBaseConfError
import json
import tarfile

parser = ArgumentParser(description="A pipeline to produce sets of bed files (bedsets) from bedbase")

parser.add_argument("-q", "--JSON-query-path", help="path to JSON file containing the query", type=str) # path to JSON with query
parser.add_argument("-c", "--bedbase-config", type=str, required=False, default=None,
                    help="path to the bedbase configuration file")
parser.add_argument("-b", "--bedset-name", help="name assigned to queried bedset", default=str )
parser.add_argument("-f", "--output-folder", help="path to folder where tar file and igd database will be stored", default=str )

# add pypiper args to make pipeline looper compatible
parser = pypiper.add_pypiper_args(parser, groups=["pypiper", "looper"],
                                  required=["--JSON-query-path", "--bedset_name", "--output-folder"])

args = parser.parse_args()

# SET OUTPUT FOLDER
# use output parent argument from looper to place pipeline stats (live separately from bedset results)
out_parent = args.output_parent

bbc = bbconf.BedBaseConf(filepath=bbconf.get_bedbase_cfg())


def JSON_to_dict(file_name):
    with open(file_name) as f_in:
        return json.load(f_in)


def main():
    pm = pypiper.PipelineManager(name="bedbuncher", outfolder=out_parent, args=args)

    # Establish Elasticsearch connection and check status using bbconf
    bbc.establish_elasticsearch_connection()

    # Use bbconf method to look for files in the ES index
    q = JSON_to_dict(args.JSON_query_path)
    search_results = bbc.search_bedfiles(query=q)
    nhits = len(search_results)
    if nhits < 1:
        raise BedBaseConfError("No BED files match the query: {}".format(q))
    print("{} BED files match the query".format(nhits))

    # check for presence of the output folder and create it if needed
    output_folder = os.path.dirname(args.output_folder)
    if not os.path.exists(output_folder):
        print("Output directory does not exist. Creating: {}".format(output_folder))
        os.makedirs(output_folder)

    # Create a tar archive using the paths to the bed files provided by the bbconf search object
    tar_archive_file = os.path.abspath(os.path.join(args.output_folder, args.bedset_name + '.tar'))
    tar_archive = tarfile.open(tar_archive_file, mode="w:", dereference=True, debug=3)
    print("Creating TAR archive: {}".format(tar_archive_file))
    for files in search_results:
        bedfile_path = files[BEDFILE_PATH_KEY][0]
        tar_archive.add(bedfile_path, arcname=os.path.basename(bedfile_path), recursive=False, filter=None)
    tar_archive.close()

    # Create df with bedfiles metadata: gc_content, num_regions, mean_abs_tss_dist, genomic_partitions
    bedstats_df = pd.DataFrame(columns=[JSON_ID_KEY] + JSON_NUMERIC_KEY_VALUES)

    # transform individual stats from dictionary into floats to perform calculations, if needed
    def make_float(es_element):
        float(es_element[0])

    # Access elements in search object produced by bbc.search
    print("Reading individual BED file statistics from Elasticsearch")
    for bed_file in search_results:
        bid = bed_file[JSON_ID_KEY][0]
        data = {JSON_ID_KEY: bid}
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
    bedfiles_stats_path = os.path.join(args.output_folder, args.bedset_name + '_bedstat.csv')
    print("Saving bedfiles statistics to: {}".format(bedfiles_stats_path))
    bedstats_df.to_csv(bedfiles_stats_path, index=False)

    # Save bedset_df  as csv file into the user defined output_folder
    means_df = pd.DataFrame(bedfiles_means, columns=["Mean"])
    stdv_df = pd.DataFrame(bedfiles_stdv, columns=["Standard Deviation"])
    bedset_df = pd.concat([means_df, stdv_df], axis=1)
    bedset_stats_path = os.path.join(args.output_folder, args.bedset_name + '_summaryStats.csv')
    print("Saving bedset statistics to: {}".format(bedset_stats_path))
    bedset_df.to_csv(bedset_stats_path)

    
    print("Creating iGD database")
    # IGD DATABASE
    # Need a .txt file with the paths to the queried bed files as input to the igd create command
    txt_bed_path = os.path.join(args.output_folder, args.bedset_name + '.txt')
    txt_file = open(txt_bed_path, "a")
    for files in search_results:
        bedfile_path = files[BEDFILE_PATH_KEY][0]
        bedfile_target = os.readlink(bedfile_path) \
            if os.path.islink(bedfile_path) else bedfile_path
        txt_file.write("{}\r\n".format(bedfile_target))
    txt_file.close()
    pm.clean_add(txt_bed_path)
    # iGD database
    igd_folder_name = args.bedset_name + "_igd"
    igd_folder_path = os.path.join(args.output_folder, igd_folder_name)
    if not os.path.exists(igd_folder_path):
        os.makedirs(igd_folder_path)

    # Command templates for IGD database construction
    igd_template = "igd create {bed_source_path} {igd_folder_path} {database_name} -f"
    gzip_template = "gzip -r {dir}"
    cmd1 = igd_template.format(bed_source_path=txt_bed_path, igd_folder_path=igd_folder_path, database_name=args.bedset_name)
    cmd2 = gzip_template.format(dir=igd_folder_path)
    cmd = [cmd1, cmd2]
    pm.run(cmd, target=os.path.join(igd_folder_path, args.bedset_name + ".igd"))

    # create a nested dictionary with avgs,stdv, paths to tar archives, bedset csv file and igd database.
    bedset_summary_info = {JSON_ID_KEY: args.bedset_name,
                           JSON_BEDSET_MEANS_KEY: means_dictionary,
                           JSON_BEDSET_SD_KEY: stdv_dictionary,
                           JSON_BEDSET_TAR_PATH_KEY: [tar_archive_file],
                           JSON_BEDSET_BEDFILES_GD_STATS_KEY: [bedfiles_stats_path],
                           JSON_BEDSET_GD_STATS: [bedset_stats_path],
                           JSON_BEDSET_IGD_DB_KEY: [igd_folder_path]}

    # Insert bedset information into BEDSET_INDEX
    bbc.insert_bedsets_data(data=bedset_summary_info)
    print("'{}' summary info was successfully inserted into {}".format(args.bedset_name, BEDSET_INDEX))
    pm.stop_pipeline()


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Pipeline aborted.")
        sys.exit(1)
