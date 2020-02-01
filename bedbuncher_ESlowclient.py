#!/usr/bin/env python
"""
bedset paths and info generating pipeline

"""

__author__ = ["Jose Verdezoto"]
__email__ = "jev4xy@virginia.edu"
__version__ = "0.0.1"

from argparse import ArgumentParser
import pypiper, os, sys
import pandas as pd
from elasticsearch import Elasticsearch
import bbconf
from bbconf.const import *
import shutil, json
import tarfile

parser = ArgumentParser(description="A pipeline to produce sets of bed files (bedsets) from bedbase")

parser.add_argument("-q", "--JSON-query", help="path to JSON file containing the query", type=str) # path to JSON with query
parser.add_argument("-c", "--bbconfig-path", dest="bedbase_config", type=str, required=False, default=None,
                    help="path to the bedbase configuration file")
parser.add_argument("-b", "--bedset-name", help="name assigned to queried bedset", default=str )
parser.add_argument("-o", "--output-folder", help="path to folder where tar file and igd database will be stored", default=str )
#parser.add_argument("-p", "--port", help="port number to set connection to elasticsearch", type=str)

# add pypiper args to make pipeline looper compatible
parser = pypiper.add_pypiper_args(parser, groups=["pypiper", "looper"],
                                            required=["--query", "--bedset_name"])
args = parser.parse_args()

# SET OUTPUT FOLDER
# use output parent argument from looper
out_parent = args.output_parent

bbc = bbconf.BedBaseConf(filepath=bbconf.get_bedbase_cfg(args.bb_config))

def main():
    pm = pypiper.PipelineManager(name="bedbuncher-pipeline", outfolder=out_parent, args=args)
    # Open connection to the elastic cluster;
 	#try:
	# 	es = Elasticsearch([{"host":args.dbhost, "port":args.port}])
	# 	print("Connected to elasticsearch cluster", es.info())
	# except elasticsearch.ConnectionError as Connection_error:
	# 	print("Error:", Connection_error)
	# # Perform bedbase search 	
	# search_result = es.search(index=BED_INDEX, body=args.query)
	# if 'hits' in search_result and 'total' in search_result['hits'] and int(search_result['hits']['total']['value']) > 0: 
	# 	print("the query {} returned {} hits".format(args.query, search_result['hits']['total']['value']))				
	# 	# need to iterate through the returned dictionary to find files paths
	# 	search_hits = search_result['hits']['hits']
	
	def JSON_to_dict(file_name):
		with open(filename) as f_in:
			return(json.load(f_in))

	# Establish Elasticsearch connection and check status using bbconf
	bbc.establish_elasticsearch_connection()

	# Use bbconf method to look for files in the es index
	es = bbc.search_bedfiles(query=JSON_to_dict(args.JSON_query))

	# check for prescence of the output folder and create it if needed
	output_folder = os.path.dirname(args.output_file)
		if not os.path.exists(output_folder):
    		print("Output directory does not exist. Creating: {}".format(output_folder))
    		os.makedirs(output_folder)

	# Create a tar archive using the paths to the bed files provided by the bbconf search object
	# find path through es[i]['path']
	tar_archive_file = os.path.join(args.output_folder, args.bedset_name + '.tar.gz') #should provide a path instead of just the name
	tar_archive = tarfile.open(tar_archive_file, mode=w:gz) #w:gz open for gzip cmpressed writing
   	for files in es:
      	bedfile_path = files[BEDFILE_PATH_KEY][0] 
		tar_archive.add(bedfile_path, arcname=os.path.basename(bedfile_path)) # see if setting arc name is default
	tar_archive.close()

	# Create df with bedfiles metadata: gc_content, num_regions, mean_abs_tss_dist, genomic_partitions
	bedstats_df = pd.DataFrame(columns=['BEDfile_id', 'GC_Content', 'Regions_number', 'Distance_from_feature', 
						'exon_frequency', 'exon_percentage', 
						'intergenic_frequency', 'intergenic_percentage',
						'intron_frequency', 'intron_percentage',
						'promoterCore_frequency', 'promoterCore_percentage'
						'promoterProx_frequency', ' promoterProx_percentage'])
	
	# transform individual stats from dictionary into floats to perform calculations
	def make_float(es_element):
		float(es_element[0])

	# How to access elements in search object produced by bbconf.search_index : es[i]['path']

	# iterate through the ['hits']['hits']['_source'] attribute of the bedset
	for bed_meta in es:
		file_id = bed_meta["id"] # 'id': ['3']
		gc_content = bed_meta["GC_content"]
		regions_number = bed_meta["number_of_regions"]
		feat_distance = bed_meta["mean_absolute_TSS_distance"]		
		bedstats_df = bedstats_df.append({'BEDfile_id': file_id, 
						'GC_Content': make_float(gc_content), 
						'Regions_number': make_float(regions_number),
						'Distance_from_feature':make_float(feat_distance)},
						ignore_index=True)
		# iterate through the genomic partitions list to get features like exon and intron with stats 
		genomic_partitions = bed_meta["genomic_partitions"]
		for item in genomic_partitions:
			partition_id = item["partition"]
			partition_frequency = item["Freq"]
			partition_percent = item["Perc"] 
			bedstats_df = bedstats_df.append({partition_id + '_frequency': partition_frequency, 
							partition_id + '_percentage': partition_percent}, ignore_index=True)		
			
	# Calculate mean and stdv statistics
	avg_stats = bedstats_df.mean(axis=0) 
	stdv_stats = bedstats_df.std(axis=0)
	avg_dictionary = dict(avg_stats)
	stdv_dictionary = dict(stdv_stats) 
	
	# Save bedstats_df as csv file and put it into the user defined output_folder	
	bedset_stats = os.path.join(args.output_folder, args.bedset_name + '.csv')
	bedstats_df.to_csv(bedset_stats, index=False)

		
	# CREATE THE IGD DATABASE
	# Create a .txt file with the paths to the queried bed files as input to igd the command
	txt_bed_path = os.path.join(args.output_folder, args.bedset_name + '.txt')
	txt_file = open(txt_bed_path, "wr") # need to double check mode
	for files in es:
		bedfile_path = files[BEDFILE_PATH_KEY][0]
		txt_file.write(bedfile_path)
	txt_file.close()
	pm.clean_add(txt_bed_path)

	# define CML template to create iGD database
	igd_template = "igd create {bed_source_path} {igd_folder_path} {database_name}" # put contents into igd folder, zip and provide path to zipped file.
	igd_folder_name = args.bedset_name + "_igd"
	igd_folder_path = os.path.join(args.output_folder, igd_folder_name)
	os.makedirs(igd_folder_path)
	print("Directory {} succesfully created".format(igd_folder_name))

	cmd = igd_template.format(bed_source_path=txt_bed_path, igd_folder_path=igd_folder_path, database_name=args.bedset_name)
	pm.run(cmd, target=os.path.join(igd_folder_path, args.bedset_name + ".igd"))
	
	# create a nested dictionary with avgs,stdv, and paths to tar archives, bedset csv file and igd database. 
	bedset_summary_info = {'bedset_means': avg_dictionary, 'bedset_stdv': stdv_dictionary, "tar_archive": [tar_archive_path],	
							'bedset_df': [bedset_stats], 'igd_path': [igd_folder_path]}
	
	# Insert bedset information into BEDSET_INDEX
	bbc.insert_bedsets_data(data=bedset_summary_info)
	print("{} summary info was succesfully inserted into the {}".format(args.bedset_name, BEDSET_INDEX))


	pm.stop_pipeline()


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Pipeline aborted.")
        sys.exit(1)




		

