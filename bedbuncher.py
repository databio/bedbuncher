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

parser.add_argument("-q", "--JSON-query-path", help="path to JSON file containing the query", type=str) # path to JSON with query
parser.add_argument("-c", "--bbconfig", type=str, required=False, default=None,
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

def main():
	pm = pypiper.PipelineManager(name="bedbuncher-pipeline", outfolder=out_parent, args=args)
	
	def JSON_to_dict(file_name):
		with open(file_name) as f_in:
			return(json.load(f_in))

	# Establish Elasticsearch connection and check status using bbconf
	bbc.establish_elasticsearch_connection()

	# Use bbconf method to look for files in the es index
	es = bbc.search_bedfiles(query=JSON_to_dict(args.JSON_query_path))

	# check for prescence of the output folder and create it if needed
	output_folder = os.path.dirname(args.output_folder)
	if not os.path.exists(output_folder):
		print("Output directory does not exist. Creating: {}".format(output_folder))
		os.makedirs(output_folder)

	
	# Determine if search object has duplicate files and remove them to avoid downstream complications.  
	# es = [i for n, i in enumerate(es) if i not in es[:n]]
	def remove_duplicates(in_list):
		new_list = []
		for i in in_list:
			if i not in new_list:
				new_list.append(i)
		return(new_list)

	if len([j for j in es if es.count(j) > 1]) != 0:
		es = remove_duplicates(es)

	# Create a tar archive using the paths to the bed files provided by the bbconf search object
	tar_archive_file = os.path.join(args.output_folder, args.bedset_name + '.tar.gz') 
	tar_archive = tarfile.open(tar_archive_file, mode="w:gz") #w:gz open for gzip cmpressed writing
   	
	for files in es:
   		bedfile_path = files[BEDFILE_PATH_KEY] 
   		tar_archive.add(bedfile_path, arcname=os.path.basename(bedfile_path), recursive=False, filter=None) 
	tar_archive.close()

	# Create df with bedfiles metadata: gc_content, num_regions, mean_abs_tss_dist, genomic_partitions
	bedstats_df = pd.DataFrame(columns=['BEDfile_id', 'GC_Content', 'Regions_number', 'Distance_from_feature', 
						'exon_frequency', 'exon_percentage', 
						'intergenic_frequency', 'intergenic_percentage',
						'intron_frequency', 'intron_percentage',
						'promoterCore_frequency', 'promoterCore_percentage',
						'promoterProx_frequency', 'promoterProx_percentage'])
	
	# transform individual stats from dictionary into floats to perform calculations, if needed
	def make_float(es_element):
		float(es_element[0])

	# Access elements in search object produced by bbc.search
	for bed_meta in es:
		file_id = bed_meta["id"][0] # 'id': ['3']
		gc_content = bed_meta["gc_content"][0]
		regions_number = bed_meta["num_regions"][0]
		feat_distance = bed_meta["mean_abs_tss_dist"][0]		
		bedstats_df = bedstats_df.append({'BEDfile_id': file_id, 
						'GC_Content': gc_content, 
						'Regions_number': regions_number,
						'Distance_from_feature':feat_distance},
						ignore_index=True)
		# iterate through the genomic partitions list to get features like exon and intron with stats 
		genomic_partitions = bed_meta["genomic_partitions"]
		for item in genomic_partitions:
			partition_id = item["partition"]
			partition_frequency = item["Freq"]
			partition_percent = item["Perc"] 
			bedstats_df = bedstats_df.append({partition_id + '_frequency': partition_frequency, 
							partition_id + '_percentage': partition_percent}, ignore_index=True)		
			
	# Calculate bedset statistics
	avg_stats = bedstats_df.mean(axis=0) 
	stdv_stats = bedstats_df.std(axis=0)
	avg_dictionary = dict(avg_stats)
	#print(avg_dictionary)
	stdv_dictionary = dict(stdv_stats) 
	#print(stdv_dictionary)
	
	# Save bedstats_df as csv file into the user defined output_folder	
	bedset_stats = os.path.join(args.output_folder, args.bedset_name + '.csv')
	bedstats_df.to_csv(bedset_stats, index=False)
	
	# IGD DATABASE
	# Need a .txt file with the paths to the queried bed files as input to the igd create command
	txt_bed_path = os.path.join(args.output_folder, args.bedset_name + '.txt')
	txt_file = open(txt_bed_path, "a") 
	for files in es:
		bedfile_path = files[BEDFILE_PATH_KEY]
		print(bedfile_path)
		txt_file.write("{}\r\n".format(bedfile_path))
	txt_file.close()
	pm.clean_add(txt_bed_path)

	# iGD database
	igd_folder_name = args.bedset_name + "_igd" 
	igd_folder_path = os.path.join(args.output_folder, igd_folder_name)
	os.makedirs(igd_folder_path)
	print("Directory {} succesfully created".format(igd_folder_name))

	# create a temp file to untar bed file, use them to create the igd database and then clean them 
	#bed_temp = args.bedset_name + "_temp"
	#bed_temp_path = os.path.join(args.output_folder, bed_temp)
	#os.makedirs(bed_temp_path)

	# Command templates for IGD database construction
	igd_template = "igd create {bed_source_path} {igd_folder_path} {database_name} -f"
	gzip_template = "gzip -r {dir}"
	cmd1 = igd_template.format(bed_source_path=txt_bed_path, igd_folder_path=igd_folder_path, database_name=args.bedset_name)
	cmd2 = gzip_template.format(dir=igd_folder_path)
	cmd = [cmd1, cmd2]
	pm.run(cmd, target=os.path.join(igd_folder_path, args.bedset_name + ".igd"))
	
	# create a nested dictionary with avgs,stdv, paths to tar archives, bedset csv file and igd database. 
	bedset_summary_info = {'bedset_means': avg_dictionary, 'bedset_stdv': stdv_dictionary, "tar_archive_path": [tar_archive_file],	
							'bedset_df': [bedset_stats], 'igd_database__path': [igd_folder_path]}
	
	print("{} summary info and additional output files:{}".format(args.bedset_name, bedset_summary_info))
	
	# Insert bedset information into BEDSET_INDEX
	bbc.insert_bedsets_data(data=bedset_summary_info)
	print("{} summary info was succesfully inserted into the {} =D".format(args.bedset_name, BEDSET_INDEX))


	pm.stop_pipeline()


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Pipeline aborted.")
        sys.exit(1)




		

