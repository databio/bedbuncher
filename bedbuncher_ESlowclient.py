#!/usr/bin/env python
from argparse import ArgumentParser
import pypiper, os, sys
import pandas as pd
from elasticsearch import Elasticsearch
from elasticsearch import helpers
from elasticsearch_dsl import MultiSearch, Search
import tarfile

# steps needed to create the pipeline
# 1. Runs the filter to create a list of a subset of BEDfiles in the database. These will be the bed files included in the bedset. 
# 2. tars the subset to make an archive for this bedset
# 3. Reads all the metadata from each bed in the bedset and creates a matrix with all statistics produced by bedstat for individual files. We then average these statistics 
# 4. Inserts the bedset statistics into the bedbase under an index for bedsets

parser = ArgumentParser(
    description="A pipeline to produce sets of bed files (bedsets) from bedbase")

parser.add_argument("-i", "--index", help="name of ES database index where files are stored", type=str)
parser.add_argument("-q", "--query", help="what variable to perform to search in", type=str)
parser.add_argument("-d", "--dbhost", help="this should be the database host address we need to connect to", default="localhost" )
parser.add_argument("-p", "--port", help="port number to set connection to elasticsearch", type=str)

# add pypiper args to make pipeline looper compatible
parser = pypiper.add_pypiper_args(parser, groups=["pypiper", "looper"],
                                            required=["--index", "--query"])
args = parser.parse_args()

# SET OUTPUT FOLDER
# use output parent argument from looper
out_parent = args.output_parent


def main():
    pm = pypiper.PipelineManager(name="bedbuncher-pipeline", outfolder=out_parent, args=args)

# ESTABLISH CONTACT WITH ELASTIC SEARCH AND PROVIDE THE FILTER CRITERIA TO GET THE SET OF BED FILES
# bed files were already inserted into the database using es.index
# Perform the search/retrieval using the low-level elasticsearch python client

# Open connection to the elastic cluster;
	try:
		es = Elasticsearch([{"host":args.dbhost, "port":args.port}])
		if es is not None:
			print("Connected to elasticsearch cluster", es.info())
			search_result = es.search(index=args.index, body=args.query)
			if not search_result['timed_out']: 
				print("the query {} returned {} hits".format(args.query, search_result['hits']['total']['value']))				
				# need to iterate through the returned dictionary to find files paths
				search_hits = search_result['hits']['hits']
				# Alternative to tar using the CML: TAR files using the tarfile module
	            tar_archive = tarfile.open(args.bedset_name, mode=w:gz)
	            for files in search_hits:
	                # need to get access to bed json file to get the paths ['_source']
	            	bed_source = files['_source'].get("sample_name") # this should point me to the raw bed file
					tar_archive.add(bed_source)
				tar_archive.close()
			else:
				print("The provided query doesn't match the database record")
	except Exception as ex:
		print("error while setting connection", ex) 

	# Create df with bedfiles metadata: gc_content, num_regions, mean_abs_tss_dist, genomic_partitions
	bedstats_df = pd.DataFrame(columns=['BEDfile_id', 'GC_Content', 'Regions_number', 'Distance_from_feature', 
						'Exon_frequency', 'Exon_percentage', 
						'Intergenic_frequency', 'Intergenic_percentage',
						'Intron_frequency', 'Intron_percentage',
						'PromoterCore_frequency', 'PromoterCore_percentage'
						'PromoterProx_frequency', ' PromoterProx_percentage'])
	
	# transform individual stats from dictionary into floats to perform calculations
	def make_float(es_element):
		float(es_element[0])

	# iterate through the ['hits']['hits']['_source'] attribute of the bedset
	for bed_file in search_hits:
		source = bed_file['_source']
		# get GenomicDIstributions data for each bed file as described in JSON file
		file_id = source.get("id")
		gc_content = source.get("gc_content")
		regions_number = source.get("num_regions")
		feat_distance = source.get("mean_abs_tss_dist")		
		bedstats_df = bedstats_df.append({'BEDfile_id': make_float(file_id), 
						'GC_Content': make_float(gc_content), 
						'Regions_number': make_float(regions_number),
						'Distance_from_feature':make_float(feat_distance)},
										ignore_index=True)
		# iterate through the genomic partitions list to get feaures like exon and intron with stats (each list has several dictionaries)
		genomic_partitions = source.get("genomic_partitions")
		for item in genomic_partitions:
			partition_id = item.get("partition")
			partition_frequency = item.get("Freq")
			partition_percent = item.get("Perc")
			if partition_id == 'exon':
			    bedstats_df = bedstats_df.append({'Exon_frequency': partition_frequency, 
								'Exon_percentage': partition_percent}, ignore_index=True)		
			elif partition_id == 'intergenic': 
				bedstats_df = bedstats_df.append({'Intergenic_frequency': partition_frequency, 
								'Intergenic_percentage': partition_percent}, ignore_index=True)				
			elif partition_id == 'intron':
				bedstats_df = bedstats_df.append({'Intron_frequency': partition_frequency, 
								'Intron_percentage': partition_percent}, ignore_index=True)
			elif partition_id == 'promoterCore':
				bedstats_df = bedstats_df.append({'PromoterCore_frequency': partition_frequency, 
								'PromoterCore_percentage': partition_percent}, ignore_index=True)
			else:
				partition_id == 'promoterProx':
				bedstats_df = bedstats_df.append({'PromoterCore_frequency': partition_frequency, 
								'PromoterCore_percentage': partition_percent}, ignore_index=True)

	# Calculate average values for gc_content, num_regions and mean_abs_tss_dist 
	avg_gc_content = bedstats_df["GC_Content"].mean()
	avg_num_regions = bedstats_df["Regions_number"].mean()
	avg_distance = bedstats_df["Distance_from_feature"].mean()
		

