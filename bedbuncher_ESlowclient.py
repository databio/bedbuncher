#!/usr/bin/env python
from argparse import ArgumentParser
import pypiper, os, sys
import pandas as pd
from elasticsearch import Elasticsearch
from elasticsearch import helpers
from bbconf.const import *
import tarfile

# steps needed to create the pipeline
# 1. Runs the filter to create a list of a subset of BEDfiles in the database. These will be the bed files included in the bedset. 
# 2. tars the subset to make an archive for this bedset
# 3. Reads all the metadata from each bed in the bedset and creates a matrix with all statistics produced by bedstat for individual files. We then average these statistics 
# 4. Inserts the bedset statistics into the bedbase under an index for bedsets

parser = ArgumentParser(description="A pipeline to produce sets of bed files (bedsets) from bedbase")

parser.add_argument("-q", "--query", help="what variable to perform to search in", type=dict)
parser.add_argument("-d", "--dbhost", help="this should be the database host address we need to connect to", default="localhost" )
parser.add_argument("-b", "--bedset-name", help="name assigned to queried bedset", default=str )
#parser.add_argument("-f", "--tar-folder", help="name of output folder to store tar bedset", default=str)
#parser.add_argument("-r", "--raw-folder", help="name of output folder for raw bed files", default=str)
parser.add_argument("-p", "--port", help="port number to set connection to elasticsearch", type=str)

# add pypiper args to make pipeline looper compatible
parser = pypiper.add_pypiper_args(parser, groups=["pypiper", "looper"],
                                            required=["--query", "--bedset_name"])
args = parser.parse_args()

# SET OUTPUT FOLDER
# use output parent argument from looper
out_parent = args.output_parent


def main():
    pm = pypiper.PipelineManager(name="bedbuncher-pipeline", outfolder=out_parent, args=args)
    # Open connection to the elastic cluster;
    try:
		es = Elasticsearch([{"host":args.dbhost, "port":args.port}])
		print("Connected to elasticsearch cluster", es.info())
	except elasticsearch.ConnectionError as Connection_error:
		print("Error:", Connection_error)
	# Perform bedbase search 	
	search_result = es.search(index=BED_INDEX, body=args.query)
	if 'hits' in search_result and 'total' in search_result['hits'] and int(search_result['hits']['total']['value']) > 0: 
		print("the query {} returned {} hits".format(args.query, search_result['hits']['total']['value']))				
		# need to iterate through the returned dictionary to find files paths
		search_hits = search_result['hits']['hits']
		
		# Alternative to tar using the CML: TAR files using the tarfile module
        tar_archive = tarfile.open(args.bedset_name + '.tar.gz', mode=w:gz) #w:gz open for gzip cmpressed writing
       	for files in search_hits:
            # need to get access to bed json file to get the paths ['_source']
        	bedfile_path = files['_source']['_id'] # path should be sourced by bedstat?
			tar_archive.add(bedfile_path, arcname=os.path.basename(bedfile_path)) 
		tar_archive.close()
	else:
		raise elasticsearch.NotFoundError("The provided query doesn't match the database record")
	

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

	# iterate through the ['hits']['hits']['_source'] attribute of the bedset
	for bed_file in search_hits:
		source = bed_file['_source']
		# get GenomicDIstributions data for each bed file as described in JSON file
		file_id = source["id"] # 'id': ['3']
		gc_content = source["gc_content"]
		regions_number = source["num_regions"]
		feat_distance = source["mean_abs_tss_dist"]		
		bedstats_df = bedstats_df.append({'BEDfile_id': file_id, 
						'GC_Content': make_float(gc_content), 
						'Regions_number': make_float(regions_number),
						'Distance_from_feature':make_float(feat_distance)},
										ignore_index=True)
		# iterate through the genomic partitions list to get feaures like exon and intron with stats (each list has several dictionaries)
		genomic_partitions = source["genomic_partitions"]
		for item in genomic_partitions:
			partition_id = item["partition"]
			partition_frequency = item["Freq"]
			partition_percent = item["Perc"] 
			bedstats_df = bedstats_df.append({partition_id + '_frequency': partition_frequency, 
								partition_id + '_percentage': partition_percent}, ignore_index=True)		
			
	# Calculate mean and stdv statistics
	# axis = 0 calculates the column wise mean of the dataframe
	avg_stats = bedstats_df.mean(axis=0) 
	stdv_stats = bedstats_df.std(axis=0)
	avg_dictionary = dict(avg_stats)
	std_dictionary = dict(stdv_stats) 

	# create a nested dictionary with avgs and std values 
	bedset_stats = {'bedset_means': avg_dictionary, 'bedset_stdv': std_dictionary}
	# create BEDSET_INDEX
	try:
		es.index(index=BEDSET_INDEX, body=bedset_stats) # index name will be sourced from bbconf
		print("{} was succesfully created".format(BEDSET_INDEX))
	except elasticsearch.ElasticsearchException as ind_ex:
		print("Error: Bedset index could not be created", ind_ex)

	# Pending: 
	#	convert bedset_stats into a JSON file
	#	create igd database for files in the bedset


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Pipeline aborted.")
        sys.exit(1)




		

