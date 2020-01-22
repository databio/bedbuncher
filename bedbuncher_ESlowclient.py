#!/usr/bin/env python
from argparse import ArgumentParser
import pypiper, os, sys
from elasticsearch import Elasticsearch
from elasticsearch import helpers
from elasticsearch_dsl import MultiSearch, Search
import json, gzip, shutil, yaml
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
                                            required=["--field", "--value"])
args = parser.parse_args()

# SET OUTPUT FOLDER
# use output parent argument from looper
out_parent = args.output_parent


def main():
	pm = pypiper.PipelineManager(name="bedbuncher-pipeline", outfolder=out_parent, args=args)

# ESTABLISH CONTACT WITH ELASTIC SEARCH AND PROVIDE THE FILTER CRITERIA TO GET THE SET OF BED FILES
#  In Elasticsearch you index, search,sort and filter documents.
# bed files were already inserted into the database using es.index
# Perform the search/retrieval using the low-level python elasticsearch_dsl

# Open connection to the elastic cluster;
	try:
		es = Elasticsearch([{"host":args.dbhost, "port":args.port}])
		if es is not None:
			print("Connected to elasticsearch cluster", es.info())
			search_result = es.search(index=args.index, body=args.query)
			if search_result.success() == True: 
				print("the query {} returned {} hits".format(args.query, search_result['hits']['total']['value']))
				
				# need to iterate through the returned dictionary to find files paths
				search_hits = search_result['hits']['hits']

				# Alternative to tar using the CML: TAR files using the tarfile module
	            tar_archive = tarfile.open(args.bedset_name, mode=w:gz)
	            for files in search_hits:
	            	bed_source = files['_source'].get("sample_name")
					tar_archive.add(bed_source)
				tar_archive.close()
			else:
				print("The provided query doesn't match the database record")

	except Exception as ex:
		print("error while setting connection", ex) 

	# Need to create matrix with bedfiles metadata
