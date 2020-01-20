#!/usr/bin/env python
from argparse import ArgumentParser
import pypiper, os, sys
from elasticsearch import Elasticsearch
from elasticsearch import helpers
from elasticsearch_dsl import Search
from elasticsearch.serializer import JSONSerializer
import json, gzip, shutil, yaml
import tarfile

# steps needed to create the pipeline
# 1. Runs the filter to create a list of a subset of BEDfiles in the database. These will be the bed files included in the bedset. 
# 2. tars the subset to make an archive for this bedset
# 3. Reads all the metadata from each bed in the bedset and creates a matrix with all statistics produced by bedstat for individual files. We then average these statistics 
# 4. Inserts the bedset statistics into the bedbase under an index for bedsets

parser = ArgumentParser(
    description="A pipeline to produce sets of bed files (bedsets) from bedbase")

parser.add_argument("-f", "--field", help="what variable to perform to search in", type=str)
parser.add_argument("-v", "--value", help="key value we're looking for e.g protocol==ATAC-seq", type=str)
parser.add_argument("-d", "--dbhost", help="this should be the database host address we need to connect to", default="localhost" )


# add pypiper args to make pipeline looper compatible
parser = pypiper.add_pypiper_args(parser, groups=["pypiper", "looper"],
                                            required=["--protocol"])
args = parser.parse_args()

# SET OUTPUT FOLDER
# use output parent argument from looper
out_parent = args.output_parent

pm = pypiper.PipelineManager(name="bedbuncher-pipeline", outfolder=out_parent, args=args)


# ESTABLISH CONTACT WITH ELASTIC SEARCH AND PROVIDE THE FILTER CRITERIA TO GET THE SET OF BED FILES
#  In Elasticsearch you index, search,sort and filter documents.

# Open connection to the elastic cluster;
try:
	es = Elasticsearch([{"host": args.dbhost}, "port":9200])
	print("Connected to elasticsearch cluster", es.info())

	# bed files were already inserted into the database using es.index(index="bedstat_bedfiles", body=data)
	# now we need to retrieve the files using the filter criteria 
	#files_set = es.search(index = "bedstat_bedfiles", body={'query':{'match': {args.field: args.value}}}) 
	
	#Alternatively perform the search/retrieval using the high-level python elasticsearch_dsl
	s = Search(using=es, index="bedstat_bedfiles") 
	s = s.query('match', args.field=args.value)
	file_count = s.count()
	print("the query {} returned {} files".format(args.value, file_count))
	retrieved_files = s[0:file_count].execute()

	# Need to TAR the files retrieved by ES using the CML
	tar_files = "tar -czvf bedset.tar.gz {bedset}"
	cmd = tar_files.format(bedset=retrieved_files) 
	pm.run(cmd, target=target)

	# ALternative: TAR files using the tarfile module
	tf = tarfile.open("bedset.tar.gz, mode="w:gz)
	bed_files = os.listdir(".")
	for f in bed_files
		tf.add(f)
	tf.close()

except Exception as ex:
	print("error while setting connection", ex) 

