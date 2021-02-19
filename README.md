# bedbuncher
Pipeline designed to create **bedsets** (sets of BED files) that will be retrieved from bedbase.

Examples of bedsets:
- BED files with regions count greater than 1000 (query: `regions_no>%s` query_val: `1000`)
- BED files produced in ChiP-seq experiments

## To run the pipeline 
1. Clone the repository
2. Install required python packages via 
```
pip install -r requirements/requirements-all.txt --user
```
3. Install additional dependencies:

- `iGD` -- a command line executable, which builds a database that integrates genomic sets from one or more data sources and minimizes the search space for a specific query. Visit the [`iGD` repository](https://github.com/databio/iGD) for more information and installation details.

4. Run a [PostgreSQL](https://www.postgresql.org/) database instance

For example, to run an instnace in a Docker container use:
```
docker run -d --name bedbase-postgres -p 5432:5432 -e POSTGRES_PASSWORD=bedbasepassword -e POSTGRES_USER=postgres -e POSTGRES_DB=postgres -v postgres-data:/var/lib/postgresql/data postgres
```
The DB login credentials need to match what's specified in the bedbase configuration file.

5. Submit the pipeline with [`looper`](https://looper.readthedocs.io/en/latest/)

**Input:** A [PEP](http://pep.databio.org/en/latest/) with one sample per bedset. Each sample needs to have an attribute specifying a query that will restrict the database search to create a bedset.

```
looper run project/cfg.yaml
```

## Pipeline outputs
`bedbuncher` inserts bedstet metadata into the PostgreSQL database and generates the following files:
- TAR ball containing the BED files that match the query criteria.
- Dataframe where rows represent individual BED files and columns show statistics generated by GenomicDistributions (for ease of user-needed calculations).
- iGD database created from the bedset.
- Bedset statistics (currenty means and standard deviations).
- PEP for a specific bedset created using the pipeline.
- trackHub directory for the BED set that can be viewed on the UCSC Genome Browser. 


