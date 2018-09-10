# the first debugged version of our semi-automatic search and fetch program

from Bio import Entrez
from Bio import SeqIO
import xml
import os
import sys
import io

# set contact email for the whole project
Entrez.email = "nemo0806@gmail.com"

#start a handle
handle = Entrez.einfo()

#fetch contents from Entrez database as xml file
einfo = handle.read()  #einfo is str

with open('einfo.xml', 'w', encoding='utf-8') as f:
    f.write(einfo)

contents = Entrez.read(open('einfo.xml', 'r', encoding='utf-8'))

# EPost raw realization

# id_list = []
# search_results = Entrez.read(Entrez.epost("pubmed", id=",".join(id_list)))
# webenv = search_results["WebEnv"]
# query_key = search_results["QueryKey"]

# from cookbook, define the annotation function using EPost API


def retrieve_annotation(id_list):
    """Annotates Entrez Gene IDs using Bio.Entrez, in particular epost (to
    submit the data to NCBI) and esummary to retrieve the information.
    Returns a list of dictionaries with the annotations."""

    request = Entrez.epost("gene", id=",".join(id_list))
    try:
        result = Entrez.read(request)
    except RuntimeError as e:
        #FIXME: How generate NAs instead of causing an error with invalid IDs?
        print("An error occurred while retrieving the annotations.")
        print("The error returned was %s" % e)
        sys.exit(-1)

    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    data = Entrez.esummary(db="gene", webenv=webEnv, query_key=queryKey)
    annotations = Entrez.read(data)

    print("Retrieved %d annotations for %d genes" % (len(annotations),
                                                     len(id_list)))

    return annotations


def print_data(annotation):
    for gene_data in annotation:
        gene_id = gene_data["Id"]
        gene_symbol = gene_data["NomenclatureSymbol"]
        gene_name = gene_data["Description"]
        print("ID: %s - Gene Symbol: %s - Gene Name: %s" %
              (gene_id, gene_symbol, gene_name))


def ext(Retrieve_type):
    ''' add extention names to the out file according to Retrieve_type'''
    if Retrieve_type == 'fasta':
        return '.' + Retrieve_type
    elif Retrieve_type == 'gbwithparts':
        return '.gbk'
    else:
        raise RuntimeError('retrieve type error')

    # TODO use the einfo type to emunerate that


# # start a query handle to the database

# query_handle = Entrez.egquery(term="Opuntia AND rpl16")
# query_record = Entrez.read(query_handle)
# # count the number of returned results for checking
# for row in query_record["eGQueryResult"]:
#     if row["DbName"] == "nuccore":
#         print(row["Count"])
# query_handle.close()

# send request to fetch ID from title input

#let us encapsule our search
search_db = 'nucest'
search_term = "human AND GLUT"

#fetch step
batch_size = 5
retrieve_type = 'fasta'
retrieve_mode = 'text'
save_to_file_name = search_term

# def search(search_db, search_term, batch_size):

search_handle = Entrez.esearch(db=search_db, term=search_term, usehistory="y")
search_record = Entrez.read(search_handle)
gi_list = search_record["IdList"]
count = int(search_record["Count"])
search_handle.close()
print(count)
# search using fetched IDlist, and save history for fetching

search_results = Entrez.read(Entrez.epost(search_db, id=",".join(gi_list)))
webenv = search_results["WebEnv"]
query_key = search_results["QueryKey"]
# fetching sequences in batches via history

with open(save_to_file_name + ext(retrieve_type), "a") as f:
    for start in range(0, count, batch_size):
        end = min(count, start + batch_size)
        print("Going to download record %i to %i" % (start + 1, end))
        fetch_handle = Entrez.efetch(
            db=search_db,
            rettype=retrieve_type,
            retmode=retrieve_mode,
            retstart=start,
            retmax=batch_size,
            webenv=webenv,
            query_key=query_key)

        data = fetch_handle.read()
        fetch_handle.close()
        print(data)
        f.write(data)