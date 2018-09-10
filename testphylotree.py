from Bio import Entrez
from Bio import SeqIO
import xml
import os
import sys

# set contact email for the whole project
Entrez.email = "nemo0806@gmail.com"

#start a handle
handle = Entrez.einfo()
print(handle)
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


# start a query handle to the database

queryhandle = Entrez.egquery(term="Opuntia AND rpl16")
queryrecord = Entrez.read(queryhandle)
# count the number of returned results for checking
for row in queryrecord["eGQueryResult"]:
    if row["DbName"] == "nuccore":
        print(row["Count"])

# send request to fetch ID from title input

searchhandle = Entrez.esearch(
    db="nuccore", term="Opuntia AND rpl16", usehistory="y")
searchrecord = Entrez.read(handle)
gi_list = searchrecord["IdList"]

# search using fetched IDlist, and save history for fetching
search_results = Entrez.read(Entrez.epost("pubmed", id=",".join(gi_list)))
webenv = search_results["WebEnv"]
query_key = search_results["QueryKey"]

#fetching a sequence via its ID

filename = "*.gbk"
seqid = ''

if not os.path.isfile(filename):
    # Downloading...
    net_handle = Entrez.efetch(
        db="nucleotide", id=seqid, rettype="gb", retmode="text")
    out_handle = open(filename, "w")
    out_handle.write(net_handle.read())
    out_handle.close()
    net_handle.close()
    print("Saved")

print("Parsing...")
record = SeqIO.read(filename, "genbank")
print(record)

#

# read a sequence
# for seq_record in SeqIO.parse("ls_orchid.fasta", "fasta"):
#     print seq_record.id
#     print repr(seq_record.seq)
#     print len(seq_record)
