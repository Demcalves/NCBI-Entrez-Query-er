import requests
import time
import xml.etree.ElementTree as ElementTree
import json 
from Bio import Entrez

""" This is an elaborate get request from NCBI, might be useful for querying faster results. Best practice would be to validate your results.
    A user should have an API-key from NCBI before initializing this, a config.json file with the users api-key can be stored in there
    to call in this function (from its root directory).

    Part of having an API key is because NCBI will throttle web requests to only 3 (Assuming the databases are cached), so the returns in
    theory are instant depending on your internet speed, however doing 3+ requests will cut your connection (HTTP Error 429). Having an 
    API key with NCBI extends your max call requests to 10 calls in an instant! In addition, by using the time package, you can thread 
    your calls so that it happens within a timeframe that is acceptable to the host server. 

    Headers are added to the get request as an additional precaution in case your client is flagged for any reason.

    To understand what is being queried: a user of this script uses a gene_id code (ex. Gene_Id 630 refers to the Human AQP7 gene)
    From this gene_id inquiry, the NCBI database will return all information referenced to the gene as metadata, with 
    corresponding elements that act as keys with value points like in JSON objects. From this metadata, the script then
    looks for all uid (protein) elements within the dataset. At this stage uid is a key for the database to handle the object, but
    other databases like RefSeq or Ensembl will also have their own syntax and id for that exact same uid. So now that we have a uid,
    we can launch another query to get the metadata tied to it, adjusting search criteria for whatever sequence desired. This script 
    currently looks for RefSeq (i.e NP..., XP...) ids, returning the first one's fasta sequence that it finds. 

    An additional quality of life that has been added is a simple write method to add the fasta sequence to a text file, saved in fasta format.

    In addition, I included a method that 

"""
"""  
    gene name is a string object, know the name of your gene!
    queryGeneIds returns a list of ten species identified in the gene.
    We are using homologene to return gene ids for a orthologs of the 
    gene we are trying to compare.

    You will need to download biopyton from your command line (cmd, powershell or linux shell)
    type:

    pip install biopython

    Im not familiar with Mac interfaces 
"""

class Multiplex_FASTA_Extractor:

    def __init__(self, gene_name, seq_count):
        self.gene_name = gene_name
        self.seq_count = seq_count
    """  
        gene name is a string object, know the name of your gene!
        queryGeneIds returns a list of ten species identified in the gene.
        We are using homologene to return gene ids for a orthologs of the 
        gene we are trying to compare.
    """

    def getGeneNameResponse(self, gene_name):
        try:
            handle = Entrez.esearch(db="gene", term=gene_name, retmax=1)
            record = Entrez.read(handle)
            return record["IDList"]
        except Exception as e:
            print("Gene name was entered incorrectly and resulted in an error response, please reenter your gene name")
        return None

    def queryGeneIds(self, gene_name, seq_count):
        search_term = f"{gene_name}[Gene Name]"
        handle = Entrez.esearch(db="homologene", term=search_term, retmax=seq_count)  # Limit to 10 results, can change output
        record = Entrez.read(handle)
        homologene_ids = record["IdList"]
        return homologene_ids
    #generates a fasta to write results
    def fetchFastaFromEntrez(self, gene_name, seq_count, api_key=None, db_from='gene', db_to='protein', linkname='gene_protein'):
        gene_id_list = self.queryGeneIds(gene_name, seq_count)
        base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
        fasta_results = {}
        file_name = input("enter your_file_name: ")
        
        # CHANGE THIS to where in your directory you want to save your fasta file
        with open(f"C:/Users/<Pathway to where you want to save your fasta/{file_name}.fasta", "w") as fasta_file: 


            headers = {'User-Agent': 'Mozilla/5.0'}  # Avoid being flagged as a bot when sending get request query
            job_count = 0
            for gene_id in gene_id_list:
                retry_count = 3  # Max retries for 429 errors
                job_count += 1
                while retry_count > 0:
                    # Step 1: Use ELink to get linked protein ID for the gene
                    elink_url = f"{base_url}elink.fcgi?dbfrom={db_from}&db={db_to}&id={gene_id}&linkname={linkname}&cmd=neighbor_history"
                    if api_key:
                        elink_url += f"&api_key={api_key}"

                    response = requests.get(elink_url, headers=headers)
                    if response.status_code == 429:
                        print(f"Rate limit hit for Gene ID {gene_id}. Retrying after delay...")
                        time.sleep(5)  # Exponential backoff
                        retry_count -= 1
                        continue  # Retry the request
                    elif response.status_code != 200:
                        print(f"ELink request failed for Gene ID {gene_id}: {response.status_code}")
                        break  # Stop retries for this gene

                    root = ElementTree.fromstring(response.text)
                    web_env = root.find(".//WebEnv")
                    query_key = root.find(".//QueryKey")

                    if web_env is None or query_key is None:
                        print(f"Failed to retrieve WebEnv or QueryKey for Gene ID {gene_id}")
                        break

                    web_env = web_env.text
                    query_key = query_key.text

                    # Step 2: Use ESummary to get metadata for the linked protein
                    esummary_url = f"{base_url}esummary.fcgi?db={db_to}&query_key={query_key}&WebEnv={web_env}" + f"&api_key={api_key}"

                    summary_response = requests.get(esummary_url, headers=headers)
                    if summary_response.status_code == 429:
                        print(f"Rate limit hit for Gene ID {gene_id} (ESummary). Retrying after delay...")
                        time.sleep(5)
                        retry_count -= 1
                        continue
                    elif summary_response.status_code != 200:
                        print(f"ESummary request failed for Gene ID {gene_id}: {summary_response.status_code}")
                        break

                    root_summary = ElementTree.fromstring(summary_response.text)
                    protein_id = None
                    

                    # Step 3: Use EFetch to get the FASTA sequence
                    # of all the protein ids that are generated we need to check there values
                    
                    # the webaddress .//Docsum/{address} has the following items that you can extract
                    # Protein Id as shown. Other ones are xml objects (Item) that need their name specified in the search:
                    # TaxId, Organism, AccesionVersion, Title (Definition), MolecularType, SourceDB, Length, Chromosome (for some entries)
                    #since we are working with a particular gene, we can get the speices name as an example below
                    #species_name = root_summary.find(".//DocSum/Item[@Name='Organism']")

                    protein_ids = [uid.text for uid in root_summary.findall(".//DocSum/Id") if uid is not None]
                    #print(protein_ids) print statement for trouble shooting or if you want to visualize just how many protein ids there are.

                    for protein_id in protein_ids:
                        if protein_id is None:
                            print(f"No protein ID found for Gene ID {gene_id}")
                            break
                        else:
                            efetch_url = f"{base_url}efetch.fcgi?db={db_to}&id={protein_id}&rettype=fasta&retmode=text"+f"&api_key={api_key}"
                            efetch_response = requests.get(efetch_url, headers=headers)

                            if efetch_response.status_code == 429:
                                print(f"Rate limit hit for Gene ID {gene_id} (EFetch). Retrying after delay...")
                                time.sleep(5)
                                retry_count -= 1
                                continue
                            elif efetch_response.status_code == 200:
                                fasta_sequence = efetch_response.text
                                # Prioritize the canonical refseq NP <-> XP.
                                if fasta_sequence[:3] == ">NP" or fasta_sequence[:3] ==">XP":
                                    fasta_results[gene_id] = fasta_sequence
                                    break # Success, move to the next gene, rather than iterate over other refseq or experimental isoforms
                                    #Python is interesting since its a dynamic language, inner scoped variables can be globally accessed. Java would have deallocated local stack memory of this variable
                                else:
                                    fasta_results[gene_id] = fasta_sequence #update and continue iterating
                            else:
                                print(f"Failed to fetch FASTA for Gene ID {gene_id}: {efetch_response.status_code}")
                                break
                    #break the loop with our fasta_sequence value derived from the iterative loop
                    if fasta_sequence[:3] == ">NP" or fasta_sequence[:3] ==">XP":
                        fasta_file.write(f"{fasta_sequence}") #writes fastas that can be opened with a text editor
                        break # Success, move to the next gene
                print(f"Job done for Gene Id: {gene_id}, {job_count} / {len(gene_id_list)} completed \n" ) #for visualization in the terminal 
                print(f"Gene ID species:")
                print(f"{fasta_sequence}\n")      
                time.sleep(0.35)  # Ensure we don't exceed 3 requests per second
        fasta_file.close()
        return fasta_results

# Example usage with an API key
# Load API key from JSON file
# the filepath should be wherever the API key is stored
# example: replace <element> entries with appropriate directory names
# C:/Users/Sally/.../<ncbi_config.json>

# CHANGE THIS
Entrez.email = "your_email@wherever.com"

# CHANGE THIS TO WHEREVER YOU SAVED YOUR KEY
with open("C:/Users/<your_user_name>/.../config.json", "r") as config_file:
    config = json.load(config_file)

api_key = config.get("NCBI_API_KEY")

if not api_key:
    raise ValueError("NCBI API key not found in config.json!")
gene_query_response = None
while gene_query_response == None:
    gene_query = input("Enter exact gene name: ").upper()
    gene_query_response = Multiplex_FASTA_Extractor.getGeneNameResponse(gene_query)

# if you enter a non integer primitive the script will crash :<
try:
    seq_count = int(input("Enter number of seqs (probably don't go over 20): "))
except TypeError as e:
    print("Non integer number added, default script will run 10 sequences")

#script can still crash if you enter your Gene_ID incorrectly
fasta_data = Multiplex_FASTA_Extractor.fetchFastaFromEntrez(gene_query, seq_count=10, api_key=api_key)