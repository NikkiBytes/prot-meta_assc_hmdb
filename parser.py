import os
#import pandas as pd
import numpy as np
import xml.etree.ElementTree as ET
from lxml import etree as etree_lxml
from biothings import config
from biothings.utils.dataload import dict_convert, dict_sweep
logging = config.logger

# Association Parser    
# Association Parser    
def load_hmdb_data(data_folder):
        
    # --- Helper Methods --- 
    # --- Enter subject into document --
    def enter_subject(data, tags):
        # setup subject info       
        #uniprot_id, uniprot_name, genbank_protein_id, hgnc_id, genbank_gene_id, and gene_name.        
        uniprot_id = tags.find("{http://www.hmdb.ca}uniprot_id")
        uniprot_id = uniprot_id.text
        data["subject"]["uniprot_id"]=uniprot_id

        uniprot_name= tags.find("{http://www.hmdb.ca}uniprot_name")
        uniprot_name = uniprot_name.text
        data["subject"]["uniprot_name"]=uniprot_name

        genbank_protein_id= tags.find("{http://www.hmdb.ca}genbank_protein_id")
        data["subject"]["genbank_protein_id"]=genbank_protein_id.text

        hgnc_id= tags.find("{http://www.hmdb.ca}hgnc_id")
        data["subject"]["hgnc_id"]=hgnc_id.text

        genbank_gene_id=tags.find("{http://www.hmdb.ca}genbank_gene_id")
        data["subject"]["genbank_gene_id"]=genbank_gene_id.text

        gene_name = tags.find("{http://www.hmdb.ca}gene_name")
        data["subject"]["gene_name"]=gene_name.text

        return data;
 

    # --- Create a dictionary to hold our metabolite mapping values from the metabolite XML ---
    def make_metbolite_dict():
            # --- Load in the metabolites XML --- 
        xml_as_bytes = open(meta_xml, 'rb').read()
        metabolite_tree = etree_lxml.fromstring(xml_as_bytes)
        metabolite = metabolite_tree.findall('{http://www.hmdb.ca}metabolite', {})
        mapping_dict={}

        for meta in metabolite:
            accession=meta.find('{http://www.hmdb.ca}accession')
            kegg=meta.find('{http://www.hmdb.ca}kegg_id')
            chemspider=meta.find('{http://www.hmdb.ca}chemspider_id')
            chebi=meta.find('{http://www.hmdb.ca}chebi_id')
            pubchem=meta.find('{http://www.hmdb.ca}pubchem_compound_id')

            mapping_dict.setdefault(accession.text, {
                                                        "kegg_id":kegg.text,
                                                        "chemspider_id": chemspider.text,
                                                        "chebi_id": chebi.text,
                                                        "pubchem_compound_id": pubchem.text

        })

        return mapping_dict;



    # --- Enter mapping IDs into the document --- 
    def enter_mapping_ids(mapping_dict, text, data):
        # get the extra IDs from the metabolite xml
        # {'kegg_id': 'C01092', 'chemspider_id': '4578', 'chebi_id': '127029', 'pubchem_compound_id': '4740'}
        data["object"]["kegg_id"]=mapping_dict[text]["kegg_id"]
        data["object"]["chemspider_id"]=mapping_dict[text]["chemspider_id"]
        data["object"]["chebi_id"]=mapping_dict[text]["chebi_id"]
        data["object"]["pubchem_compound_id"]=mapping_dict[text]["pubchem_compound_id"]

        return mapping_dict;

    # ------------------------------------------------- 

    def construct_rec(protein_tree, mapping_dict):
        records=[]
        # --- Iterate over the root ---
        for tags in protein_tree[0].iter():
            metabolites = [] # holder for metabolite_associations.metabolite, the associations w/o references
            
            if tags.tag == "{http://www.hmdb.ca}protein": # get the protein nodes only
                
                # we need the first accession number, this is main protein _id in our doc
                _id = tags.find("{http://www.hmdb.ca}accession")
                _id = _id.text   
                protein_type=tags.find("{http://www.hmdb.ca}protein_type")
                ct=1 # setup counter for the associations
                
                
                # ---------- Metabolite associations with references ------------  
                for m in tags.findall("{http://www.hmdb.ca}metabolite_references"):
                    for ref in m:
                        # set the main accession id 
                        _id2=_id+"_%s"%ct
                        ct+=1 # update accession counter
                        # create dictionary document
                        data={}
                        data={ "_id": _id2, "pmid": None, "subject": { "protein_type": protein_type.text}, "object":{}}
            
                        # pull out the reference tags and get the pubmed_id
                        for met_ref in ref.findall("{http://www.hmdb.ca}reference"):
                            for refs in met_ref:
                                if "pubmed_id" in refs.tag:
                                    data['pmid']=refs.text

                        for met in ref.findall("{http://www.hmdb.ca}metabolite"):
                            for info in met:
                                tag=info.tag.split("}")[1]
                                text=info.text
                                data["object"][tag]=text

                                if "accession" in tag:
                                    enter_mapping_ids(mapping_dict, text, data)
                                
                        # Call enter_subject method to fill in subject data 
                        data=enter_subject(data,tags)  
                        #data=dict_convert(data,keyfn=process_key)
                        #data = dict_sweep(data,vals=[np.nan])
                        records.append(data)

                # ---------- Metabolite associations without references ------------         
                # find the metabolite_association tags and extract the information
                for met_assc in tags.findall("{http://www.hmdb.ca}metabolite_associations"):
                    for met_assc_ in met_assc.findall("{http://www.hmdb.ca}metabolite"):
                        for met_assc_id in met_assc_.findall("{http://www.hmdb.ca}accession"):
                            
                            # --- Check for duplicate ID, if found, skip making document --- 
                            # if the metabolite association was already present above (in metabolite_refereces)
                            # we want to pass adding id to dict to avoid making a duplicate document 
                            pass_assc=False # set bool 
                            for elem in data_list:
                                if met_assc_id.text == elem['object']['accession']:                            
                                    pass_assc = True
                                    
                            # if bool is True pass making duplicate doc       
                            if pass_assc==True: 
                                pass

                            # else create the document 
                            else:
                                # create data dict for association accession 
                                data={"_id": _id+"_%s"%ct, 'pmid': 'Unknown', 'subject':{}, 'object':{'accession': met_assc_id.text} }
                                ct+=1 # update the id counter 

                                enter_mapping_ids(mapping_dict, met_assc_id.text, data) # add the mapping ids for this accession

                                for met_assc_name in met_assc_.findall("{http://www.hmdb.ca}name"):
                                    data["object"]['name'] = met_assc_name.text

                                # Call enter_subject method to fill in subject data 
                                data=enter_subject(data,tags)  
                                #data=dict_convert(data,keyfn=process_key)
                                #data = dict_sweep(data,vals=[np.nan])
                                records.append(data)        

        return records;

    # --- Execute Program ---
    process_key = lambda k: k.replace(" ","_").lower()

    # --- Set input XML file path ---
    protein_xml = os.path.join(data_folder, "hmdb_proteins.xml")
    meta_xml = os.path.join(data_folder, "hmdb_metabolites.xml")

    # --- Load XML Data --- 
    xml_data = open(protein_xml, 'r', encoding='UTF-8').read()  # Read file
    protein_tree = ET.XML(xml_data)  # Parse XML
    data_list=[] # final data holder 
    mapping_dict=make_metbolite_dict() # load metabolite file and get the mapping ids 
    records=construct_rec(protein_tree, mapping_dict)
    #print(records)
    
    if(records):
        for record in records:
            yield record #print(json.dumps(record, sort_keys=False, indent=4))


    