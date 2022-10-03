import os
from lxml import etree

"""
The cardinality relationships between the elements of our interest in hmdb_proteins.xml are described below:

    <protein>
        <metabolite_references>                    <!-- ONLY ONE for each <protein> -->
            <metabolite_reference>                 <!-- MULTIPLE for each <metabolite_references> -->
                <metabolite>                       <!-- ONLY ONE for each <metabolite_reference> -->
                    <name> XXX </name>             <!-- ONLY ONE for each <metabolite> -->
                    <accession> XXX </accession>   <!-- ONLY ONE for each <metabolite> -->
                </metabolite>
                <reference>                        <!-- ONLY ONE for each <metabolite_reference> -->
                    <pubmed_id> XXX </pubmed_id>   <!-- ONLY ONE for each <reference> -->
                </reference>
            </metabolite_reference>

            <metabolite_reference>
                <!-- ignored -->
            </metabolite_reference>
        </metabolite_references>

        <metabolite_associations>                  <!-- ONLY ONE for each <protein> -->
            <metabolite>                           <!-- MULTIPLE for each <metabolite_associations> -->
                <name> XXX </name>                 <!-- ONLY ONE for each <metabolite> -->
                <accession> XXX </accession>       <!-- ONLY ONE for each <metabolite> -->
            </metabolite>

            <metabolite>
                <!-- ignored -->
            </metabolite>
        </metabolite_associations>
    </protein>

For those 1-to-1 relationships, lxml `find()` method will suffice to locate the element;
for those 1-to-N relationships, lxml `find_all()` method should be used.
"""

# Seen accession IDs
# When creating documents from <metabolite_associations>, we need to check if the accession IDs found inside are not used yet
OBJECT_ACCESSIONS = set()


def make_metabolite_dict(metabolites_xml_path: str):
    """
    # Create a dictionary to hold our metabolite mapping values
    # from the metabolite XML
    """
    metabolite_dict = {}  # initialize dictionary

    # --- Load in the metabolites XML ---
    root = etree.parse(metabolites_xml_path)
    for metabolite in root.findall('{http://www.hmdb.ca}metabolite'):
        # all 'accession.text's are not None
        # other text fields may contain None
        accession = metabolite.find('{http://www.hmdb.ca}accession')
        kegg = metabolite.find('{http://www.hmdb.ca}kegg_id')
        chemspider = metabolite.find('{http://www.hmdb.ca}chemspider_id')
        chebi = metabolite.find('{http://www.hmdb.ca}chebi_id')
        pubchem = metabolite.find('{http://www.hmdb.ca}pubchem_compound_id')

        metabolite_dict[accession.text] = {
            "kegg_id": kegg.text or "",
            "chemspider_id": chemspider.text or "",
            "chebi_id": chebi.text or "",
            "pubchem_compound_id": pubchem.text or ""
        }

    return metabolite_dict


def enter_subject(data: dict, protein: etree._Element):
    """
    Fill subject fields into document
    """
    # protein_type, uniprot_id, uniprot_name, genbank_protein_id, hgnc_id, genbank_gene_id, and gene_name.
    protein_type = protein.find("{http://www.hmdb.ca}protein_type")
    data["subject"]["protein_type"] = protein_type.text or ""

    uniprot_id = protein.find("{http://www.hmdb.ca}uniprot_id")
    data["subject"]["uniprot_id"] = uniprot_id.text or ""

    uniprot_name = protein.find("{http://www.hmdb.ca}uniprot_name")
    data["subject"]["uniprot_name"] = uniprot_name.text or ""

    genbank_protein_id = protein.find("{http://www.hmdb.ca}genbank_protein_id")
    data["subject"]["genbank_protein_id"] = genbank_protein_id.text or ""

    hgnc_id = protein.find("{http://www.hmdb.ca}hgnc_id")
    data["subject"]["hgnc_id"] = hgnc_id.text or ""

    genbank_gene_id = protein.find("{http://www.hmdb.ca}genbank_gene_id")
    data["subject"]["genbank_gene_id"] = genbank_gene_id.text or ""

    gene_name = protein.find("{http://www.hmdb.ca}gene_name")
    data["subject"]["gene_name"] = gene_name.text or ""

    return data


def enter_object(data: dict, metabolite_dict: dict, accession: str):
    """
    Fill object fields into document
    """
    # get the extra IDs from the metabolite xml
    # {'kegg_id': 'C01092', 'chemspider_id': '4578', 'chebi_id': '127029', 'pubchem_compound_id': '4740'}

    data["object"]["kegg_id"] = metabolite_dict[accession]["kegg_id"]
    data["object"]["chemspider_id"] = metabolite_dict[accession]["chemspider_id"]
    data["object"]["chebi_id"] = metabolite_dict[accession]["chebi_id"]
    data["object"]["pubchem_compound_id"] = metabolite_dict[accession]["pubchem_compound_id"]

    return data


def construct_rec(protein: etree._Element, metabolite_dict: dict):
    """
    # Construct Record
    # This is essentially the main controller method, it contains the major loops and data entry
    """
    protein_accession = protein.find("{http://www.hmdb.ca}accession")  # main accession id
    ct = 1  # setup counter for the associations

    # ---------- Metabolite associations with references ------------
    # every protein element has only one "metabolite_references", which contains multiple "metabolite_reference"s
    metabolite_references = protein.find("{http://www.hmdb.ca}metabolite_references")
    for metabolite_reference in metabolite_references.findall("{http://www.hmdb.ca}metabolite_reference"):
        metabolite = metabolite_reference.find("{http://www.hmdb.ca}metabolite")
        metabolite_accession = metabolite.find("{http://www.hmdb.ca}accession")
        metabolite_name = metabolite.find("{http://www.hmdb.ca}name")

        # make sure that `enter_object` won't raise a KeyError
        if metabolite_accession.text not in metabolite_dict:
            continue

        reference = metabolite_reference.find("{http://www.hmdb.ca}reference")
        reference_pubmed_id = reference.find("{http://www.hmdb.ca}pubmed_id")

        _id = protein_accession.text + "_%s" % ct  # set the main accession id
        ct += 1  # update accession counter

        data = {
            "_id": _id,
            "pmid": None,
            "subject": {},
            "object": {}
        }

        data['pmid'] = reference_pubmed_id.text
        data["object"]["accession"] = metabolite_accession.text
        data["object"]["name"] = metabolite_name.text

        data = enter_object(data, metabolite_dict, metabolite_accession.text)
        data = enter_subject(data, protein)

        OBJECT_ACCESSIONS.add(metabolite_accession.text)
        yield data

    # ---------- Metabolite associations without references ------------
    # find the metabolite_association tags and extract the information
    metabolite_associations = protein.find("{http://www.hmdb.ca}metabolite_associations")
    for metabolite in metabolite_associations.findall("{http://www.hmdb.ca}metabolite"):
        metabolite_accession = metabolite.find("{http://www.hmdb.ca}accession")
        metabolite_name = metabolite.find("{http://www.hmdb.ca}name")

        # --- Check for duplicate ID, if found, skip making document ---
        # if the metabolite association was already present above (in metabolite_refereces)
        # we want to pass adding id to dict to avoid making a duplicate document
        if metabolite_accession.text in OBJECT_ACCESSIONS:
            continue

        # make sure that `enter_object` won't raise a KeyError
        if metabolite_accession.text not in metabolite_dict:
            continue

        _id = protein_accession.text + "_%s" % ct  # set the main accession id
        ct += 1  # update accession counter

        data = {
            "_id": _id,
            'pmid': 'Unknown',
            'subject': {},
            'object': {}
        }

        data["object"]["accession"] = metabolite_accession.text
        data["object"]['name'] = metabolite_name.text

        data = enter_object(data, metabolite_dict, metabolite_accession.text)
        data = enter_subject(data, protein)

        OBJECT_ACCESSIONS.add(metabolite_accession.text)
        yield data


def load_hmdb_data(data_folder):
    """
    # HMDB Association Data Load Method Main Method
    """
    # --- Set input XML file path ---
    protein_xml = os.path.join(data_folder, "hmdb_proteins.xml")
    meta_xml = os.path.join(data_folder, "hmdb_metabolites.xml")

    # load metabolite XML file and get the mapping ids
    metabolite_dict = make_metabolite_dict(meta_xml)

    # Parse protein XML file
    protein_tree = etree.parse(protein_xml)
    for protein in protein_tree.findall("{http://www.hmdb.ca}protein"):
        yield from construct_rec(protein, metabolite_dict)


# if __name__ == "__main__":
#     for doc in load_hmdb_data("."):
#         print(doc)
