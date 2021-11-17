# HMDB Protein Metabolite Associations Data Plugin
**Data from Human Metbolite Database [HMDB](https://hmdb.ca/downloads)**
  
## Notes  

- The input metabolite file, `hmdb_metabolites.xml`, the current version 5.0 is corrupt and produces error `BAD CRC-32`, reference file [here]() for more details. Version 4.0 of the file is being used, *may want to update this in the future.*  

## <u> Document Structure </u>

```
[
    {
        "_id": "HMDBP00001_1",
        "pmid": "11752352",
        "subject": {
            "protein_type": "Unknown",
            "uniprot_id": "P21589",
            "uniprot_name": "5NTD_HUMAN",
            "genbank_protein_id": "23897",
            "hgnc_id": "HGNC:8021",
            "genbank_gene_id": "X55740",
            "gene_name": "NT5E"
        },
        "object": {
            "name": "Pentoxifylline",
            "accession": "HMDB0014944",
            "kegg_id": "C07424",
            "chemspider_id": "4578",
            "chebi_id": "127029",
            "pubchem_compound_id": "4740"
        }
    },
    {
        "_id": "HMDBP00001_2",
        "pmid": "16426349",
        "subject": {
            "protein_type": "Unknown",
            "uniprot_id": "P21589",
            "uniprot_name": "5NTD_HUMAN",
            "genbank_protein_id": "23897",
            "hgnc_id": "HGNC:8021",
            "genbank_gene_id": "X55740",
            "gene_name": "NT5E"
        },
        "object": {
            "name": "Pentoxifylline",
            "accession": "HMDB0014944",
            "kegg_id": "C07424",
            "chemspider_id": "4578",
            "chebi_id": "127029",
            "pubchem_compound_id": "4740"
        }
    },
    .
    .
    .
]
```