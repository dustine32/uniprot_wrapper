import json
import requests
from pandas import read_csv
from io import StringIO
from os import path

def make_uniprot_call(uniprot_ids, current_map=None):
    request_min = 0
    while request_min < len(uniprot_ids):
        for field in dbs_to_lookup + other_fields_to_fetch:
            r = requests.get('http://www.uniprot.org/uploadlists/?from=ACC&to=' + field + '&format=tab&query=' + " ".join(uniprot_ids[request_min:request_min+500]))
            uniprot_results = read_csv(StringIO(r.text), delimiter='\t')

            for index, row in uniprot_results.iterrows():
                print(row[0] + " - " + field)
                current_map[row[0]][field] = row[1]

        request_min += 500 # For some reason requesting >1000 results in 400 error
    return current_map

def lookup_uniprot(uniprot_ids, current_map=None, isoform_check=True):
    if current_map is None:
        current_map = {}
        for uid in uniprot_ids:
            current_map[uid] = {}
    current_map = make_uniprot_call(uniprot_ids, current_map)

    # Adjust for isoforms
    if isoform_check:
        redo_ids = []
        for k in current_map:
            if current_map[k] == {}:
                redo_id = k.split("-")[0]
                redo_ids.append(redo_id)
        for rid in redo_ids:
            current_map[rid] = {}
        current_map = make_uniprot_call(redo_ids, current_map)

    return current_map

def get_field_for_id(current_map, field, uniprot_id):
    if field in current_map[uniprot_id]: 
        return str(current_map[uniprot_id][field])

def get_noctua_gene_id(current_map, dbs, uniprot_id):
    noctua_gene_id = None
    for db in dbs:
        if db in current_map[uniprot_id]:
            noctua_gene_id = get_field_for_id(current_map, db, uniprot_id)
    return noctua_gene_id

filename = "SynGO_export_2018-03-17.json"
with open(filename) as f:
    data = json.load(f)

extensions = []
ext_genes = []
ext_sets = []
id_map = {}

for a in data["SynGO"]:
    for m in a['models']:
        uniprot_id = m["uniprot"]
        if uniprot_id not in id_map:
            id_map[uniprot_id] = {}

dbs_to_lookup = ["MGI_ID", "RGD_ID", "FLYBASE_ID", "WORMBASE_ID", "HGNC_ID"]
other_fields_to_fetch = ["GENENAME"]
id_map = lookup_uniprot(list(id_map.keys()))

# <http://rgd.mcw.edu/rgdweb/report/gene/main.html?id=620107>
# <http://www.informatics.jax.org/accession/MGI:MGI:95615>

for a in data["SynGO"]:
    models = []
    for m in a['models']:
        uniprot_id = m['uniprot']
        if id_map[uniprot_id] == {}:
            uniprot_id = uniprot_id.split("-")[0] # Adjust for isoforms
        noctua_gene_id = get_noctua_gene_id(id_map, dbs_to_lookup, uniprot_id)
        if noctua_gene_id is None:
            print(m['uniprot'] + " - " + a['combi_id'] + " - " + str(get_field_for_id(id_map, "GENENAME", uniprot_id)))
            noctua_gene_id = m["uniprot"]
        m["noctua_gene_id"] = noctua_gene_id
        for field in other_fields_to_fetch:
            if field in id_map[uniprot_id]:
                m[field] = get_field_for_id(id_map, field, uniprot_id)
        ext_relations = []
        for e in m['extensions']:
            for k in e.keys():
                go_terms = []
                uberon_terms = []
                other_terms = []
                if k not in extensions:
                    extensions.append(k)
                if k not in ext_relations:
                    ext_relations.append(k)
                for t in e[k]:
                    if ":" not in t and t not in ext_genes:
                        ext_genes.append(t)
                    if t.startswith("GO:"):
                        go_terms.append(t)
                    elif t.startswith("UBERON:"):
                        uberon_terms.append(t)
                    else:
                        other_terms.append(t)
                # if len(go_terms) > 1 or len(uberon_terms) > 1 or len(other_terms) > 0:
                #     ext_sets.append([a["combi_id"],e[k]])
        if len(set(ext_relations)) > 1:
            ext_sets.append(a["combi_id"])
        models.append(m)
    a['models'] = models

with open(path.splitext(filename)[0] + "_updated" + path.splitext(filename)[1], "w") as wf:
    json.dump(data, wf, indent=4)