from uniprot_wrapper import UniprotWrapper
import json

filename = "SynGO_export_2018-03-17.json"
with open(filename) as f:
    data = json.load(f)

wrapper = UniprotWrapper()

id_map = {}

for a in data["SynGO"]:
    for m in a['models']:
        uniprot_id = m["uniprot"]
        if uniprot_id not in id_map:
            id_map[uniprot_id] = ""

counter = 0
with open("results.txt", 'w') as wf:
    for k in id_map:
        id_map[k] = wrapper.one_off_call(k)
        print(k + ": " + id_map[k])
        wf.write(k + ": " + id_map[k] + "\n")
        counter += 1
        # if counter > 10:
        #     break
