"""
Displays the GO terms linked to specific SOIs.
"""

from goreverselookup import ReverseLookup

# data_file = "C:\\Aljosa\\Development\\Environments\\goreverselookup_test\\test-two-tailed\\two_tailed_test-rhart\\twotail-results\\results\\data.json"
# data_file = "C:\\Aljosa\\Development\\GOReverseLookup-Env\\research_models\\chr_infl_cancer-v6_FINAL\\IEA-\\ind_ann,p=0.05,IEA- (149)\\results\\data.json"
# data_file = "C:\\Aljosa\\Development\\GOReverseLookup\\research_models\\chronic-inflammation_cancer\\IEA-\\ind_ann,p=0.05,IEA- (149)\\input.txt"
data_file = "F:\\Aljoša\\Dropbox\\Raziskovalni članki\\GOReverseLookup\\article - Computers in biology and medicine\\CBM_submission\\Supplemental_Table_S3.txt"

model = ReverseLookup.from_input_file(data_file)
#model = ReverseLookup.load_model(data_file)

SOIs = {}

# find which GO terms belong to the SOIs
for goterm in model.goterms:
    soi_name = goterm.SOIs[0].get("SOI") # eg chronic_inflammation
    soi_direction = goterm.SOIs[0].get("direction") # eg '-'
    soi = f"{soi_name}{soi_direction}" # eg chronic_inflammation-
    if soi not in SOIs:
        SOIs[soi] = [goterm]
    else:
        SOIs[soi].append(goterm)
        
for SOI,goterms in SOIs.items():
    print(f"SOI '{SOI}' is associated with {len(goterms)} GO terms:")
    i = 0
    for goterm in goterms:
        i += 1
        print(f"  - [{i}]: {goterm.id}\t{goterm.name}")