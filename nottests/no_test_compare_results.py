# Determines the Jaccard similarity index between two different models.

from goreverselookup import ReverseLookup

if1 = "C:\\Aljosa\Development\\GOReverseLookup-Env\\research_models\\chr_infl_cancer-v5\\IEA-\\ind_ann,p=0.05,IEA- (184)\\results\\data.json"
if2 = "C:\\Aljosa\\Development\\GOReverseLookup-Env\\research_models\\chr_infl_cancer-v5\\IEA+\\ind_ann,p=0.05,IEA+ (184)\\results\\data.json"

if3 = "C:\\Aljosa\\Development\\GOReverseLookup-Env\\research_models\\chr_infl_cancer-v5\\IEA-\\ind_ann,p=5e-8,IEA- (164)\\results\\data.json"
if4 = "C:\\Aljosa\\Development\\GOReverseLookup-Env\\research_models\\chr_infl_cancer-v5\\IEA+\\ind_ann,p=5e-8,IEA+ (163)\\results\\data.json"

if5 = "C:\\Aljosa\\Development\\GOReverseLookup-Env\\research_models\\chr_infl_cancer-v5\\IEA-\\no_ind_ann,p=0.05_IEA- (28)\\results\\data.json"
if6 = "C:\\Aljosa\\Development\\GOReverseLookup-Env\\research_models\\chr_infl_cancer-v5\\IEA+\\no_ind_ann,p=0.05,IEA+ (23)\\results\\data.json"

model_IEAneg = ReverseLookup.load_model(if5)
model_IEAplus = ReverseLookup.load_model(if6)

stat_rev_key = "chronic_inflammation+:cancer+"
results_IEAneg = model_IEAneg.statistically_relevant_products[stat_rev_key]
results_IEAplus = model_IEAplus.statistically_relevant_products[stat_rev_key]

results_IEAneg_set = set()
results_IEAplus_set = set()

for result_IEAneg in results_IEAneg:
    results_IEAneg_set.add(result_IEAneg.get("genename"))

for result_IEAplus in results_IEAplus:
    results_IEAplus_set.add(result_IEAplus.get("genename"))

intersection = results_IEAneg_set.intersection(results_IEAplus_set)
union = results_IEAneg_set.union(results_IEAplus_set)

jaccard_index = float(len(intersection) / len(union))

print(f"results IEA negative: {len(results_IEAneg_set)}")
print(f"results IEA positive: {len(results_IEAplus_set)}")
print(f"Jaccard similarity index: {len(intersection)} / {len(union)} = {jaccard_index}")
