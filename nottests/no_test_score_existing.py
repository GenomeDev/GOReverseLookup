from goreverselookup import ReverseLookup, fisher_exact_test

SOIs = {}

data_file = "results/data.json"
model = ReverseLookup.load_model(data_file)

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

model.model_settings.indirect_annotations_max_depth = 5
model.model_settings.include_indirect_annotations = True

fisher_score = fisher_exact_test(model)
model.score_products(score_classes=[fisher_score])


model.perform_statistical_analysis(
	test_name="fisher_test", 
	filepath="results/results.json", 
)
print(f"Number of statistically significant genes: {len(model.statistically_relevant_products['chronic_inflammation+:cancer+'])}")
model.save_model("results/data.json")