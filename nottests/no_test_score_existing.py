from goreverselookup import ReverseLookup, fisher_exact_test

SOIs = {}

# data_file = "results/data.json"
# data_file = "C:\\Aljosa\\Development\\Environments\\goreverselookup_test\\test-two-tailed\\two_tailed_test-chronic_infl_cancer\\twotail\\results\\data.json"
data_file = "C:\\Aljosa\\Development\\Environments\\goreverselookup_test\\test-two-tailed\\two_tailed_test-rhart\\twotail-results\\results\\data.json"
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

# change any of the model settings
# model.model_settings.indirect_annotations_max_depth = 5
# model.model_settings.include_indirect_annotations = True
# model.model_settings.two_tailed = True
# model.model_settings.exclude_opposite_regulation_direction_check = True

fisher_score = fisher_exact_test(model)
# model.score_products(score_classes=[fisher_score]) # re-enable this if you want to re-score products

model.perform_statistical_analysis(
	test_name="fisher_test", 
	filepath="results/results.json",
    two_tailed=model.model_settings.two_tailed,
    use_dest_dir=True
)

num_relevant_genes = 0
for SOIs_label,relevant_genes in model.statistically_relevant_products.items():
    num_relevant_genes += len(relevant_genes)
    print(f"{SOIs_label} : {len(relevant_genes)} genes")
print(f"total relevant genes: {num_relevant_genes}")

#print(f"Number of statistically significant genes: {len(model.statistically_relevant_products['chronic_inflammation+:cancer+'])}")
# model.save_model("results/data.json")