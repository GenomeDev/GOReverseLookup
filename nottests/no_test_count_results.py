# Counts statistically relevant GO terms and displays them.

from goreverselookup import JsonUtil

# results_filepath = "research_models\\test_models\\two_tailed_test-rhart\\twotail-results\\results.json"
results_filepath = "C:\\Aljosa\\Development\\Environments\\goreverselookup_test\\test-two-tailed\\two_tailed_test-chronic_infl_cancer\\twotail\\results\\results.json"
json = JsonUtil.load_json(results_filepath)

num_relevant_genes = 0
for SOIs_label,relevant_genes in json.items():
    num_relevant_genes += len(relevant_genes)
    print(f"{SOIs_label} : {len(relevant_genes)} genes")
print(f"total relevant genes: {num_relevant_genes}")