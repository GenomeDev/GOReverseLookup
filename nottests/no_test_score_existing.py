from goreverselookup import ReverseLookup, fisher_exact_test

data_file = "results/data.json"
model = ReverseLookup.load_model(data_file)

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