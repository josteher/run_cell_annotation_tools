from sklearn.metrics import classification_report
import pandas as pd

preds = pd.read_csv("prediction_summary.csv")

for i in preds.columns:
	pd.DataFrame(classification_report(preds.true_test_label, preds[i],output_dict = True)).T.to_csv("/data/gruen/herman/scDatasets/run_CellAnnotation_tools/run_Benchmarking_TM/classification_reports/" + i + ".csv")

