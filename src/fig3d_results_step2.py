#%%
import pandas as pd
from pycaret.classification import *
from imblearn.over_sampling import *
from sklearn.model_selection import RepeatedKFold
import matplotlib.pyplot as plt
#%%
n_jobs = 60
X = pd.read_csv("X.csv")
y = pd.read_csv("y.csv")
dataset = pd.concat([X,y], axis=1)
feature_table = pd.read_csv("features_table.csv")

# %%
compare_df = pd.DataFrame()
max_mcc = 0
for i in sorted(list(set(feature_table["Freq"])),reverse=True):
    if i < 50:
        continue
    features = list(feature_table[feature_table["Freq"] >=i]["Feature"]) + ["group"]
    if len(features) <= 3:
        continue
    rkf = RepeatedKFold(n_splits=10, n_repeats=5, random_state=2021)
    adasyn = ADASYN(sampling_strategy='minority',random_state=2021)
    clf = setup(dataset[features], target = 'group', session_id=2021,log_experiment=True, experiment_name=str(i),silent=True,fix_imbalance = True,fix_imbalance_method=adasyn,fold_strategy=rkf,n_jobs=n_jobs)
    model= create_model('et')
    tuned_model = tune_model(model,n_iter = 200,search_library="scikit-optimize",optimize = "MCC")
    re = pull()
    mean_cv = pd.DataFrame(re.loc["Mean"])
    mean_cv.columns= [i]
    compare_df = mean_cv if compare_df.empty else pd.concat([compare_df, mean_cv],axis=1)
    if mean_cv[i]["MCC"] > max_mcc :
        max_mcc = mean_cv[i]["MCC"]
# %%
compare_df.to_csv("model_cv.csv")