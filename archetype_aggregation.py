import pandas as pd

""" import data """
# results

NaiveBayes_results = pd.read_csv("")
NaiveBayes_shap_values = pd.read_csv("")

NearestCentroid_results = pd.read_csv("results/nearest_centroid_results.csv")
NearestCentroid_shap_values = pd.read_csv("results/nearest_centroid_shap_values.csv")
KNeighborsClassifier_results = pd.read_csv("results/KNeighborsClassifier_results.csv")
KNeighborsClassifier_shap_values = pd.read_csv("results/KNeighborsClassifier_shap_values.csv")
RandomForestClassifier_results = pd.read_csv("results/RandomForestClassifier_results.csv")
RandomForestClassifier_shap_values = pd.read_csv("results/RandomForestClassifier_shap_values.csv")
ExtraTreesClassifier_results = pd.read_csv("results/ExtraTreesClassifier_results.csv")
ExtraTreesClassifier_shap_values = pd.read_csv("results/ExtraTreesClassifier_shap_values.csv")
DecisionTreeClassifier_results = pd.read_csv("results/decisionTreeClassifier_results.csv")
DecisionTreeClassifier_shap_values = pd.read_csv("results/DecisionTreeClassifier_shap_values.csv")
LGBMClassifier_results = pd.read_csv("results/LGBMClassifier_results.csv")
LGBMClassifier_shap_values = pd.read_csv("results/LGBMClassifier_shap_values.csv")
XGBClassifier_results = pd.read_csv("results/XGBClassifier_results.csv")
XGBClassifier_shap_values = pd.read_csv("results/XGBClassifier_shap_values.csv")

# shap values


""" Aggregate archetype results """
# Aggregating NearestNeighbour results
NearestNeighbour_results = pd.concat([NearestCentroid_results,
                                      KNeighborsClassifier_results,
                                      RandomForestClassifier_results])
NearestNeighbour_shap_values = pd.concat([NearestCentroid_shap_values,
                                          KNeighborsClassifier_shap_values,
                                          RandomForestClassifier_shap_values])

# Aggregating AveragingMethods results
AveragingMethods_results = pd.concat([ExtraTreesClassifier_results,
                                      DecisionTreeClassifier_results])
AveragingMethods_shap_values = pd.concat([ExtraTreesClassifier_shap_values,
                                          DecisionTreeClassifier_shap_values])

# NaiveBayes
NaiveBayes_results
NaiveBayes_shap_values

# Gradient Boosting results
GradientBoosting_results = pd.concat([XGBClassifier_results,
                                      LGBMClassifier_results
                                      ])
GradientBoosting_shap_values = pd.concat([XGBClassifier_shap_values,
                                          LGBMClassifier_shap_values
                                          ])

""" Save output """
NearestNeighbour_results.to_csv("NearestNeighbour_results.csv", index=False)
NearestNeighbour_shap_values.to_csv("NearestNeighbour_shap_values.csv", index=False)

AveragingMethods_results.to_csv("results/AveragingMethods_results.csv", index=False)
AveragingMethods_shap_values.to_csv("results/AveragingMethods_shap_values.csv", index=False)

GradientBoosting_results.to_csv("results/GradientBoosting_results.csv", index=False)
GradientBoosting_shap_values.to_csv("results/GradientBoosting_shap_values.csv", index=False)