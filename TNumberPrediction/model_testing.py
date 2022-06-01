import pandas as pd
import pickle
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.preprocessing import LabelBinarizer
import seaborn as sns
import matplotlib.pyplot as plt

#Load dataframe and sequence embeddings
encapsulin_data = pd.read_csv("TNumberPrediction/data/all_seqs.csv")
X = pickle.load(open("TNumberPrediction/data/encapsulin_t_number_embeddings.pkl", "rb"))

# #Preprocess data and split into training and test sets
# y = pd.DataFrame(encapsulin_data["T"])
# preprocessor = LabelBinarizer()
# y = preprocessor.fit_transform(y)


# X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=14)

# #Set up model and hyperparameter search
# clf = RandomForestClassifier()
# param_search = {
#     "n_estimators": [50, 100, 250, 500],
#     "max_features": [0.8, 0.5, 0.3],
# }

# #Perform grid search for hyperparameters with cross-validation
# output = GridSearchCV( clf, param_search, cv=10, return_train_score=False)
# output.fit(X_train, y_train)

# #Show best parameters and scores
# print(output.best_params_)
# print(output.best_score_)

# #Evaluate performance on the test set
# print(output.score(X_test, y_test))