#Imports
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.preprocessing import QuantileTransformer, RobustScaler
from sklearn.compose import ColumnTransformer

#Import encapsulin data
all_encapsulins = pd.read_csv("Clustering Analysis/all_encapsulins_processed.csv")
all_encapsulins = all_encapsulins.drop(["Unnamed: 0"], axis=1)

#Choose features

all_features = ['Length', 'Molecular Weight', 'pI', 'A_frequency',
       'C_frequency', 'D_frequency', 'E_frequency', 'F_frequency',
       'G_frequency', 'H_frequency', 'I_frequency', 'K_frequency',
       'L_frequency', 'M_frequency', 'N_frequency', 'P_frequency',
       'Q_frequency', 'R_frequency', 'S_frequency', 'T_frequency',
       'V_frequency', 'W_frequency', 'Y_frequency']

quantile_features = ['Length', 'Molecular Weight', 'pI', 'A_frequency',
       'D_frequency', 'E_frequency', 'F_frequency',
       'G_frequency', 'H_frequency', 'I_frequency', 'K_frequency',
       'L_frequency', 'M_frequency', 'N_frequency', 'P_frequency',
       'Q_frequency', 'R_frequency', 'S_frequency', 'T_frequency',
       'V_frequency', 'W_frequency', 'Y_frequency']

log_features = ['C_frequency']

#Preprocess data
preprocessor = ColumnTransformer(transformers=[

        ("quantile", QuantileTransformer(output_distribution="normal"), quantile_features),
        ("log", RobustScaler(), log_features)

], verbose=True, verbose_feature_names_out=True)

all_encapsulins_processed = pd.DataFrame(preprocessor.fit_transform(all_encapsulins), columns=preprocessor.get_feature_names_out())

fig, axs = plt.subplots(5, 5)

i = 0

for row in axs:
    for col in row:
        try:
            sns.histplot(all_encapsulins_processed.iloc[:, i], ax=col)
            i += 1
        except IndexError:
            break
        

plt.show()