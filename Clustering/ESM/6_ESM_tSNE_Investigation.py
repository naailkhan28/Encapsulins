import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#Load t-SNE clustered dataframe
all_encapsulins = pd.read_csv("Clustering Analysis/all_encapsulins_ESM_tSNE.csv")

#Plot t-SNE embeddings
sns.scatterplot(all_encapsulins["tSNE-x"], all_encapsulins["tSNE-y"], hue=all_encapsulins["Family"])
plt.show()