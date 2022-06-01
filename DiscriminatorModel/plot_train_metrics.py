import seaborn as sns
import matplotlib.pyplot as plt
from custom_cpu_unpickler import CPU_Unpickler
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import pickle
import numpy as np
import pandas as pd

contents = CPU_Unpickler(open("DiscriminatorModel/data/scores/esm_discriminator_model_0_30_testing.pkl", "rb")).load()

loss, accuracy = contents

fig, axs = plt.subplots(1, 2)
sns.lineplot(x=[i for i in range(len(loss))], y=loss, ax=axs[0])
axs[0].set_xlabel("Checkpoint")
axs[0].set_ylabel("Loss")

sns.lineplot(x=[i for i in range(len(accuracy))], y=accuracy, ax=axs[1], color="orange")
axs[1].set_xlabel("Checkpoint")
axs[1].set_ylabel("Accuracy")
plt.show()


# sns.lineplot(x=[i+1 for i in range(len(train_accuracies))], y=train_accuracies, label="Train")
# sns.lineplot(x=[i+1 for i in range(len(val_accuracies))], y=val_accuracies, label="Validate")
# plt.xlabel("Epoch")
# plt.ylabel("Accuracies")
# plt.show()

# y_pred, y_true = pickle.load(open("DiscriminatorModel/data/scores/esm_discriminator_model_0_21_labels.pkl", "rb"))
# matrix = confusion_matrix(y_true, y_pred)
# sns.heatmap(matrix, annot=True, fmt=".1f", xticklabels=["Natural", "Generated"], yticklabels=["Natural", "Generated"])
# plt.xlabel("Predicted")
# plt.ylabel("Actual")
# plt.show()
