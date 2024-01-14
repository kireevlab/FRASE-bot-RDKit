#!/usr/bin/env python3
# coding: utf-8

# In[2]:


import json, os, shutil
import numpy as np
import tensorflow as tf
from tensorflow import keras
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from tensorflow.keras.optimizers import SGD
from sklearn.metrics import confusion_matrix
import seaborn as sns
import matplotlib.pyplot as plt


# In[51]:


target_file = 'IF_target.json'
FRASEdb_file = 'IF_FRASEdb.json'

with open(target_file, 'r') as file1:
    if_list_target = json.load(file1)
print('The length from target is',len(if_list_target))

with open(FRASEdb_file, 'r') as file2:
    if_list_frasedb = json.load(file2)
print('The length from FRASE database is',len(if_list_frasedb))


# In[3]:


# Combine the frase and predicted data and create labels

# Create labels (0 for predicted data, 1 for frase data)
if_list_frasedb_labels = np.ones(len(if_list_frasedb))
if_list_target_labels = np.zeros(len(if_list_target))

# Combine the data into a single dataset
combined_data = np.concatenate((if_list_frasedb, if_list_target), axis=0)
combined_labels = np.concatenate((if_list_frasedb_labels, if_list_target_labels), axis=0)

print(combined_data.shape)
print(combined_labels.shape)

# Standardize the data using StandardScaler
scaler = StandardScaler()
combined_data = scaler.fit_transform(combined_data)


# Convert combined_data and combined_labels to NumPy arrays
combined_data = np.array(combined_data,dtype='float')
combined_labels = np.array(combined_labels)


#Hyperparameters
ts = 0.2
EPOCHS = 30
LearnRate = 0.001
drop = 0.3

# split the input dataset
(trainX, testX, trainY, testY) = train_test_split(combined_data,combined_labels,test_size=ts, random_state=10)

# Define the model
model = keras.Sequential([
    keras.layers.Dense(128, activation='relu', input_shape=(264,)),
    keras.layers.Dropout(drop),
    keras.layers.Dense(64, activation='relu'),
    keras.layers.Dropout(drop),
    keras.layers.Dense(32, activation='relu'),
    keras.layers.Dropout(drop),
    keras.layers.Dense(1, activation='sigmoid')  # Output layer with sigmoid activation
])



# Compile the model
model.compile(
    loss='binary_crossentropy',  # Binary cross-entropy for binary classification
    optimizer=SGD(learning_rate=LearnRate),
    metrics=['accuracy']
)

# Train the model
print('Start to train the model')
H = model.fit(trainX, trainY, validation_data=(testX,testY), epochs=EPOCHS)  # Adjust the number of epochs

# evaluate the model
print('Start to evaluate the model')
model.evaluate(trainX, trainY)

# after training dataset, draw the loss and validation plots
#N = np.arange(0, EPOCHS)
#plt.style.use("ggplot")
#plt.figure()
#plt.plot(N, H.history["loss"], label="train_loss")
#plt.plot(N, H.history["val_loss"], label="val_loss")
#plt.plot(N, H.history["accuracy"], label="train_acc")
#plt.plot(N, H.history["val_accuracy"], label="val_acc")
#plt.title("Training Loss and Accuracy for Neural Network")
#plt.xlabel("Epoch #")
#plt.ylabel("Loss/Accuracy")
#plt.legend()
#plt.savefig(f'NN-loss-Acc-lr{LearnRate}-epochs{EPOCHS}-testsize{ts}-dropout{drop}.png')


# ## predict class for target:

# In[36]:


# Make predictions on new data (optional)
predicted_class_target = model.predict(if_list_target)
#print(f"Predicted Class for validation set: {predicted_class_target}")

score_ranges = [(0.0, 0.5), (0.5, 0.6), (0.6, 0.7), (0.7, 0.8), (0.8, 0.9), (0.9, 1.0)]
score_lists = [[] for _ in range(len(score_ranges))]

for i in predicted_class_target:
    for idx, (lower, upper) in enumerate(score_ranges):
        if lower < i <= upper:
            score_lists[idx].append(i)
            
print("predict class for target:")
for idx, (lower, upper) in enumerate(score_ranges):
    print(f'score({lower}~{upper}): {len(score_lists[idx])}')

### last section
#################
threshold = 0.8 #
#################
scores = [score for index, score in enumerate(predicted_class_target) if score >= threshold]
indexs = [index+1 for index, score in enumerate(predicted_class_target) if score >= threshold]

with open('IF_target_sort.txt', 'r') as file3:
    lines = file3.readlines()
frase_index = []
frase_name = []
for line in lines:
    columns = line.split()
    frase_index.append(columns[0])
    frase_name.append(columns[1])
frase_com = list(zip(frase_index,frase_name))

extract = []
for i in indexs:
    for j in frase_com:
        int_j = int(j[0])
        if i == int_j:
            extract.append(j[1])

folder_name = "screen_ligFragments"
if not os.path.exists(folder_name):
    os.mkdir(folder_name)
    print(f"Folder '{folder_name}' sucessfully created.")
else:
    print(f"Folder '{folder_name}' already exists.")

targetML = "target_ML"

for file_name in extract:
    shutil.move(f"{targetML}/{file_name}", folder_name)

print("FRASE database screening is completed now. Cheers ... ")

# ## predict class for FRASE database:

# In[37]:


# Make predictions on new data (optional)
predicted_class_frasedb = model.predict(if_list_frasedb)
#print(f"Predicted Class for validation set: {predicted_class_frasedb}")

score_lists = [[] for _ in range(len(score_ranges))]

for i in predicted_class_frasedb:
    for idx, (lower, upper) in enumerate(score_ranges):
        if lower < i <= upper:
            score_lists[idx].append(i)

print("predict class for FRASE database:")
for idx, (lower, upper) in enumerate(score_ranges):
    print(f'score({lower}~{upper}): {len(score_lists[idx])}')


# ## predict class for Train set:

# In[38]:


# Make predictions on new data (optional)
predicted_class_trainX = model.predict(trainX)
#print(f"Predicted Class for training set: {predicted_class_trainX}")

score_lists = [[] for _ in range(len(score_ranges))]

for i in predicted_class_trainX:
    for idx, (lower, upper) in enumerate(score_ranges):
        if lower < i <= upper:
            score_lists[idx].append(i)

print("predict class for Train set:")
for idx, (lower, upper) in enumerate(score_ranges):
    print(f'score({lower}~{upper}): {len(score_lists[idx])}')


# ## predict class for Test set:

# In[39]:


# Make predictions on new data (optional)
predicted_class_testX = model.predict(testX)
#print(f"Predicted Class for testing set: {predicted_class_testX}")

score_lists = [[] for _ in range(len(score_ranges))]

for i in predicted_class_testX:
    for idx, (lower, upper) in enumerate(score_ranges):
        if lower < i <= upper:
            score_lists[idx].append(i)

print("predict class for Test set:")
for idx, (lower, upper) in enumerate(score_ranges):
    print(f'score({lower}~{upper}): {len(score_lists[idx])}')


# ## Confusion matrix:

# In[47]:


# Make predictions on new data
predicted_class_combine = model.predict(combined_data)
#print(f"Predicted Class for testing set: {predicted_class_testX}")

score_lists = [[] for _ in range(len(score_ranges))]

for i in predicted_class_combine:
    for idx, (lower, upper) in enumerate(score_ranges):
        if lower < i <= upper:
            score_lists[idx].append(i)

for idx, (lower, upper) in enumerate(score_ranges):
    print(f'score({lower}~{upper}): {len(score_lists[idx])}')
    
# Assuming you have already trained your model and made predictions
predicted_labels = predicted_class_combine # Replace with your test data
true_labels = combined_labels  # Replace with your true labels

# Calculate the confusion matrix
conf_matrix = confusion_matrix(true_labels, (predicted_labels > threshold).astype(int))  # Adjust the threshold as needed
reversed_conf_matrix = [list(reversed(row)) for row in reversed(conf_matrix)]
print(reversed_conf_matrix)

# Create a heatmap for the confusion matrix
sns.set(font_scale=1.2)
plt.figure(figsize=(8, 6))
ax = sns.heatmap(reversed_conf_matrix, annot=True, fmt='d', cmap='Blues', cbar=True, square=True, linewidths=0.5)
# Set the x-labels at the top
ax.xaxis.tick_top()
ax.set_xticklabels(reversed(ax.get_xticklabels()))
ax.set_yticklabels(reversed(ax.get_yticklabels()))
# Add labels and title
plt.ylabel('Actual')
plt.title('Predict')

# save the plot
plt.savefig('confusion_matrix.png')
