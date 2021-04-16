#################
# Preprocessing #
#################
def balanceGrp(X, y, Xids, q=0.5):
    cnts = np.unique(y, return_counts=True)
    target_n = int(np.quantile(cnts[1], q))
    min_n = int(target_n/4)
    frames = np.zeros((len(cnts[0])*target_n, 1))
    for ctype in cnts[0]:
        fill_idx = int(np.argwhere(ctype==cnts[0])) * target_n
        ctype_idx = np.argwhere(y==ctype)
        np.random.seed(seed=1234)
        
        if ctype_idx.shape[0] < min_n:
            print(str(ctype) + " [Remove]: " + str(ctype_idx.shape[0])
                + " < " + str(min_n))
            sample_idx = np.zeros(target_n) - 1
        elif ctype_idx.shape[0] < target_n:
            print(str(ctype) + " [Upsample]: " + str(ctype_idx.shape[0])
                + " < " + str(target_n))
            sample_idx = np.random.choice(ctype_idx[:,0],
                                          size=target_n, replace=True)
        elif ctype_idx.shape[0] > target_n:
            print(str(ctype) + " [Downsample]: " + str(ctype_idx.shape[0])
                + " > " + str(target_n))
            sample_idx = np.random.choice(a=ctype_idx[:,0], size=target_n,
                                          replace=False)
        else:
            sample_idx = ctype_idx[:,0]
        frames[fill_idx:fill_idx+target_n,0] = sample_idx
    
    
    frames = frames[frames[:,0] >= 0,:]
    Xids_bal=np.asarray(Xids)[frames[:,0].astype('int'),:]
    X_bal=X[frames[:,0].astype('int'),:,:,:]
    y_bal=y[frames[:,0].astype('int'),:]
    return X_bal,y_bal,Xids_bal

##############
# CNN Models #
##############
class CNN:
    def __init__(self, y=0, width=256, height=256, channel=3, model=0,
        lr=0.01, fine_tune_at=0, l2_loss_lambda=0.1, y_class='multi'):
        self.y=y
        self.width=width
        self.height=height
        self.channel=channel
        self.model=model
        self.img_size=[width, height, channel]
        self.lr=lr
        self.fine_tune_at=fine_tune_at
        self.l2_loss_lambda=l2_loss_lambda
        self.y_class=y_class
        
    def model_two(self):
        print("CNN model 2")
        model = Sequential()
        model.add(Conv2D(32, (3, 3), activation='relu', padding='same', name='conv_1',
                         input_shape=self.img_size))
        model.add(MaxPooling2D((2, 2), name='maxpool_1'))
        model.add(Conv2D(64, (3, 3), activation='relu', padding='same', name='conv_2'))
        model.add(MaxPooling2D((2, 2), name='maxpool_2'))
        model.add(Conv2D(128, (3, 3), activation='relu', padding='same', name='conv_3'))
        model.add(MaxPooling2D((2, 2), name='maxpool_3'))
        model.add(Conv2D(128, (3, 3), activation='relu', padding='same', name='conv_4'))
        model.add(MaxPooling2D((2, 2), name='maxpool_4'))
        model.add(Flatten())
        model.add(Dropout(rate=0.25))
        
        model.add(Dense(512, activation='relu', name='dense_1'))
        #model.add(Dropout(rate=0.20))
        model.add(Dense(units=self.y.max()+1, activation='softmax', name='out'))
        
        model.compile(optimizer='adam',
                      loss='categorical_crossentropy',
                      metrics=['accuracy'])
        #optimizer='adam'
        self.model=model
    
    def model_four(self):
        print("CNN model 4")
        model = Sequential()
        model.add(Conv2D(32, (3, 3), activation='relu', padding='same', name='conv_1',
                         input_shape=self.img_size))
        model.add(Conv2D(64, (3, 3), activation='relu', padding='same', name='conv_2'))
        model.add(Conv2D(128, (3, 3), activation='relu', padding='same', name='conv_3'))
        model.add(Conv2D(128, (3, 3), activation='relu', padding='same', name='conv_4'))
        model.add(Dropout(rate=0.20))
        model.add(MaxPooling2D((2, 2), name='maxpool_1'))
        model.add(Flatten())
        
        model.add(Dense(512, activation='relu', name='dense_1'))
        model.add(Dropout(rate=0.20))
        model.add(Dense(512, activation='relu', name='transfer_1'))
        
        if self.y_class == 'multi':
            model.add(Dense(units=self.y.max()+1, activation='softmax', name='out'))
            model.compile(loss='categorical_crossentropy',
                          optimizer=Adam(learning_rate=self.lr),
                          metrics=['accuracy'])
        elif self.y_class == 'binary':
            model.add(Dense(units=1, activation='sigmoid', name='out'))
            model.compile(loss='binary_crossentropy',
                          optimizer=Adam(learning_rate=self.lr),
                          metrics=['accuracy'])
        elif self.y_class == 'regression':
            model.add(Dense(units=1, name='out'))
            model.compile(loss='mean_squared_error',
                          optimizer=Adam(learning_rate=self.lr))
        else:
            print("y_class must be either 'regression' or 'multi' or 'binary'")
        
        self.model=model
    
    def transfer(self):
        print("Transfer learning at layer " + str(self.fine_tune_at))
        regression_model = Sequential()
        for layer in self.model.layers[:-1]: # just exclude last layer from copying
            regression_model.add(layer)
        
        if self.y_class == 'regression':
            regression_model.add(Dense(units=1, name='out'))
        elif self.y_class == 'multi':
            regression_model.add(Dense(self.y.max()+1, activation='softmax'))
        elif self.y_class == 'binary':
            regression_model.add(Dense(units=1, activation='sigmoid', name='out'))
        self.model = regression_model
        
        # Freeze all layers
        self.model.trainable = True
        
        # Make top layers trainable
        for layer in self.model.layers[:self.fine_tune_at]:
          layer.trainable =  False
        
        if self.y_class == 'multi':
            #model.add(Dense(units=self.y.max()+1, activation='softmax', name='out'))
            self.model.compile(loss='categorical_crossentropy',
                          optimizer=Adam(learning_rate=self.lr),
                          metrics=['accuracy'])
        elif self.y_class == 'binary':
            self.model.compile(loss='binary_crossentropy',
                          optimizer=Adam(learning_rate=self.lr),
                          metrics=['accuracy'])
        elif self.y_class == 'regression':
            self.model.compile(loss='mean_squared_error',
                          optimizer=Adam(learning_rate=self.lr))
        else:
            print("y_class must be either 'regression' or 'classification'")


#################
# Visualization #
#################
def plotX(X, outfile):
    plt.imshow(X.reshape(IMG_SIZE, IMG_SIZE,3))
    plt.savefig(outfile)
    plt.close("all")

def plot_loss_accuracy(hist, outfile):
    acc = hist.history['accuracy']
    val_acc = hist.history['val_accuracy']
    loss = hist.history['loss']
    val_loss = hist.history['val_loss']
    
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    ax1.plot(acc, label='Training Accuracy')
    ax1.plot(val_acc, label='Validation Accuracy')
    ax1.legend(loc='lower right')
    ax1.set_ylabel('Accuracy')
    ax1.set_ylim(0, 1.0) #([min(plt.ylim()),1])
    
    ax2.plot(loss, label='Training Loss')
    ax2.plot(val_loss, label='Validation Loss')
    ax2.legend(loc='upper right')
    ax2.set_ylabel('Cross Entropy')
    ax2.set_xlabel('epoch')
    fig.show()
    fig.savefig(outfile)

def plot_confusion_matrix(model, X, y, range_y, outfile):
    y_pred = list(model.predict_classes(X, verbose=0))
    y = list(y)
    plt.figure(figsize=(8, 6))
    
    # Ensure that all the classes are represented
    y.extend(list(range_y))
    y_pred.extend(list(range_y))
    cm=confusion_matrix(y, y_pred)
    np.fill_diagonal(cm, list(cm.diagonal()-1))
    
    # Make the heatmaps!!!
    sns.heatmap(pd.DataFrame(cm), annot=True, fmt='d', cmap='YlGnBu', alpha=0.8, vmin=0)
    plt.savefig(outfile)
    return cm

def plot_class_cm(cm, outfile,metric ='f1'):
    (p, r, f1) = get_F1score(cm)
    if (metric is 'f1'):
        class_frac=f1
    elif (metric is 'p'):
        class_frac = p
    else:
        class_frac = r
    
    class_frac=class_frac.round(2)
    class_id=list(range(0, len(class_frac)))
    class_df = pd.DataFrame({'ID':class_id, 'Frac':class_frac,
                            'Cnt':cm.diagonal(), 'Total':cm.sum(0)})
    
    #b1=plt.bar("ID", "Total", data=class_df,class_id, class_frac)
    b2=plt.bar("ID", "Frac", data=class_df)
    plt.ylim((0,1))
    #plt.rcParams["figure.figsize"] = [6,2]
    plt.xlabel("Cancer_Types")
    plt.ylabel(metric)
    plt.subplots_adjust(bottom=0.2, top=0.8)
    plt.xticks(class_id, rotation=90)
    plt.savefig(outfile)
    
    return class_df

def get_F1score(cm):
    precision = cm.diagonal() / cm.sum(axis=0)
    recall = cm.diagonal() / cm.sum(axis=1)
    f1 = 2 * ((precision * recall)/(precision + recall))
    return (precision, recall, f1)

def increase_brightness(img, value=30):
    hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
    h, s, v = cv2.split(hsv)
    
    lim = 255 - value
    v[v > lim] = 255
    v[v <= lim] += value
    
    final_hsv = cv2.merge((h, s, v))
    img = cv2.cvtColor(final_hsv, cv2.COLOR_HSV2BGR)
    return img

#############
# Load Data #
#############
import imutils
import cv2
import csv
import pandas as pd
import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import pickle
import sys
import seaborn as sns
from itertools import compress

from sklearn.metrics import confusion_matrix, classification_report, r2_score
from sklearn.model_selection import train_test_split

import tensorflow.keras
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import Dense, Dropout, Flatten, Conv2D, MaxPooling2D
from tensorflow.keras.optimizers import Adam

#os.chdir("/cluster/home/quever/git/cnvML")
#from pyimagesearch.gradcam import GradCAM
#from occlusioncnn.occlusion import Occlusion

# VARIABLES
PDIR='/cluster/projects/pughlab/projects/cancer_cell_lines'
SFC='sweep'     # sweep or hilbert
CNTYPE='ASCN'   # TCN or ASCN
model_type='model4'
lr=0.0001

model_type=sys.argv[1]    # 'model4'
lr=sys.argv[2]            # 0.0001
SFC=sys.argv[3]           # sweep or hilbert
CNTYPE=sys.argv[4]        # TCN or ASCN

DATADIR = os.path.join(PDIR, "TCGA", "data", SFC, CNTYPE)
OUTDIR = os.path.join(PDIR, "TCGA", "models", SFC, CNTYPE)
IMG_SIZE=300
CATEGORIES = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD",
              "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC",
              "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",
              "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
              "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM",
              "UCEC", "UCS", "UVM", "Normal"]

# Load data from TCGA Hilberts
pickle_X = open(os.path.join(DATADIR, "X.pickle"), "rb")
pickle_Xids = open(os.path.join(DATADIR, "Xids.pickle"), "rb")
pickle_y = open(os.path.join(DATADIR, "y.pickle"), "rb")
X = pickle.load(pickle_X)
Xids = pickle.load(pickle_Xids)
y = pickle.load(pickle_y)
X,y,Xids = balanceGrp(X, y, Xids, q=0.5)
#plt.imshow(X[2,].reshape(IMG_SIZE, IMG_SIZE,3))
#plt.savefig(os.path.join(OUTDIR, model_type, "class_25.png"))
#plt.close("all")

#X = np.array(X).reshape(-1, IMG_SIZE, IMG_SIZE, 3)

x_train,x_test,y_train,y_test=train_test_split(X,y,test_size=0.2, random_state=1234)

# One-hot encoding of y
y_train_one_hot = tensorflow.keras.utils.to_categorical(y_train,
    num_classes=len(CATEGORIES))
y_test_one_hot = tensorflow.keras.utils.to_categorical(y_test,
    num_classes=len(CATEGORIES))

# Format
x_train = x_train.astype('float32')
x_test = x_test.astype('float32')
x_train = x_train / 255
x_test = x_test / 255

########################
# Build/Train ConvNet #
########################
if not os.path.exists(os.path.join(OUTDIR, model_type, 'my_tcga_model_layer2.h5')):
    M=CNN(y, width=IMG_SIZE, height=IMG_SIZE, channel=3)
    if model_type=='model1':
        M.model_one()
    elif model_type=='model2':
        M.model_two()
    elif model_type=='model3':
        M.model_three()
    elif model_type=='model4':
        M.model_four()
    else:
        print("no model selected")
        
    hist = M.model.fit(x_train, y_train_one_hot, batch_size=32, epochs=10, validation_split=0.2)
    M.model.evaluate(x_test, y_test_one_hot)[1]
    M.model.save(os.path.join(OUTDIR, model_type, 'my_tcga_model_layer2.h5'))
    plot_loss_accuracy(hist, os.path.join(OUTDIR, model_type, "cnn_performance.png"))
    M=M.model
else:
    print("Loading existing model...")
    M = load_model(os.path.join(OUTDIR, model_type, 'my_tcga_model_layer2.h5'))

####################
# Spot-check model #
####################
y_test_class = np.argmax(y_test_one_hot, axis=1)
cm = plot_confusion_matrix(M, x_test, y_test_class, range(0, len(CATEGORIES)),
                           os.path.join(OUTDIR, model_type, "cnn_confusion-matrix.png"))
plt.close("all")

f1 = plot_class_cm(cm, os.path.join(OUTDIR, model_type, "cnn_cm_barplot.png"))
plt.close("all")
np.savetxt(os.path.join(OUTDIR, model_type, 'F1_test.csv'),
    f1, fmt='%5.2f', delimiter=",", header='ID,Frac,Cnt,Total')
print(np.nanmean(f1.Frac)) #model2: 0.749264630372789 ; model3: 0.7793333333333333; model4: 0.7796666666666666
