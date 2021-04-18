import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import confusion_matrix
import seaborn as sns
import pandas as pd
import cv2
from pycnvML import anal


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
    plt.close("all")
    return cm

def plot_class_cm(cm, outfile,metric ='f1'):
    (p, r, f1) = anal.get_F1score(cm)
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
    plt.close("all")
    
    return class_df

def increase_brightness(img, value=30):
    hsv = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
    h, s, v = cv2.split(hsv)
    
    lim = 255 - value
    v[v > lim] = 255
    v[v <= lim] += value
    
    final_hsv = cv2.merge((h, s, v))
    img = cv2.cvtColor(final_hsv, cv2.COLOR_HSV2BGR)
    return img
