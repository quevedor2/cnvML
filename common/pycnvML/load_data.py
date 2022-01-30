import os
import pickle
from itertools import compress
from sklearn.model_selection import train_test_split
import tensorflow.keras

from pycnvML.anal import *

def readPickle(PDIR, DATASET, SFC, CNTYPE, CATEGORIES, IMG_SIZE=300, CCL_DATASET='', CCL_TRAIN=''):
    DATADIR = os.path.join(PDIR, DATASET, "data", SFC, CNTYPE, CCL_DATASET)
    OUTDIR = os.path.join(PDIR, DATASET, "models", SFC, CNTYPE, CCL_DATASET)
    TCGADIR = os.path.join(PDIR, 'TCGA', "models", SFC, CNTYPE, '')
    CCLDIR = os.path.join(PDIR, 'CCL', "models", SFC, CNTYPE, CCL_DATASET)
    CCLTRAIN = os.path.join(PDIR, 'CCL', "models", SFC, CNTYPE, CCL_TRAIN)
    
    try:
        os.makedirs(TCGADIR)
        os.makedirs(CCLDIR)
        os.makedirs(CCLTRAIN)
    except OSError as error:
        print(error)
    IMG_SIZE=IMG_SIZE
    
    # Load data from TCGA Hilberts
    pickle_X = open(os.path.join(DATADIR, "X.pickle"), "rb")
    pickle_Xids = open(os.path.join(DATADIR, "Xids.pickle"), "rb")
    pickle_y = open(os.path.join(DATADIR, "y.pickle"), "rb")
    X = pickle.load(pickle_X)
    Xids = pickle.load(pickle_Xids)
    y = pickle.load(pickle_y)
    
    return(X, Xids, y, DATADIR, TCGADIR, CCLDIR, CCLTRAIN, OUTDIR)

def balanceAndFormatData(X, y, Xids, CATEGORIES):
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
    
    return(x_train, x_test, y_train_one_hot, y_test_one_hot)

def matchCategories(ccl_dir, X, Xids, y, CATEGORIES):
    ccl_CAT=[]
    for o in os.listdir(ccl_dir):
        if os.path.isdir(os.path.join(ccl_dir, o)):
            ccl_CAT.append(o)
    
    # Find overlaps between CATEGORIES and the listed categories in ccl_dir
    ccl_y_verbose = [ccl_CAT[int(i)] for i in y]
    ccl_tcga_ov = [elem in CATEGORIES for elem in ccl_y_verbose]
    
    # Reduce to overlapping categories
    X = X[ccl_tcga_ov,:]
    Xids = Xids.loc[ccl_tcga_ov,:]
    yT = list(compress(ccl_y_verbose, ccl_tcga_ov))
    y = [ CATEGORIES.index(x) for x in yT ]
    y = np.c_[y]
    return(X, Xids, y, ccl_tcga_ov)
