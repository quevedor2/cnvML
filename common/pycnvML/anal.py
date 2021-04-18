import numpy as np
import seaborn as sns
import pandas as pd

from sklearn.metrics import confusion_matrix, classification_report, r2_score
from sklearn.model_selection import train_test_split

from pycnvML import viz

def get_F1score(cm):
    precision = cm.diagonal() / cm.sum(axis=0)
    recall = cm.diagonal() / cm.sum(axis=1)
    f1 = 2 * ((precision * recall)/(precision + recall))
    return (precision, recall, f1)

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

def spotcheckModel(M, x_test, y_test_one_hot, CATEGORIES, outpath):
    y_test_class = np.argmax(y_test_one_hot, axis=1)
    
    cm = viz.plot_confusion_matrix(M, x_test, y_test_class, range(0, len(CATEGORIES)),
                                  os.path.join(outpath, "cnn_confusion-matrix.png"))
    plt.close("all")
    
    f1 = viz.plot_class_cm(cm, os.path.join(outpath, "cnn_cm_barplot.png"))
    plt.close("all")
    
    np.savetxt(os.path.join(outpath, 'F1_test.csv'),
        f1, fmt='%5.2f', delimiter=",", header='ID,Frac,Cnt,Total')
    
    m_perf = np.nanmean(f1.Frac)
    print(m_perf) #model2: 0.749264630372789 ; model3: 0.7793333333333333; model4: 0.7796666666666666
    return m_perf
