#!/usr/bin/python

'''
DESCRIPTION:
	Takes the pickle files outputted from pickle_hilbert.py. It will use Keras
    tensorflow to build a CNN using these pickle files. The helper functions are
    found in the pycnvML modules, while the CNN models can be found in the
    pycnvML/model.py module. This script will run the pycnvML package to
    both build the model and do a basic performance assessment of Oncocode
    predictions.

USAGE:
	build_cnn_model.py -m <model> -l <lr> -s <sfc> -c <cntype> -d <dir>

EXAMPLES:
    # runBuildCnnModel.sh:
    python build_cnn_model.py -m <model_type> -l <lr> -s <SFC> -c <CNTYPE> -d <DATASET>
    python ~/git/cnvML/common/build_cnn_model.py -m $1 -l $2 -s $3 -c $4 -d $5
    
    # queueBuildCnnModel.sh
	sbatch -J tcga_sweep_ascn_model4_0.0005 runBuildCnnModel.sh model4 0.0005 sweep ASCN TCGA
    sbatch -J tcga_hilbert_tcn_model4_0.0005 runBuildCnnModel.sh model4 0.0005 hilbert TCN TCGA
    sbatch -J tcga_hilbert_ascn_model4_0.0005 runBuildCnnModel.sh model4 0.0005 hilbert ASCN TCGA
'''


import os, sys, getopt
from pycnvML import load_data, model, anal
#from pyimagesearch.gradcam import GradCAM
#from occlusioncnn.occlusion import Occlusion

def main(argv):
    PDIR='/cluster/projects/pughlab/projects/cancer_cell_lines'
    IMG_SIZE=300
    
    SFC='sweep'     # sweep or hilbert
    CNTYPE='ASCN'   # TCN or ASCN
    model_type='model4'
    lr=0.0001
	EPOCHS=10
    
    model_type=sys.argv[1]    # 'model4'
    lr=sys.argv[2]            # 0.0001
    SFC=sys.argv[3]           # sweep or hilbert
    CNTYPE=sys.argv[4]        # TCN or ASCN
    DATASET=sys.argv[5]       # TCGA or CCL
	EPOCHS=sys.argv[6]       # TCGA or CCL
    
    try:
        opts, args = getopt.getopt(argv,"hm:l:s:c:d:e:",["model=","lr=","sfc=","cntype=","dataset=","epochs="])
    except getopt.GetoptError:
        print('build_cnn_model.py -m <model> -l <lr> -s <sfc> -c <cntype> -d <dir> -e <epochs>')
        sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-h':
            print('build_cnn_model.py -m <model> -l <lr> -s <sfc> -c <cntype> -d <dir> -e <epochs>')
            sys.exit()
        elif opt in ("-m", "--model"):
            model_type = arg
        elif opt in ("-l", "--lr"):
            lr = arg
        elif opt in ("-s", "--sfc"):
            SFC = arg
        elif opt in ("-c", "--cntype"):
            CNTYPE = arg
        elif opt in ("-d", "--dataset"):
            DATASET = arg
		elif opt in ("-e", "--epochs"):
            EPOCHS = arg
    
    CATEGORIES = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD",
        "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC",
        "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",
        "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
        "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM",
        "UCEC", "UCS", "UVM", "Normal"]
    
    (X, Xids, y, DATADIR, OUTDIR) = load_data.readPickle(PDIR, DATASET, SFC, CNTYPE, CATEGORIES, IMG_SIZE=IMG_SIZE)
    
    (x_train, x_test, y_train_one_hot, y_test_one_hot) = load_data.balanceAndFormatData(X, y, Xids, CATEGORIES)
    
    model_path = os.path.join(OUTDIR, model_type)
    M = model.buildModel(y, IMG_SIZE, lr, model_type, x_train, y_train_one_hot, EPOCHS,
    	x_test, y_test_one_hot, model_path)
    m_perf = anal.spotcheckModel(M, x_test, y_test_one_hot, CATEGORIES, model_path)


if __name__ == "__main__":
    main(sys.argv[1:])
