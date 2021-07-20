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


import os, sys, getopt, argparse
from itertools import compress
from pathlib import Path
from pycnvML import load_data, model, anal
#from pyimagesearch.gradcam import GradCAM
#from occlusioncnn.occlusion import Occlusion

def main(argv):
    IMG_SIZE=300
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--model', dest="model_type",
                    required=False, type=str,
                    metavar="<model>", help="name of the model architecture to use: alexnet or model4")
    parser.add_argument('-l', '--lr', dest="lr",
                    required=False, type=float, default=0.0001,
                    metavar="<lr>", help="Learning rate")
    parser.add_argument('-s', '--sfc', dest="SFC",
                    required=False, type=str, default='hilbert',
                    metavar="<sfc>", help="Space-filling curve (hilbert, morton, random, sweep, scan, diagonal)")
    parser.add_argument('-c', '--cntype', dest="CNTYPE",
                    required=False, type=str, default='ASCN',
                    metavar="<cntype>", help="ASCN or TCN")
    parser.add_argument('-d', '--dataset', dest="DATASET",
                    required=False, type=str, default='TCGA',
                    metavar="<dataset>", help="Either TCGA or CCL")
    parser.add_argument('-e', '--epochs', dest="EPOCHS",
                    required=False, type=int, default=70,
                    metavar="<epochs>", help="Number of epochs to run for training")
    parser.add_argument('-w', '--pdir', dest="PDIR",
                    required=False, type=str, default='/cluster/projects/pughlab/projects/cancer_cell_lines',
                    metavar="<pdir>", help="Parent directory containing the CCL and TCGA directories")
    parser.add_argument('-a', '--analysis', dest="ANALYSIS",
                    required=False, type=str, default='naive',
                    metavar="<analysis>", help="Whether to run a 'naive' training or 'transfer' from an existing model")
    parser.add_argument('-p', '--cclds', dest="CCL_DATASET",
                    required=False, type=str, default='',
                    metavar="<ccl_dataset>", help="What CCL dataset to use: GDCS, CCLE or GNE")
    args = parser.parse_args()
    
    #args.CCL_DATASET='GDSC'
    #args.ANALYSIS='transfer'
    #args.DATASET='CCL'
    
    #args.CCL_DATASET=''
    #args.ANALYSIS='naive'
    #args.DATASET='TCGA'
    
    CATEGORIES = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD",
        "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC",
        "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",
        "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
        "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM",
        "UCEC", "UCS", "UVM", "Normal"]


    (X, Xids, y, DATADIR, TCGADIR, CCLDIR, OUTDIR) = load_data.readPickle(args.PDIR, args.DATASET, args.SFC,
        args.CNTYPE, CATEGORIES, IMG_SIZE=IMG_SIZE, CCL_DATASET=args.CCL_DATASET)
    if args.DATASET == 'CCL':
        ccl_dir=os.path.join(args.PDIR, args.DATASET, 'data', args.SFC, args.CNTYPE, args.CCL_DATASET)
        (X, Xids, y, ov_idx) = load_data.matchCategories(ccl_dir, X, Xids, y, CATEGORIES)
    
    (x_train, x_test, y_train_one_hot, y_test_one_hot) = load_data.balanceAndFormatData(X, y, Xids, CATEGORIES)
    
    
    # OUTDIR = os.path.join(PDIR, DATASET, "models", SFC, CNTYPE, CCL_DATASET)
    # OUTDIR = os.path.join("/cluster/projects/pughlab/projects/cancer_cell_lines",
    #    "CCL", "models", "hilbert", "ASCN", "GDSC")
    model_path = os.path.join(OUTDIR, args.model_type)
    tcga_model_path = os.path.join(TCGADIR, args.model_type)
    ccl_model_path = os.path.join(CCLDIR, args.model_type)
    if args.ANALYSIS=='naive':
        M = model.buildModel(y, IMG_SIZE, args.lr, args.model_type, x_train, y_train_one_hot,
            args.EPOCHS, x_test, y_test_one_hot, model_path)
    elif args.ANALYSIS=='transfer':
        M = model.transferModel(y, IMG_SIZE, args.lr, x_train, y_train_one_hot, args.EPOCHS,
            x_test, y_test_one_hot, tcga_model_path, ccl_model_path, CATEGORIES, layer=9)
    
    m_perf = anal.spotcheckModel(M, x_test, y_test_one_hot, CATEGORIES, model_path)


if __name__ == "__main__":
    main(sys.argv[1:])
