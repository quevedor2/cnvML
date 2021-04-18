import os, sys, getopt
from pycnvML import load_data, model, anal
#from pyimagesearch.gradcam import GradCAM
#from occlusioncnn.occlusion import Occlusion

def main(argv):
    PDIR='/cluster/projects/pughlab/projects/cancer_cell_lines'
    SFC='sweep'     # sweep or hilbert
    CNTYPE='ASCN'   # TCN or ASCN
    model_type='model4'
    lr=0.0001
    
    model_type=sys.argv[1]    # 'model4'
    lr=sys.argv[2]            # 0.0001
    SFC=sys.argv[3]           # sweep or hilbert
    CNTYPE=sys.argv[4]        # TCN or ASCN
    DATASET=sys.argv[5]       # TCGA or CCL
    
	try:
		opts, args = getopt.getopt(argv,"hm:l:s:c:d:",["model=","lr=","sfc=","cntype=","dataset="])
	except getopt.GetoptError:
		print('build_cnn_model.py -m <model> -l <lr> -s <sfc> -c <cntype> -d <dir>')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print('build_cnn_model.py -m <model> -l <lr> -s <sfc> -c <cntype> -d <dir>')
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
    
    CATEGORIES = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD",
                  "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC",
                  "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",
                  "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
                  "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM",
                  "UCEC", "UCS", "UVM", "Normal"]
    
    (X, Xids, y) = load_data.readPickle(PDIR, DATASET, SFC, CNTYPE, CATEGORIES)
    
    (x_train, x_test, y_train_one_hot, y_test_one_hot) = load_data.balanceAndFormatData(X, y, Xids, CATEGORIES)
    
    model_path = os.path.join(OUTDIR, model_type)
    M = model.buildModel(y, IMG_SIZE, lr, model_type, x_train, y_train_one_hot,
                         x_test, y_test_one_hot, model_path)
    m_perf = anal.spotcheckModel(M, x_test, y_test_one_hot, CATEGORIES, model_path)


if __name__ == "__main__":
	main(sys.argv[1:])
