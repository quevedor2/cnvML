#!/usr/bin/python

'''
DESCRIPTION:
	Takes the outputted directory structure from assignOncocode.R, which
	uses the oncocode metadata for each file to assign each SFC png to a
	separate directory based on its oncocode. This code will then create a
	pickle file to be passed into downstream deep learning applications.

USAGE:
	pickle_hilbert.py -a <analysis> -s <sfc> -c <cntype> -d <dir>

EXAMPLES:
	pickle_hilbert.py -a TCGA -s sweep -c ASCN -d /cluster/projects/pughlab/projects/cancer_cell_lines/
	pickle_hilbert.py -a TCGA -s hilbert -c ASCN -d /cluster/projects/pughlab/projects/cancer_cell_lines/
	pickle_hilbert.py -a TCGA -s hilbert -c TCN -d /cluster/projects/pughlab/projects/cancer_cell_lines/

	pickle_hilbert.py -a CCL -s sweep -c ASCN -d /cluster/projects/pughlab/projects/cancer_cell_lines/
	pickle_hilbert.py -a CCL -s hilbert -c ASCN -d /cluster/projects/pughlab/projects/cancer_cell_lines/
	pickle_hilbert.py -a CCL -s hilbert -c TCN -d /cluster/projects/pughlab/projects/cancer_cell_lines/
'''

import os, sys, getopt
import cv2, random, pickle, re
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

def pickleit(DATADIRS, CATEGORIES, IMG_SIZE, seed, analysis):
	for DATADIR in DATADIRS:
		print(DATADIR)
		training_data = []
		
		def create_training_data(CANCER_TYPES):
			for category in CANCER_TYPES :
				path = os.path.join(DATADIR, category)
				class_num = CANCER_TYPES.index(category)
				for img in os.listdir(path):
					try :
						img_array = cv2.imread(os.path.join(path, img), cv2.IMREAD_COLOR) # (300,300,3)
						new_array = cv2.resize(img_array, (IMG_SIZE, IMG_SIZE), interpolation = cv2.INTER_AREA)
						id = re.sub(".png$", "", img)
						training_data.append([new_array, class_num, id, category])
					except Exception as e:
						pass
		
		if(str(analysis) == 'CCL'):
			CATEGORIES=[]
			for o in os.listdir(DATADIR):
				if os.path.isdir(os.path.join(DATADIR, o)):
					CATEGORIES.append(o)
		
		create_training_data(CATEGORIES)
		random.Random(seed).shuffle(training_data,)
		
		X = [] #features
		y = [] #labels
		Xids = [] #Sample IDs
		yids = [] #Category IDs (descriptive labels)
		
		for features, label, sample, category in training_data:
			X.append(features)
			y.append(label)
			Xids.append(sample)
			yids.append(category)
		
		X = np.array(X).reshape(-1, IMG_SIZE, IMG_SIZE, 3)
		y = np.array(y).reshape(-1, 1)
		Xids = np.array(Xids).reshape(-1, 1)
		yids = np.array(yids).reshape(-1, 1)
		
		img = plt.imshow(X[1])
		plt.savefig(os.path.join(DATADIR, "test.png"))
		
		# Creating the files containing all the information about your model
		pickle_out = open(os.path.join(DATADIR, "X.pickle"), "wb")
		pickle.dump(X, pickle_out)
		pickle_out.close()
		
		pickle_out = open(os.path.join(DATADIR, "y.pickle"), "wb")
		pickle.dump(y, pickle_out)
		pickle_out.close()
		
		Xdf = pd.DataFrame(data=Xids, columns=['samples'])
		ydf = pd.DataFrame(data=yids, columns=['cancer_type'])
		Xdf.to_pickle(os.path.join(DATADIR, "Xids.pickle"))
		ydf.to_pickle(os.path.join(DATADIR, "yids.pickle"))

def main(argv):
	analysis = '' 	# 'TCGA' or 'CCL'
	sfc = '' 		# 'sweep' or 'hilbert'
	cntype = ''		# 'TCN' or 'ASCN' (TCN is only for hilbert)
	dir='' 			# '/cluster/projects/pughlab/projects/cancer_cell_lines/'
	
	try:
		opts, args = getopt.getopt(argv,"ha:s:c:d:",["analysis=","sfc=","cntype=","dir="])
	except getopt.GetoptError:
		print('pickle_hilbert.py -a <analysis> -s <sfc> -c <cntype> -d <dir>')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print('pickle_hilbert.py -a <analysis> -s <sfc> -c <cntype> -d <dir>')
			sys.exit()
		elif opt in ("-a", "--analysis"):
			analysis = arg
		elif opt in ("-s", "--sfc"):
			sfc = arg
		elif opt in ("-c", "--cntype"):
			cntype = arg
		elif opt in ("-d", "--dir"):
			dir = arg
	
	PDIR=os.path.join(dir, analysis)
	
	if(str(analysis) == 'CCL'):
		DATADIRS=[os.path.join(PDIR, "data", sfc, cntype, x) for x in ['GDSC', 'CCLE', 'GNE']]
		CATEGORIES=[]
	else:
		DATADIRS=[os.path.join(PDIR, "data", sfc, cntype)]
		# All the categories you want your neural network to detect
		CATEGORIES = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD",
					  "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC",
					  "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",
					  "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
					  "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM",
					  "UCEC", "UCS", "UVM", "Normal"]

	# The size of the images that your neural network will use
	IMG_SIZE = 300
	seed=1234
	pickleit(DATADIRS, CATEGORIES, IMG_SIZE, seed, analysis)

if __name__ == "__main__":
	main(sys.argv[1:])
