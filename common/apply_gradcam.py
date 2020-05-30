# USAGE
# python apply_gradcam.py --image images/space_shuttle.jpg
# python apply_gradcam.py --image images/beagle.jpg
# python apply_gradcam.py --image images/soccer_ball.jpg --model resnet

# import the necessary packages
from tensorflow.keras.applications import ResNet50, VGG16, imagenet_utils
from tensorflow.keras.preprocessing.image import img_to_array, load_img
from tensorflow.keras.models import load_model

from sklearn.model_selection import train_test_split
import numpy as np
import argparse
import imutils
import cv2
import os
import pickle


os.chdir("/cluster/home/quever/git/cnvML")
from pyimagesearch.gradcam import GradCAM

# construct the argument parser and parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--image", required=True,
	help="path to the input image")
ap.add_argument("-m", "--model", type=str, default="vgg",
	choices=("vgg", "resnet"),
	help="model to be used")
ap.add_argument("-d", "--pdir", type=str, default="./",
	help="Parent directory containing the 'models' directory")
ap.add_argument("-j", "--class_j", type=int, default=int(0),
	choices=range(0,32), metavar="[0-32]",
	help="Cancer type class to pick (int)")
args = vars(ap.parse_args())
# LIMPS_p_NCLE_DNA2N_GenomeWideSNP_6_G12_246776.png - MDAMB231_BREAST breast carcinoma
# /cluster/home/quever/xfer/ARLES_p_NCLE_DNAAffy2_S_GenomeWideSNP_6_A10_256020.png - CaOV-3 high grade serous ovarian
args = {'image':'/cluster/home/quever/xfer/LIMPS_p_NCLE_DNA2N_GenomeWideSNP_6_G12_246776.png',
  'pdir':'/cluster/projects/pughlab/projects/cancer_cell_lines/TCGA',
  'model':'model1', 'class_j':19}

## Initialize the model
OUTDIR = os.path.join(args['pdir'], "models")
#Model=load_model(os.path.join(OUTDIR, "binary", args['model'], str(args['class_j']), 'my_tcga_model_layer2.h5'))
Model=load_model(os.path.join(OUTDIR, args['model'], 'my_tcga_model_layer2.h5'))

## Load data from a pickle file of TCGA Hilberts
IMG_SIZE=300
#DATADIR = os.path.join(args['pdir'], "data")
#pickle_X = open(os.path.join(DATADIR, "X.pickle"), "rb")
#pickle_y = open(os.path.join(DATADIR, "y.pickle"), "rb")
#X = pickle.load(pickle_X)
#y = pickle.load(pickle_y)

#tumor_idx=np.argwhere(y==int(args['class_j']))[:,0]
#Xsub=X[tumor_idx,:]
#Xsub = Xsub.astype('float32')
#Xsub = Xsub / 255

orig = cv2.imread(args["image"])
resized = cv2.resize(orig, (IMG_SIZE, IMG_SIZE))

# load the input image from disk (in Keras/TensorFlow format) and
# preprocess it
image = load_img(args["image"], target_size=(IMG_SIZE, IMG_SIZE))
image = img_to_array(image)
image = np.expand_dims(image, axis=0)
image = image.astype('float32')
image = image / 255

preds = Model.predict(image)
i = np.argmax(preds[0])

# initialize our gradient class activation map and build the heatmap
for each_i in list(set([i, int(args['class_j'])])):
	print(str(each_i) + '-' + str(args['class_j']))
	cam = GradCAM(Model, each_i)
	heatmap = cam.compute_heatmap(image)
	
	# resize the resulting heatmap to the original input image dimensions
	# and then overlay heatmap on top of the image
	heatmap = cv2.resize(heatmap, (IMG_SIZE, IMG_SIZE))
	(heatmap, output) = cam.overlay_heatmap(heatmap, orig, alpha=0.5)
	
	# draw the predicted label on the output image
	#cv2.rectangle(output, (0, 0), (340, 40), (0, 0, 0), -1)
	#cv2.putText(output, str(i), (10, 25), cv2.FONT_HERSHEY_SIMPLEX,
	#	0.8, (255, 255, 255), 2)
	
	# display the original image and resulting heatmap and output image
	# to our screen
	output = np.vstack([orig, heatmap, output])
	output = imutils.resize(output, height=700)
	#cv2.imshow("Output", output)
	cv2.imwrite(os.path.join(OUTDIR, args['model'], 'img_' + str(each_i) + '-' + str(args['class_j']) + '.jpg'), output)
