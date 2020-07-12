# import the necessary packages
from tensorflow.keras.models import Model
import tensorflow as tf
import numpy as np
import cv2

class Occlusion:
    def __init__(self, model, img, IMG_SIZE):
        self.model = model
        self.img = img
        self.IMG_SIZE = IMG_SIZE
        self.occ_array = None
        self.occ_prob = None
        self.occ_idx = None
        self.occ_viz = None
        self.step=None
    
    def create_occlusion_dataset(self, step, mem='high'):
        occ_idx=list(range(0, self.IMG_SIZE, step))
        if(mem=='low'):
            # Does predictiosn in real time rather than as an array at the end
            occ_ret = np.zeros(((len(occ_idx)*len(occ_idx)),
                int(self.model.output_shape[1])))
        else:
            #initializes occlusion matrix to predict on at the end
            occ_ret=np.zeros((len(occ_idx)*len(occ_idx),
                self.IMG_SIZE, self.IMG_SIZE, 3))
        
        arr_ind = 0
        for ypos in occ_idx:
            for xpos in occ_idx:
                occ_tmp=self.img.copy()
                # Blank out box area
                occ_tmp[xpos:(xpos+step-1), ypos:(ypos+step-1),] = 1
                if(mem=='low'):
                    occ_ret[arr_ind,:] = self.model.predict_proba(
                        np.expand_dims(occ_tmp, axis=0), verbose=False)
                else:
                    occ_ret[arr_ind,:,:,:] = occ_tmp
                    arr_ind=arr_ind+1
        
        self.occ_idx = occ_idx
        self.step = step
        self.occ_array = occ_ret
    
    def overlay_heatmap(self, heatmap, image, alpha=0.5, colormap=cv2.COLORMAP_VIRIDIS):
        heatmap = cv2.applyColorMap(heatmap, colormap)
        output = cv2.addWeighted(image, alpha, heatmap, 1 - alpha, 0)
        return (heatmap, output)
    
    def predict_occlusion(self):
        occ_prob=self.model.predict_proba(self.occ_array, verbose=False)
        self.occ_prob = occ_prob
    
    def fill_occ_array(self, ctype):
        # Initialize zeroes np.array
        occ_viz=np.zeros((self.IMG_SIZE, self.IMG_SIZE,3))
        arr_ind=0
        for ypos in self.occ_idx:
            for xpos in self.occ_idx:
                yidx=ypos if ypos == 0 else ypos - 1
                xidx=xpos if xpos == 0 else xpos - 1
                occ_viz[xidx:(xpos+self.step-1),
                    yidx:(ypos+self.step-1),] = self.occ_prob[arr_ind, ctype]
                arr_ind=arr_ind+1
        self.occ_viz = occ_viz
    
    def winsorize_heatmap(self, winsor_lo=0.05, winsor_hi=0.95):
        numer = self.occ_viz - np.quantile(self.occ_viz, winsor_lo)
        numer = np.where(numer < 0, 0, numer)
        numer = np.abs(1-numer)
        
        denom = (np.quantile(self.occ_viz, winsor_hi) - np.quantile(self.occ_viz, winsor_lo)) + 1e-8
        denom = np.where(denom > 1, 1, denom)
        
        heatmap = cv2.resize(((numer / denom)*255).astype('uint8'), (self.IMG_SIZE, self.IMG_SIZE))
        return(heatmap[:,:,0])
    
    def quantile_match(self, hmap1, hmap2):
        percentiles = np.argsort(np.argsort(c1)) * 100. / (len(c1) - 1)
        c0 = np.quantile(hmap2, np.round(percentiles/100, 2))
        hmap0 = c0.reshape((IMG_SIZE, IMG_SIZE))
        hmap0 = hmap0.astype("uint8")
        return hmap0
