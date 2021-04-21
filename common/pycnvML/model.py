import os
import tensorflow.keras
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import Dense, Dropout, Flatten, Conv2D, MaxPooling2D, BatchNormalization
from tensorflow.keras.optimizers import Adam
from pycnvML import viz

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
    
    def model_alexnet(self):
        print("AlexNet-like model")
        model = Sequential()
        model.add(Conv2D(filters=96, kernel_size=(11, 11), strides=(4,4), activation='relu',
                         padding='same', name='conv_1', input_shape=self.img_size))
        model.add(BatchNormalization())
        model.add(MaxPooling2D(pool_size=(3,3), strides=(2,2)))
        model.add(Conv2D(filters=256, kernel_size=(5,5), strides=(1,1), activation='relu',
                         padding='same', name='conv_2'))
        model.add(BatchNormalization())
        model.add(MaxPooling2D(pool_size=(3,3), strides=(2,2)))
        model.add(Conv2D(filters=384, kernel_size=(3,3), strides=(1,1), activation='relu',
                         padding='same', name='conv_3'))
        model.add(BatchNormalization())
        model.add(Conv2D(filters=384, kernel_size=(3,3), strides=(1,1), activation='relu',
                         padding='same', name='conv_4'))
        model.add(BatchNormalization())
        model.add(Conv2D(filters=256, kernel_size=(3,3), strides=(1,1), activation='relu',
                         padding='same', name='conv_5'))
        model.add(MaxPooling2D((3, 3), strides=(2,2), name='maxpool_1'))
        model.add(Flatten())
        model.add(Dense(4096, activation='relu', name='dense_1'))
        model.add(Dropout(rate=0.5))
        model.add(Dense(512, activation='relu', name='transfer_1'))
        model.add(Dropout(rate=0.5))
        
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


def buildModel(y, IMG_SIZE, lr, model_type, x_train, y_train_one_hot, epochs,
               x_test, y_test_one_hot, outpath):
    ########################
    # Build/Train ConvNet #
    ########################
    if not os.path.exists(os.path.join(outpath, 'model_' + str(lr) + '.h5')):
        M=CNN(y, width=IMG_SIZE, height=IMG_SIZE, channel=3, lr=float(lr))
        if model_type=='model1':
            M.model_one()
        elif model_type=='model2':
            M.model_two()
        elif model_type=='alexnet':
            M.model_alexnet()
        elif model_type=='model4':
            M.model_four()
        else:
            print("no model selected")
        
        hist = M.model.fit(x_train, y_train_one_hot, batch_size=32, epochs=epochs, validation_split=0.2)
        M.model.evaluate(x_test, y_test_one_hot)[1]
        M.model.save(os.path.join(outpath, 'model_' + str(lr) + '.h5'))
        print("Model saved: " + os.path.join(outpath, 'model_' + str(lr) + '.h5'))
        viz.plot_loss_accuracy(hist, os.path.join(outpath, "cnn_performance.png"))
        M=M.model
    else:
        print("Loading existing model...")
        M = load_model(os.path.join(outpath, 'model_' + str(lr) + '.h5'))
    return M
