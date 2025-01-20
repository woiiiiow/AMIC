# -*- coding: utf-8 -*-

# CHEN Junwei, EAP LAP, UC3M, junwei.chen@uc3m.es
# using autoencoder to denoise the PIV fields
# 240510 v4 reduce model to prevent GPU RAM exploding
# 240513 v5 refine
# 240513 v6 transplant to google colab
# 240514 v7 not subtract by mean
# 240514 AEv10
# 240515 change to 2D wing experiment
# import pandas as pd
import tensorflow as tf
from tensorflow import keras
import numpy as np
import time
import random
import h5py
import os
from tensorflow.keras import layers, losses
from tensorflow.keras.models import Model
from scipy.io import loadmat
from scipy.io import savemat

start_time0 = time.time()


# set GPU number
os.environ['CUDA_VISIBLE_DEVICES']='0'
filename, _ = os.path.splitext(os.path.basename(__file__))
SourceFolder = '../../Experiment_wing3/OUT_TRPIV/'
SaveFolder = './'
'''
# google colab utilities
from google.colab import drive
drive.mount('/content/gdrive')
prefix = 'gdrive/MyDrive/COLAB/gauss_noise/'
filename = 'AutoEncoder'
'''

# parameters
AFrame_train = range(1001, 4601)
AFrame = range(11,531)
ROI = [41, 110, 41, 150]
Mm_per_px_ratio = 0.1197
Sample_rate = 30
Vector_Spacing = 10

MyOptimizer   = keras.optimizers.Adam(learning_rate=1e-4)
MyEpochs = 1200
MyBatchSize = 256

print('loading data...')
matdata = loadmat(SourceFolder+'Grid_Wing.mat')
X_etr = matdata['X']
Y_etr = matdata['Y']
X = X_etr[ROI[0]:ROI[1]+1, ROI[2]:ROI[3]+1]
Y = Y_etr[ROI[0]:ROI[1]+1, ROI[2]:ROI[3]+1]
XT = X.transpose()

U_train = np.zeros((len(AFrame_train), np.prod(X.shape)))
V_train = np.zeros((len(AFrame_train), np.prod(X.shape)))
iCount = 0
for iFrame in AFrame_train:
    source = SourceFolder+'Wing_'+str(iFrame).zfill(6)+'.mat'
    matdata = loadmat(source)
    U = matdata['V']
    V = matdata['U']
    U_train[iCount, :] = U[ROI[0]:ROI[1]+1, ROI[2]:ROI[3]+1].transpose().flatten()
    V_train[iCount, :] = V[ROI[0]:ROI[1]+1, ROI[2]:ROI[3]+1].transpose().flatten()
    iCount += 1
U_array = np.zeros((len(AFrame), np.prod(X.shape)))
V_array = np.zeros((len(AFrame), np.prod(X.shape)))
iCount = 0
for iFrame in AFrame:
    source = SourceFolder+'Wing_'+str(iFrame).zfill(6)+'.mat'
    matdata = loadmat(source)
    U = matdata['V']
    V = matdata['U']
    U_array[iCount, :] = U[ROI[0]:ROI[1]+1, ROI[2]:ROI[3]+1].transpose().flatten()
    V_array[iCount, :] = V[ROI[0]:ROI[1]+1, ROI[2]:ROI[3]+1].transpose().flatten()
    iCount += 1

del U, V
total_sam = U_train.shape[0]*2

'''
def AEDeNoise():
    input_layer = keras.Input(shape = XT.shape+(1,))
    x0 = input_layer
    x1 = keras.layers.Conv2D(16, 3, activation='relu', padding='same')(x0)
    x2 = keras.layers.MaxPool2D(2, padding='same')(x1)
    x2 = keras.layers.Conv2D(32, 3, activation='relu', padding='same')(x2)
    x3 = keras.layers.MaxPool2D(2, padding='same')(x2)
    x3 = keras.layers.Conv2D(64, 3, activation='relu', padding='same')(x3)
    z2 = keras.layers.Conv2DTranspose( 32, 3, strides=2, activation='relu', padding='same')(x3)
    z2 = keras.layers.Conv2D(32, 3, activation='relu', padding='same')(z2)
    z1 = keras.layers.Conv2DTranspose( 16, 3, strides=2, activation='relu', padding='same')(z2)
    z1 = keras.layers.Conv2D(16, 3, activation='relu', padding='same')(z1)
    z0 = keras.layers.Conv2D( 1, 3, activation='relu', padding='same')(z1)
    z0 = z0[:,0:x0.shape[1],0:x0.shape[2],:]
    model = tf.keras.Model(input_layer,z0)
    return model
'''

def AEDeNoise():
    input_layer = keras.Input(shape = XT.shape+(1,))
    x0 = input_layer
    x1 = keras.layers.Conv2D(16, 3, activation='LeakyReLU', padding='same')(x0)
    x2 = keras.layers.MaxPool2D(2, padding='same')(x1)
    x2 = keras.layers.Conv2D(32, 3, activation='LeakyReLU', padding='same')(x2)
    z1 = keras.layers.Conv2DTranspose( 16, 3, strides=2, activation='LeakyReLU', padding='same')(x2)
    z1 = keras.layers.Conv2D(16, 3, activation='LeakyReLU', padding='same')(z1)
    z0 = keras.layers.Conv2D( 1, 3, activation='LeakyReLU', padding='same')(z1)
    z0 = z0[:,0:x0.shape[1],0:x0.shape[2],:]
    model = tf.keras.Model(input_layer,z0)
    return model
    
ae_denoise = AEDeNoise()

# ae_denoise.summary()

@tf.function
def train_stage(x):
    with tf.GradientTape() as tape:
        field    = tf.reshape(x, (MyBatchSize,)+XT.shape+(1,))
        field_ae = ae_denoise(field, training=True)
        loss = tf.reduce_sum(tf.square(field_ae-field))
    gradients = tape.gradient(loss, ae_denoise.trainable_variables)
    MyOptimizer.apply_gradients(zip(gradients,ae_denoise.trainable_variables))
    return loss/MyBatchSize

time_elapsed = time.time() - start_time0
print('preparation time of '+filename+': %.1fs' %time_elapsed)

# Commented out IPython magic to ensure Python compatibility.
# training
for epoch in range(MyEpochs):
    start_time = time.time()
    # sample_list = np.array(random.sample(list(range(total_sam)),total_sam))
    sample_list = range(total_sam)
    x0 = np.concatenate((V_train, U_train), axis=0)
    x0 = tf.cast(x0[sample_list,:], tf.float32)
    for batch in range(int(total_sam/MyBatchSize)):
        flag = batch*MyBatchSize
        x    = x0[flag:flag+MyBatchSize,:]
        loss = train_stage(x)
    print('epoch: %d/%d - time: %.1fms/epoch - loss: %f'
              %(epoch+1,MyEpochs,(time.time()-start_time)*1000, loss))

# save result
u_filt = ae_denoise(tf.reshape(U_array, (len(AFrame),)+XT.shape+(1,)))
u_out = np.reshape(u_filt, (len(AFrame),np.size(U_array, axis=1)))*Sample_rate*Mm_per_px_ratio*1e-3
v_filt = ae_denoise(tf.reshape(V_array, (len(AFrame),)+XT.shape+(1,)))
v_out = np.reshape(v_filt, (len(AFrame),np.size(U_array, axis=1)))*Sample_rate*Mm_per_px_ratio*1e-3
data_dict = {
    'AFrame': AFrame,
    'U': np.transpose(u_out),
    'V': np.transpose(v_out),
    }
savemat(SaveFolder+'V_AE.mat', data_dict)

ae_denoise.save_weights(SaveFolder+'AE_weights.h5')

time_elapsed = time.time() - start_time0
print('total time of '+filename+': %.1fs'%time_elapsed)
