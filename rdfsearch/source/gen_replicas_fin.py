# Written by In-Ho Lee, KRISS, November 19, 2019.
# nohup nice python gen_replicas.py &> replicas.out &
#
from os import listdir, getcwd
from os.path import isfile, join
import pathlib
from pathlib import Path
from datetime import datetime
import time

"""
Title: Variational AutoEncoder
Author: [fchollet](https://twitter.com/fchollet)
Date created: 2020/05/03
Last modified: 2020/05/03
Description: Convolutional Variational AutoEncoder (VAE) trained on MNIST digits.
"""

"""
## Setup
"""

import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
#from keras import optimizers
from scipy.stats import norm


"""
## Create a sampling layer
"""
class Sampling(layers.Layer):
    """Uses (z_mean, z_log_var) to sample z, the vector encoding a digit."""

    def call(self, inputs):
        z_mean, z_log_var = inputs
        batch = tf.shape(z_mean)[0]
        dim = tf.shape(z_mean)[1]
        epsilon = tf.keras.backend.random_normal(shape=(batch, dim))
        return z_mean + tf.exp(0.5 * z_log_var) * epsilon

"""
## Build the encoder
"""

original_dim =200*9
ngrid =200*9
latent_dim = 2
icase=4

if icase == 1:
    encoder_inputs = keras.Input(shape=(original_dim,))
    x = layers.Dense(100, activation="relu")(encoder_inputs)
    x = layers.Dense(50, activation="relu")(x)
    x = layers.Dense(25, activation="relu")(x)
    x = layers.Dense(latent_dim)(x)
if icase == 2:
    encoder_inputs = keras.Input(shape=(original_dim,))
    x = layers.Dense(100, activation="relu")(encoder_inputs)
    x = layers.Dense(50, activation="relu")(x)
    x = layers.Dense( 5, activation="relu")(x)
    x = layers.Dense(latent_dim)(x)
if icase == 3:
    encoder_inputs = keras.Input(shape=(original_dim,))
    x = layers.Dense(50, activation="relu")(encoder_inputs)
    x = layers.Dense(25, activation="relu")(x)
    x = layers.Dense(12, activation="relu")(x)
    x = layers.Dense(latent_dim)(x)
if icase == 4:
    encoder_inputs = keras.Input(shape=(original_dim,))
    x = layers.Dense(20, activation="relu")(encoder_inputs)
#   x = layers.Dense(20, activation="relu")(x)
    x = layers.Dense(10, activation="relu")(x)
    x = layers.Dense(latent_dim)(x)

z_mean = layers.Dense(latent_dim, name="z_mean")(x)
z_log_var = layers.Dense(latent_dim, name="z_log_var")(x)
z = Sampling()([z_mean, z_log_var])
encoder = keras.Model(encoder_inputs, [z_mean, z_log_var, z], name="encoder")
encoder.summary()

"""
## Build the decoder
"""

if icase == 1:
    latent_inputs = keras.Input(shape=(latent_dim,))
    x = layers.Dense(25, activation="relu")(latent_inputs)
    x = layers.Dense(50, activation="relu")(x)
    x = layers.Dense(100, activation="relu")(x)
if icase == 2:
    latent_inputs = keras.Input(shape=(latent_dim,))
    x = layers.Dense( 5, activation="relu")(latent_inputs)
    x = layers.Dense(50, activation="relu")(x)
    x = layers.Dense(100, activation="relu")(x)
if icase == 3:
    latent_inputs = keras.Input(shape=(latent_dim,))
    x = layers.Dense(12, activation="relu")(latent_inputs)
    x = layers.Dense(25, activation="relu")(x)
    x = layers.Dense(50, activation="relu")(x)
if icase == 4:
    latent_inputs = keras.Input(shape=(latent_dim,))
    x = layers.Dense(10, activation="relu")(latent_inputs)
#   x = layers.Dense(20, activation="relu")(x)
    x = layers.Dense(20, activation="relu")(x)
decoder_outputs = layers.Dense(original_dim, activation="sigmoid" )(x)
decoder = keras.Model(latent_inputs, decoder_outputs, name="decoder")
decoder.summary()

"""
## Define the VAE as a `Model` with a custom `train_step`
"""

class VAE(keras.Model):
    def __init__(self, encoder, decoder, **kwargs):
        super(VAE, self).__init__(**kwargs)
        self.encoder = encoder
        self.decoder = decoder
        self.total_loss_tracker = keras.metrics.Mean(name="total_loss")
        self.reconstruction_loss_tracker = keras.metrics.Mean(
            name="reconstruction_loss"
        )
        self.kl_loss_tracker = keras.metrics.Mean(name="kl_loss")

    @property
    def metrics(self):
        return [
            self.total_loss_tracker,
            self.reconstruction_loss_tracker,
            self.kl_loss_tracker,
        ]

    def train_step(self, data):
        with tf.GradientTape() as tape:
            z_mean, z_log_var, z = self.encoder(data)
            reconstruction = self.decoder(z)
#           print(reconstruction.shape,data.shape)
            reconstruction_loss = tf.reduce_mean(
                tf.reduce_sum(
#                   keras.losses.binary_crossentropy(data, reconstruction), axis=(1, 2)
                    keras.losses.binary_crossentropy(data, reconstruction), axis=-1
                )
            )
            kl_loss = -0.5 * (1 + z_log_var - tf.square(z_mean) - tf.exp(z_log_var))
#           print(kl_loss.shape)
            kl_loss = tf.reduce_mean(tf.reduce_sum(kl_loss, axis=1))
            total_loss = reconstruction_loss + kl_loss
        grads = tape.gradient(total_loss, self.trainable_weights)
        self.optimizer.apply_gradients(zip(grads, self.trainable_weights))
        self.total_loss_tracker.update_state(total_loss)
        self.reconstruction_loss_tracker.update_state(reconstruction_loss)
        self.kl_loss_tracker.update_state(kl_loss)
        return {
            "loss": self.total_loss_tracker.result(),
            "reconstruction_loss": self.reconstruction_loss_tracker.result(),
            "kl_loss": self.kl_loss_tracker.result(),
        }

"""
## Train the VAE
"""

if False:
    (x_train, _), (x_test, _) = keras.datasets.mnist.load_data()
    mnist_digits = np.concatenate([x_train, x_test], axis=0)
    mnist_digits = np.expand_dims(mnist_digits, -1).astype("float32") / 255
    mnist_digits = mnist_digits.reshape(-1,28*28*1)
    X=mnist_digits

if False:
    mypath='/home/ihlee/rdfsearch/fingerprintdiff/work1/fingers/'
    mypath_replicas='/home/ihlee/rdfsearch/fingerprintdiff/work1/replica_fingers/'
if False:
    afile=open('rdfsearch.in','r')
    for line in afile:
        if len(line.split()) > 2 :
           if line.split()[1] == 'directory' :
               print(line.split()[0])
               mypath0=line.split()[0]
               if mypath0[-1] != '/': 
                  mypath0=mypath0+'/'
    afile.close()
    mypath=mypath0+'fingers/'
    mypath_replicas=mypath0+'replica_fingers/'
if True:
    mypath=getcwd()+'/fingers/'
    mypath_replicas=getcwd()+'/replica_fingers/'
print(mypath)
print(mypath_replicas)

onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
nsamples=len(onlyfiles)
X=np.zeros((nsamples,ngrid))
Y=np.zeros((nsamples))
xgrd=[]
for i in range(nsamples):
     print(onlyfiles[i])
     afile=open(mypath+onlyfiles[i],'r')
     j=0
     for line in afile:
         if len(line.split()) == 2:
             X[i,j]=float(line.split()[1])
             j=j+1
             if i == 0 :
                xgrd.append(float(line.split()[0]))
     afile.close()
xgrd=np.array(xgrd)
elevation=np.amin(X)
for i in range(nsamples):
    X[i,:]=X[i,:]-elevation

vae = VAE(encoder, decoder)
sgd = keras.optimizers.SGD(lr=0.001, decay=1e-6, momentum=0.9, nesterov=True)
vae.compile(optimizer=sgd)
#vae.compile(optimizer=keras.optimizers.Adam())
vae.fit(X, epochs=3000, batch_size=1)

"""
## Display a grid of sampled digits
"""

# display a 2D manifold of the digits
n = 15  # figure with 15x15 digits
# linearly spaced coordinates on the unit square were transformed
# through the inverse CDF (ppf) of the Gaussian to produce values
# of the latent variables z, since the prior of the latent space
# is Gaussian
u_grid = np.dstack(np.meshgrid(np.linspace(0.05, 0.95, n), np.linspace(0.05, 0.95, n)))
z_grid = norm.ppf(u_grid)
x_decoded = decoder.predict(z_grid.reshape(n*n, 2))
x_decoded = x_decoded.reshape(-1,ngrid) 

if True:
   pathlib.Path(mypath_replicas).mkdir(parents=True, exist_ok=True)
   for i0 in range(n*n):
       mypath1=mypath_replicas+str(i0).rjust(5,'0')+'_finplt'
       print(mypath1)
       afile=open(mypath1,"w")
       jj=0
       for k0 in range(9):
           for j in range(200):
               afile.write(str(xgrd[j])+' '+str(x_decoded[i0,jj]+elevation)+'\n')
               jj=jj+1
           afile.write('&\n')
       afile.close()
print(elevation)
