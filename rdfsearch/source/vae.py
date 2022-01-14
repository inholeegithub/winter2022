# Written by In-Ho Lee, KRISS, October 24, 2020.
import random
import numpy as np
from keras import backend as K
from keras.layers import Input, Dense, Lambda, Layer, Add, Multiply
from keras.models import Model, Sequential
from sklearn.model_selection import train_test_split
from sklearn.utils import shuffle
from scipy.stats import norm
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator
def make_plot_latent2d(z_test,y_test,kode):
    plt.figure(figsize=(6,6))
    plt.rcParams['axes.linewidth'] = 1.4 #set the value globally
    ax=plt.axes()
#   ax.set_aspect('equal','box')
#   ax.set_xlabel(r'$\alpha$'+' (arb. unit)',fontsize=22)
#   ax.set_ylabel(r'$\beta$'+' (arb. unit)',fontsize=22)
    ax.set_xlabel(r'$z_1$'+' (arb. unit)',fontsize=22)
    ax.set_ylabel(r'$z_2$'+' (arb. unit)',fontsize=22)
#   majorLocator= MultipleLocator(0.01)
#   minorLocator= AutoMinorLocator()
#   majorFormatter= FormatStrFormatter('%4.2f')
#   minorFormatter= FormatStrFormatter('%4.2f')
#   ax.xaxis.set_major_locator(majorLocator)
#   ax.xaxis.set_major_formatter(majorFormatter)
#   ax.xaxis.set_minor_locator(minorLocator)
#   majorLocator= MultipleLocator(0.05)
#   minorLocator= AutoMinorLocator()
#   majorFormatter= FormatStrFormatter('%4.2f')
#   minorFormatter= FormatStrFormatter('%4.2f')
#   ax.yaxis.set_major_locator(majorLocator)
#   ax.yaxis.set_major_formatter(majorFormatter)
#   ax.yaxis.set_minor_locator(minorLocator)
#   ax.tick_params(which='major', length=2, color='black')
#   ax.tick_params(which='minor', length=4, color='brown')
    ax.set_facecolor("ivory")
#   ax.set_facecolor("beige")
    plt.grid(True)
    plt.scatter(z_test[:,0], z_test[:,1], c=y_test, alpha=.9, s=4**2, cmap=plt.cm.cool)
#   plt.scatter(z_test[:,0], z_test[:,1], c=y_test, alpha=.9, s=4**2, cmap='viridis')
#   plt.scatter(z_test[:,0], z_test[:,1], c=y_test, alpha=.9, s=4**2, cmap=plt.cm.Greens)
#   plt.scatter(z_test[:,0], z_test[:,1], c=y_test, alpha=.9, s=4**2, cmap=plt.cm.binary)
#   plt.scatter(z_test[:,0], z_test[:,1], c=y_test, alpha=.9, s=4**2, cmap=plt.cm.autumn)
#   plt.scatter(z_test[:,0], z_test[:,1], c=y_test, alpha=.9, s=4**2, cmap=plt.cm.RdGy)
    plt.colorbar()
    plt.tight_layout()
    if kode == 0:
       str1='latent2d_test.eps'
    if kode == 1:
       str1='latent2d_train.eps'
    if kode == 2:
       str1='latent2d_train_test.eps'
    plt.savefig(str1,dpi=1200)
#   plt.show()
    plt.close()
def nll(y_true, y_pred):
    """ Negative log likelihood (Bernoulli). """
    # keras.losses.binary_crossentropy gives the mean over the last axis. we require the sum
    return K.sum(K.binary_crossentropy(y_true, y_pred), axis=-1)
class KLDivergenceLayer(Layer):
    """ Identity transform layer that adds KL divergence to the final model loss. """
    def __init__(self, *args, **kwargs):
        self.is_placeholder = True
        super(KLDivergenceLayer, self).__init__(*args, **kwargs)
    def call(self, inputs):
        mu, log_var = inputs
        kl_batch = - .5 * K.sum(1+log_var-K.square(mu)-K.exp(log_var),axis=-1)
        self.add_loss(K.mean(kl_batch), inputs=inputs)
        return inputs
original_dim = 200*9
afile=open('targets',"r")
iline0=0
for line in afile:
    iline0=iline0+1
    if iline0 == 1:
       nsamples1=int(line.split()[0])
       if nsamples1 <  0 :
          nsamples1=-nsamples1
       break
afile.close()
print(nsamples1)
X=np.zeros((nsamples1,original_dim))
Y=np.zeros((nsamples1,))
afile=open('targets',"r")
iline0=0
for line in afile:
    iline0=iline0+1
    if iline0 == 1:
       iline=0
    else :
       fname=line.split()[0]
       fname=fname.strip()+'_fin'
       fname=fname.strip()
       jline=-2
       bfile=open(fname,"r")
       for line11 in bfile:
           jline=jline+1
           if jline > -1:
              X[iline,jline]=float(line11.split()[0])
#             print(X[iline,jline])
       bfile.close()
       Y[iline]=iline
       iline=iline+1
afile.close()
intermediate_dim = 1000
intermediate2_dim = 100
intermediate3_dim = 20
latent_dim = 2
nbatch = nsamples1
if nbatch > 100:
   nbatch = 100
decoder = Sequential([ Dense(intermediate3_dim, input_dim=latent_dim, activation='relu'),
    Dense(intermediate2_dim, activation='relu'),
    Dense(intermediate_dim, activation='relu'),
    Dense(original_dim, activation='sigmoid') ])
x = Input(shape=(original_dim,))
h = Dense(intermediate_dim, activation='relu')(x)
h = Dense(intermediate2_dim, activation='relu')(x)
h = Dense(intermediate3_dim, activation='relu')(x)
z_mu = Dense(latent_dim)(h)
z_log_var = Dense(latent_dim)(h)
z_mu, z_log_var = KLDivergenceLayer()([z_mu, z_log_var])
z_sigma = Lambda(lambda t: K.exp(.5*t))(z_log_var)
eps = Input(tensor=K.random_normal(stddev=1.0,shape=(K.shape(x)[0], latent_dim)))
z_eps = Multiply()([z_sigma, eps])
z = Add()([z_mu, z_eps])
x_pred = decoder(z)
vae = Model(inputs=[x, eps], outputs=x_pred)
vae.summary()
if False:
    trainable_count = int(np.sum([K.count_params(p) for p in set(vae.trainable_weights)]))
    print(trainable_count)
vae.compile(optimizer='adam', loss=nll)
if True:
   ifile=random.randint(0,1000)
   X, Y = shuffle(X,Y,random_state=ifile)
x_train, x_test, y_train, y_test = train_test_split(X, Y, test_size=0.20)
vae.fit(x_train,x_train,shuffle=True,epochs=3000,batch_size=nbatch,validation_data=(x_test,x_test))
encoder = Model(x, z_mu)
print(x_train.shape,x_test.shape)
print(nsamples1,original_dim,intermediate_dim,intermediate2_dim,intermediate3_dim,latent_dim)
# display a 2D plot of the crystal classes in the latent space
if True:
   z_test = encoder.predict(x_test, batch_size=nbatch)
   print(x_test.shape,y_test.shape,z_test.shape)
   print('latent-space, z_test, y_test',len(y_test))
   for i in range(len(y_test)):
       print(z_test[i,0],z_test[i,1],y_test[i])
   make_plot_latent2d(z_test,y_test,0)
if True:
   z_test = encoder.predict(x_train, batch_size=nbatch)
   print(x_train.shape,y_train.shape,z_test.shape)
   print('latent-space, z_test, y_train',len(y_train))
   for i in range(len(y_train)):
       print(z_test[i,0],z_test[i,1],y_train[i])
   make_plot_latent2d(z_test,y_train,1)
if True:
   z_test = encoder.predict(X, batch_size=nbatch)
   print(X.shape,Y.shape,z_test.shape)
   print('latent-space, z_test, Y',len(Y))
   for i in range(len(Y)):
       print(z_test[i,0],z_test[i,1],Y[i])
   make_plot_latent2d(z_test,Y,2)
n = 15  
# linearly spaced coordinates on the unit square were transformed
# through the inverse CDF (ppf) of the Gaussian to produce values of the latent variables z, 
# since the prior of the latent space is Gaussian
if latent_dim == 2:
   u_grid = np.dstack(np.meshgrid(np.linspace(0.05, 0.95, n),np.linspace(0.05, 0.95, n)))
if False and latent_dim == 2:
   u_grid=np.zeros((n*n,latent_dim))
   for ireplica in range(n*n):
       for j in range(latent_dim):
           u_grid[ireplica,j]=random.random()*0.9+0.05
z_grid = norm.ppf(u_grid)
x_decoded = decoder.predict(z_grid.reshape(n*n, latent_dim))
x_decoded = x_decoded.reshape(-1,original_dim)
print(x_decoded.shape,x_decoded[0,:].shape)
print('the number of known structures considered',nsamples1)
ifile=1
for i in range(nsamples1):
    fname=str(ifile).zfill(5)+'_fin' ; fname=fname.strip()
    afile=open(fname,"w")
    afile.write('# 9 1 200 \n') 
    for j in range(original_dim):
        afile.write(str(X[i,j])+'\n')
    afile.close()
    ifile=ifile+1
for i in range(n*n):
    fname=str(ifile).zfill(5)+'_fin' ; fname=fname.strip()
    afile=open(fname,"w")
    afile.write('# 9 1 200 \n') 
    for j in range(original_dim):
        afile.write(str(x_decoded[i,j])+'\n')
    afile.close()
    ifile=ifile+1
