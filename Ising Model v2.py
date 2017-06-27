
# coding: utf-8

# In[1]:

import pandas as pd
import numpy as np
import keras
from keras.models import Sequential
from keras.layers import Dense
from keras.optimizers import SGD
from keras import regularizers
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')


# In[2]:

def read_data(L):
    """
    :type L: int(system size)
    :rtype train,test: pd.dataframe(with phases labeled)
    """
    train_filename='configuration'+str(L)+'.csv'
    test_filename='test'+str(L)+'.csv'
    columns=['temperature']+['Spin'+str(i) for i in xrange(1,L*L+1)]
    train=pd.read_csv('Ising_data/'+'L'+str(L)+'/'+train_filename,names=columns)
    test=pd.read_csv('Ising_data/'+'L'+str(L)+'/'+test_filename,names=columns)
    
    # decide critical temperature
    filename='Cv'+str(L)+'.csv'
    specific_heat=pd.read_csv('Ising_data/'+'L'+str(L)+'/'+filename,names=['temperature','Cv'])
    Tc=specific_heat['temperature'][np.argmax(specific_heat['Cv'])]
    
    # add phase column
    train['FM']=[int(T<=Tc) for T in train['temperature']]
    test['FM']=[int(T<=Tc) for T in test['temperature']]
    train['PM']=[int(T>Tc) for T in train['temperature']]
    test['PM']=[int(T>Tc) for T in test['temperature']]
    
    return train, test


# In[3]:

train, test=read_data(60)


# In[4]:

def data_process(train, test):
    """
    :type train,test: pd.dataframe(with phases labeled)
    :rtype trX, trY, teX, teY
    """
    # shuffle
    train=train.sample(frac=1).reset_index(drop=True)
    
    trX, trY = train.drop(['FM','PM','temperature'],axis=1), np.array(train[['FM','PM']])
    teX, teY = test.drop(['FM','PM','temperature'],axis=1), np.array(test[['FM','PM']])
    
    return trX, trY, teX, teY


# In[5]:

trX, trY, teX, teY=data_process(train, test)


# In[6]:

def create_model(hidden_units=3, lamb=0.05,lr=0.001):
    classifier=Sequential()
    classifier.add(Dense(output_dim=hidden_units, init='uniform', kernel_regularizer=regularizers.l2(lamb),                      activation='sigmoid',input_dim=trX.shape[1]))
    classifier.add(Dense(output_dim=1, init='uniform', kernel_regularizer=regularizers.l2(lamb),                      activation='sigmoid'))
    sgd = SGD(lr=0.001, decay=1e-6)#, momentum=0.9, nesterov=True)
    classifier.compile(optimizer='adam',loss='binary_crossentropy',metrics=['accuracy'])
    return classifier


# In[7]:

classifier=create_model(hidden_units=3, lamb=0.0)


# In[ ]:

classifier.fit(np.array(trX),np.argmax(trY,axis=1),batch_size=100, nb_epoch=1,verbose=1)


# In[10]:

y_pred=classifier.predict(np.array(teX))


# In[11]:

np.mean((y_pred>0.5)[:,0]==np.argmax(teY,axis=1))


# In[12]:

Temperatures=np.array(sorted(list(set(test['temperature']))))
Pred_avgs=[]
num_T=Temperatures.shape[0]
num_test=len(test['temperature'])
batch=num_test/num_T

for i in xrange(num_T):
    start=batch*i
    end=start+batch
    batch_prediction=y_pred[start:end]
    Pred_avgs.append(np.mean(batch_prediction,axis=0))

Pred_avgs=np.array(Pred_avgs)


# In[13]:

import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')

plt.plot(Temperatures, Pred_avgs[:,0], 'bo', Temperatures, Pred_avgs[:,0], 'b-',          Temperatures, 1-Pred_avgs[:,0], 'ro', Temperatures, 1-Pred_avgs[:,0], 'r-')


# In[ ]:



