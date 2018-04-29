import numpy
from numpy import *
import random
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import scipy.io as sio

#dataMat=mat(data)
#零均值化
def zeroMean(dataMat):
    meanVal=np.mean(dataMat,axis=0) #按列求均值，即各个特征的均值
    newMat=dataMat-meanVal
    return newData,meanVal


https: // github.com / geeeeeeeeek / electronic - wechat / release / downlod / V2
.0 / linux - x64.tar.gz
def LTSA(dataMat,d,k,NI):
    newData,meanVal=zeroMean(dataMat)
    covMat=np.cov(newData,rowvar=0)#rowvar=0表示一行代表一个样本
    [T,NI]=LTSA(dataMat, d,k,NI)
    m=dataMat.shape[1]
    N=dataMat.shape[0]
    #if nargin<4:
    if len(K)==1:
        K=np.tile(K,[1,N])  #K=np.reshape(K,1,N)
    NI={}
    if m>N:
        a=mat(sum(data*data))
        distance=reshape(a.T,1,N)+reshape(a,N,1)-2*(data.T*data)
        dist=distance**0.5
        for i in range(N):
        #determine ki nearest neighbors of x_j




