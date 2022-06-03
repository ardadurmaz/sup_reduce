
# coding: utf-8

# In[10]:


import numpy as np
import pandas as pd
from sklearn.cross_decomposition import PLSRegression
import matplotlib.pyplot as plt
from matplotlib import cm
import proplot as plot
from sklearn.metrics import pairwise_distances
import seaborn as sns


# In[2]:


in_wd = '/home/axd497/wgs_ecoli/'
mut_lat = np.loadtxt('{}data/MutationLatentCoordinates.tsv'.format(in_wd), delimiter='\t')
mic_lat = np.loadtxt('{}data/MICLatentCoordinates.tsv'.format(in_wd), delimiter='\t')
mic_mat = np.loadtxt('{}data/MICMat.tsv'.format(in_wd), delimiter='\t')


# In[ ]:


# LOO Cross-validation for number of dimension selection
# Sample Code from another data++
pred_err_cv = []
exp_err_cv = []

kf = KFold(n_splits=5, shuffle=True)
for train_idx, test_idx in kf.split(expr_mat):
    x_train, x_test, y_train, y_test = [expr_mat[train_idx,], expr_mat[test_idx,], auc_mat[train_idx,], auc_mat[test_idx,]]
    local_pred_err_cv = []
    local_exp_err_cv = []
    for k in np.arange(1, 11):
        local_fit = PLSRegression(n_components=k).fit(x_train, y_train)
        local_pred_err_cv.append(np.sum((y_test-local_fit.predict(x_test))**2))
        local_exp_err_cv.append(np.sum((x_test-local_fit.inverse_transform(local_fit.transform(x_test)))**2))
    pred_err_cv.append(np.asarray(local_pred_err_cv))
    exp_err_cv.append(np.asarray(local_exp_err_cv))


# In[ ]:


# PLS
pls_fit = PLSRegression(n_components=3).fit(mut_lat, mic_mat)
pls_embedd = pls_fit.transform(mut_lat)
pls_pred = pls_fit.predict(mut_lat)


# In[12]:


cos_sim=1.0-pairwise_distances(pls_fit.y_loadings_, pls_fit.y_loadings_, metric='cosine')
np.savetxt('{}data/CosineSimilarity.tsv'.format(in_wd), X=cos_sim, delimiter='\t')
#sns.clustermap(cos_sim, cmap='vlag', xaxisticklabels=)


# In[3]:


# PLS
pls_fit = PLSRegression(n_components=2).fit(mut_lat, mic_mat)
pls_embedd = pls_fit.transform(mut_lat)
pls_pred = pls_fit.predict(mut_lat)


# In[34]:


X, Y = np.meshgrid(np.linspace(-3.0, 3.0, 100), np.linspace(-3.0, 3.0, 100))
cont_pred = np.zeros(shape=[X.shape[0], X.shape[1], mic_mat.shape[1]])
for i in range(X.shape[0]):
    for j in range(X.shape[1]):
        cont_pred[i,j,] = np.matmul(np.asarray([X[i,j], Y[i,j]]).reshape([1, -1]), pls_fit.y_loadings_.T)
cont_pred


# In[40]:


drug_ids = pd.read_csv('{}data/FeatAnnot.tsv'.format(in_wd))
drug_ids


# In[44]:


fig, ax = plot.subplots(nrows=1, ncols=9)
for i in range(9):
    ax[0,i].contourf(X, Y, cont_pred[::,::,i], colorbar='r')
    ax[0,i].format(title=drug_ids.loc[i,'id'])
ax.format(xlabel='PLS Component 1', ylabel='PLS Component 2', suptitle='PLS Embedding')
fig.savefig('{}/PLS_Embedding.pdf'.format(in_wd))
#plot.show()


# In[8]:


fig, ax = plot.subplots(nrows=1, ncols=1)
ax.scatter(pls_embedd[::,0], pls_embedd[::,1])
ax[0,0].format(xlabel='PLS Component 1', ylabel='PLS Component 2', suptitle='PLS Embedding')
plot.show()

