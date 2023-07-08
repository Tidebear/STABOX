import numpy as np
import scanpy as sc
import pandas as pd

from annoy import AnnoyIndex
import random

from scipy.sparse import issparse
from scipy.spatial import cKDTree

from sklearn.base import BaseEstimator
from sklearn.neighbors import NearestNeighbors
import json
import os
import shutil
import multiprocessing
import platform

from sklearn.metrics.pairwise import rbf_kernel, euclidean_distances
from sklearn.neighbors import NearestNeighbors
from intervaltree import IntervalTree
import operator

from scipy.sparse import issparse
from annoy import AnnoyIndex
from multiprocessing import Process, cpu_count, Queue
from collections import namedtuple
from operator import attrgetter
from tqdm import tqdm
import time
import itertools
import networkx as nx

import hnswlib

from sklearn.preprocessing import LabelEncoder
le = LabelEncoder()
from sklearn import metrics


def generator_from_index(adata, batch_name,  celltype_name=None, mask_batch=None, Y = None, k = 20, label_ratio = 0.8, k_to_m_ratio = 0.75, batch_size = 32, search_k=-1,
                         save_on_disk = True, approx = True, verbose=1):

    print('version 0.0.2. 09:00, 12/01/2020')

    # Calculate MNNs by pairwise comparison between batches
    
    cells = adata.obs_names
    
    if(verbose > 0):
        print("Calculating MNNs...")
           
    mnn_dict = create_dictionary_mnn(adata, batch_name=batch_name, k = k, save_on_disk = save_on_disk, approx = approx, verbose = verbose)

    if(verbose > 0):
        print(str(len(mnn_dict)) + " cells defined as MNNs")
        
    if celltype_name is None:
        label_dict=dict()
    else:
        
        if (verbose > 0):
            print ('Generating supervised positive pairs...')

        label_dict_original = create_dictionary_label(adata, celltype_name= celltype_name, batch_name=batch_name,  mask_batch=mask_batch,  k=k, verbose=verbose)
        num_label = round(label_ratio * len(label_dict_original))

        cells_for_label = np.random.choice(list(label_dict_original.keys()), num_label, replace = False)

        label_dict = {key: value for key, value in label_dict_original.items() if key in cells_for_label}

        if(verbose > 0):
            print(str(len(label_dict.keys())) + " cells defined as supervision triplets")

        print (len(set(mnn_dict.keys())&set(label_dict.keys())))

    if k_to_m_ratio == 0.0:
        knn_dict = dict()
    else:
        num_k = round(k_to_m_ratio * len(mnn_dict))
        # Calculate KNNs for subset of residual cells
        # 除MNN节点之外的所有节点
        cells_for_knn = list(set(cells) - (set(list(label_dict.keys())) | set(list(mnn_dict.keys()))))
        if(len(cells_for_knn) > num_k): # 如果剩余的节点过多，就只选择num_k个
            cells_for_knn = np.random.choice(cells_for_knn, num_k, replace = False)

        if(verbose > 0):
            print("Calculating KNNs...")

        cdata = adata[cells_for_knn]
        knn_dict = create_dictionary_knn(cdata, cells_for_knn, k = k, save_on_disk = save_on_disk, approx = approx)
        if(verbose > 0):
            print(str(len(cells_for_knn)) + " cells defined as KNNs")

    final_dict = merge_dict(mnn_dict, label_dict)
    final_dict.update(knn_dict)

    
    cells_for_train = list(final_dict.keys())
    print ('Total cells for training:'+ str(len(cells_for_train)))

    ddata = adata[cells_for_train]

    # Reorder triplet list according to cells
    if(verbose > 0):
        print("Reorder")
    names_as_dict = dict(zip(list(adata.obs_names), range(0, adata.shape[0]))) #建立adata.obs_names与顺序数字编号的对应关系
    def get_indices2(name):
          return([names_as_dict[x] for x in final_dict[name]])

    triplet_list = list(map(get_indices2, cells_for_train)) #把用于训练的细胞单独提出来到这个list，细胞重新顺序编号，用于找anchor
    
    batch_list = ddata.obs[batch_name]  #用于训练的ddata的obs_names还是原始的编号
    batch_indices = []
    for i in batch_list.unique(): #把三个批次的细胞分别提取出来
        batch_indices.append(list(np.where(batch_list == i)[0])) #但是这里的编号就被重新设定了

    batch_as_dict = dict(zip(list(batch_list.unique()), range(0, len(batch_list.unique()))))
    tmp = map(lambda _: batch_as_dict[_], batch_list)
    batch_list = list(tmp)

    if Y is None:
        return KnnTripletGenerator(X = ddata.obsm["X_pca"], X1 = adata.obsm['X_pca'], dictionary = triplet_list,
                               batch_list = batch_list, batch_indices = batch_indices, batch_size=batch_size)

    else:
        tmp = dict(zip(cells, Y))
        Y_new = [tmp[x] for x in cells_for_train]
        Y_new = le.fit_transform(Y_new)
        return LabeledKnnTripletGenerator(X = ddata.obsm["X_pca"], X1 = adata.obsm['X_pca'], Y = Y_new, dictionary = triplet_list,
                               batch_list = batch_list, batch_indices = batch_indices, batch_size = batch_size)


def merge_dict(x,y):
    for k,v in x.items():
                if k in y.keys():
                    y[k] += v
                else:
                    y[k] = v
    return y


def create_dictionary_mnn(adata, use_rep, batch_name, k = 50, save_on_disk = True, approx = True, verbose = 1, iter_comb = None):

    cell_names = adata.obs_names

    batch_list = adata.obs[batch_name]
    datasets = []
    datasets_pcs = []
    cells = []
    for i in batch_list.unique():
        datasets.append(adata[batch_list == i])
        datasets_pcs.append(adata[batch_list == i].obsm[use_rep])
        cells.append(cell_names[batch_list == i])

    batch_name_df = pd.DataFrame(np.array(batch_list.unique()))
    mnns = dict()

    # if len(cells) > 2:
    if iter_comb is None:
        iter_comb = list(itertools.combinations(range(len(cells)), 2))
    # for comb in list(itertools.combinations(range(len(cells)), 2)): # 返回多个批次所有可能的组合
    for comb in iter_comb:
        i = comb[0]
        j = comb[1]
        key_name1 = batch_name_df.loc[comb[0]].values[0] + "_" + batch_name_df.loc[comb[1]].values[0]
        mnns[key_name1] = {}

        if(verbose > 0):
            print('Processing datasets {}'.format((i, j)))

        new = list(cells[j])
        ref = list(cells[i])

        ds1 = adata[new].obsm[use_rep]
        ds2 = adata[ref].obsm[use_rep]
        names1 = new
        names2 = ref
        # 如果K>1，则MNN点就很有可能出现1:n,即一对多的情况
        match = mnn(ds1, ds2, names1, names2, knn=k, save_on_disk = save_on_disk, approx = approx)

        G = nx.Graph()
        G.add_edges_from(match)
        node_names = np.array(G.nodes)
        anchors = list(node_names)
        adj = nx.adjacency_matrix(G) #src 和 dst 中的points拼起来作为矩阵的列或行，这个矩阵是对称的
        #https://blog.csdn.net/Snowmyth/article/details/121280577?spm=1001.2101.3001.6661.1&utm_medium=distribute.pc_relevant_t0.none-task-blog-2%7Edefault%7ECTRLIST%7ERate-1.pc_relevant_paycolumn_v3&depth_1-utm_source=distribute.pc_relevant_t0.none-task-blog-2%7Edefault%7ECTRLIST%7ERate-1.pc_relevant_paycolumn_v3&utm_relevant_index=1
        #indptr提示的是非零数在稀疏矩阵中的位置信息。indices是具体的连接边的一个节点的编号。
        # https://www.csdn.net/tags/NtzaQgysMjgyNjEtYmxvZwO0O0OO0O0O.html
        tmp = np.split(adj.indices, adj.indptr[1:-1]) #把一个数组从左到右按顺序切分;

        for i in range(0, len(anchors)): #把src 和 dst 中所有的points都包含到字典的key中了
            key = anchors[i]
            i = tmp[i]
            names = list(node_names[i])
            # mnns这里是个字典，多个切片时，由于key是相同的
            # 最后一个切片的mnn会把前面的mnn都覆盖掉，导致最后一个切片的mnn特别多！
            mnns[key_name1][key]= names
    return(mnns)


def create_dictionary_knn(adata, use_rep, cell_subset, k = 50, save_on_disk = True, approx = True):

    # cell_names = adata.obs_names

    dataset = adata[cell_subset]
    pcs = dataset.obsm[use_rep]

    def get_names(ind):
        return np.array(cell_subset)[ind]

    if approx:
        dim = pcs.shape[1]
        num_elements = pcs.shape[0]
        p = hnswlib.Index(space='l2', dim = dim)
        p.init_index(max_elements=num_elements, ef_construction=100, M=16)
        p.set_ef(10)
        p.add_items(pcs)
        ind, distances = p.knn_query(pcs, k=k)
        ind = ind[1:] # remove self-point

        cell_subset = np.array(cell_subset)
        names = list(map(lambda x: cell_subset[x], ind))
        knns = dict(zip(cell_subset, names))

    else:
        nn_ = NearestNeighbors(n_neighbors = k, p = 2)
        nn_.fit(pcs)
        ind = nn_.kneighbors(pcs, return_distance=False)
        ind = ind[1:]  # remove self-point

        names = list(map(lambda x: cell_subset[x], ind))
        knns = dict(zip(cell_subset, names))

    return(knns)


def validate_sparse_labels(Y):
    if not zero_indexed(Y):
        raise ValueError('Ensure that your labels are zero-indexed')
    if not consecutive_indexed(Y):
        raise ValueError('Ensure that your labels are indexed consecutively')


def zero_indexed(Y):
    if min(abs(Y)) != 0:
        return False
    return True


def consecutive_indexed(Y):
    """ Assumes that Y is zero-indexed. """
    n_classes = len(np.unique(Y[Y != np.array(-1)]))
    if max(Y) >= n_classes:
        return False
    return True


def nn_approx(ds1, ds2, names1, names2, knn=50, pos_knn=None):
    dim = ds2.shape[1]
    num_elements = ds2.shape[0]
    p = hnswlib.Index(space='l2', dim=dim)
    p.init_index(max_elements=num_elements, ef_construction=100, M = 16)
    p.set_ef(10)
    p.add_items(ds2)
    ind,  distances = p.knn_query(ds1, k=knn)
    ## 遍历ind中每一行的每个KNN节点，判断当前点的空间相近KNN和跨批次的KNN有多少是重叠的, 去掉重叠过少的点
    match = set()
    for a, b in zip(range(ds1.shape[0]), ind):
        # knn_used = np.asarray([len(set(pos_knn[b[i]]).intersection(set(b))) for i in range(knn)]) > 0
        # for b_i in b[knn_used]:
        #     match.add((names1[a], names2[b_i]))
        # for b_i in b:
        #     if len(set(pos_knn[b_i]).intersection(set(b))) > 0: #去掉没有重叠的点
        #         match.add((names1[a], names2[b_i]))
        for b_i in b:
            match.add((names1[a], names2[b_i]))
    return match


def nn(ds1, ds2, names1, names2, knn=50, metric_p=2):
    # Find nearest neighbors of first dataset.
    nn_ = NearestNeighbors(knn, p=metric_p)
    nn_.fit(ds2)
    ind = nn_.kneighbors(ds1, return_distance=False)

    match = set()
    for a, b in zip(range(ds1.shape[0]), ind):
        for b_i in b:
            match.add((names1[a], names2[b_i]))

    return match


def nn_annoy(ds1, ds2, names1, names2, knn = 20, metric='euclidean', n_trees = 50, save_on_disk = True):
    """ Assumes that Y is zero-indexed. """
    # Build index.
    a = AnnoyIndex(ds2.shape[1], metric=metric)
    if(save_on_disk):
        a.on_disk_build('annoy.index')
    for i in range(ds2.shape[0]):
        a.add_item(i, ds2[i, :])
    a.build(n_trees)

    # Search index.
    ind = []
    for i in range(ds1.shape[0]):
        ind.append(a.get_nns_by_vector(ds1[i, :], knn, search_k=-1))
    ind = np.array(ind)

    # Match.
    match = set()
    for a, b in zip(range(ds1.shape[0]), ind):
        for b_i in b:
            match.add((names1[a], names2[b_i]))

    return match


def mnn(ds1, ds2, names1, names2, knn = 20, save_on_disk = True, approx = True, pos_knn1=None, pos_knn2=None):
    # Find nearest neighbors in first direction.
    if approx: #输出KNN pair; match1: (names1中节点，names2中节点), 大小为ds1.shape[0]*knn
        match1 = nn_approx(ds1, ds2, names1, names2, knn=knn, pos_knn=pos_knn1)#, save_on_disk = save_on_disk)
        # Find nearest neighbors in second direction.
        match2 = nn_approx(ds2, ds1, names2, names1, knn=knn, pos_knn=pos_knn2)#, save_on_disk = save_on_disk)
    else:
        match1 = nn(ds1, ds2, names1, names2, knn=knn)
        match2 = nn(ds2, ds1, names2, names1, knn=knn)
    # Compute mutual nearest neighbors.
    mutual = match1 & set([ (b, a) for a, b in match2 ])

    return mutual
