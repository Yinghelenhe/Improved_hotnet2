#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 17:34:04 2018

@author: YINGHELENHE
"""
"""This script contains all the function I define"""

import sys, os, json, h5py, numpy as np, scipy.io, networkx as nx, scipy as sp
from collections import defaultdict
from constants import *
from hotnet2 import component_sizes
from scipy.stats import norm
import hnio
import multiprocessing as mp
import hotnet2 as hn
import math


###############################################################################################
###### Output gene score for each subnetwork
def output_scoreVec(scoreVec,genelist, output_dir):
    output_dir = '{}/consensus/New_Results/'.format(output_dir)
    hnio.setup_output_dir(output_dir)

    with open(output_dir + '/scoreVec.tsv', 'w') as f:f.write('\n'.join('{}\t{}'.format(gene, score) for gene, score in zip(genelist, scoreVec)))

##### output subnetwork with score only
def output_scoreSubnet(score_subnet,output_dir):
    output_dir = '{}/consensus/New_Results/'.format(output_dir)
    hnio.setup_output_dir(output_dir)

    with open(output_dir + '/Score_Subnet.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in score_subnet]))

#####output permutated subnetwork, with score
def output_permuted_subnetscore(result,output_dir):
    for subscore, delta in result:
        # create output directory
        delta_out_dir = os.path.abspath(output_dir + "/delta_" + str(delta)+ "/New_Results" + "/Permuted_Results")
        if not os.path.isdir(delta_out_dir): os.mkdir(delta_out_dir)
        
        with open(os.path.abspath(delta_out_dir) + '/Permuted_Subnet_Score.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in subscore]))

def output_subnetscore_deltas(result,output_dir):
    for subscore, delta in result:
        # create output directory
        delta_out_dir = os.path.abspath(output_dir + "/delta_" + str(delta)+ "/New_Results")
        if not os.path.isdir(delta_out_dir): os.mkdir(delta_out_dir)
        
        with open(os.path.abspath(delta_out_dir) + '/Subnet_Score.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in subscore]))


def output_subnet_sig_deltas(result,output_dir):
    for subscore, delta in result:
        # create output directory
        delta_out_dir = os.path.abspath(output_dir + "/delta_" + str(delta)+ "/New_Results")
        if not os.path.isdir(delta_out_dir): os.mkdir(delta_out_dir)
        
        with open(os.path.abspath(delta_out_dir) + '/Subnetscore_with_Significant.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in subscore]))


def output_component_sig_deltas(result,output_dir):
    for subscore, delta in result:
        # create output directory
        delta_out_dir = os.path.abspath(output_dir + "/delta_" + str(delta)+ "/New_Results")
        if not os.path.isdir(delta_out_dir): os.mkdir(delta_out_dir)
        
        with open(os.path.abspath(delta_out_dir) + '/Components_with_Significant.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in subscore]))


def output_alpha_sig_deltas(result,output_dir):
    for subscore, delta in result:
        # create output directory
        delta_out_dir = os.path.abspath(output_dir + "/delta_" + str(delta)+ "/New_Results")
        if not os.path.isdir(delta_out_dir): os.mkdir(delta_out_dir)
        
        with open(os.path.abspath(delta_out_dir) + '/Alpha_Significant.tsv', 'w') as f:f.write(str(subscore))


def output_significant_component(result,output_dir):
    for subscore, delta in result:
        # create output directory
        delta_out_dir = os.path.abspath(output_dir + "/delta_" + str(delta) + "/New_Results")
        if not os.path.isdir(delta_out_dir): os.mkdir(delta_out_dir)
        
        with open(os.path.abspath(delta_out_dir) + '/Significant_Components.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in subscore]))


def output_significant_count(result,output_dir):
    for subscore, delta in result:
        # create output directory
        delta_out_dir = os.path.abspath(output_dir + "/delta_" + str(delta)+ "/New_Results")
        if not os.path.isdir(delta_out_dir): os.mkdir(delta_out_dir)
        
        with open(os.path.abspath(delta_out_dir) + '/Significant_Counts.tsv', 'w') as f:f.write(str(subscore))



def output_gene_degree(result,output_dir):
    for subscore, delta in result:
        # create output directory
        delta_out_dir = os.path.abspath(output_dir + "/delta_" + str(delta) + "/New_Results")
        if not os.path.isdir(delta_out_dir): os.mkdir(delta_out_dir)
        
        with open(os.path.abspath(delta_out_dir) + '/Gene_degree.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in subscore]))



def output_permuted_gene_degree(result,output_dir):
    for subscore, delta in result:
        # create output directory
        delta_out_dir = os.path.abspath(output_dir + "/delta_" + str(delta) + "/New_Results" + "/Permuted_Results")
        if not os.path.isdir(delta_out_dir): os.mkdir(delta_out_dir)
        
        with open(os.path.abspath(delta_out_dir) + '/Permuted_Gene_degree.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in subscore]))

def output_Degree_Weighted_Score(result,output_dir):
    for subscore, delta in result:
        # create output directory
        delta_out_dir = os.path.abspath(output_dir + "/delta_" + str(delta) + "/New_Results")
        if not os.path.isdir(delta_out_dir): os.mkdir(delta_out_dir)
        
        with open(os.path.abspath(delta_out_dir) + '/Degree_Weighted_Score.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in subscore]))


def output_Permuted_Degree_Weighted_Score(result,output_dir):
    for subscore, delta in result:
        # create output directory
        delta_out_dir = os.path.abspath(output_dir + "/delta_" + str(delta) + "/New_Results" + "/Permuted_Results")
        if not os.path.isdir(delta_out_dir): os.mkdir(delta_out_dir)
        
        with open(os.path.abspath(delta_out_dir) + '/Permuted_Degree_Weighted_Score.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in subscore]))


def output_Centrality_Weighted_Score(result,output_dir):
    for subscore, delta in result:
        # create output directory
        delta_out_dir = os.path.abspath(output_dir + "/delta_" + str(delta) + "/New_Results")
        if not os.path.isdir(delta_out_dir): os.mkdir(delta_out_dir)
        
        with open(os.path.abspath(delta_out_dir) + '/Centrality_Weighted_Score.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in subscore]))


def output_Permuted_Centrality_Weighted_Score(result,output_dir):
    for subscore, delta in result:
        # create output directory
        delta_out_dir = os.path.abspath(output_dir + "/delta_" + str(delta) + "/New_Results" + "/Permuted_Results")
        if not os.path.isdir(delta_out_dir): os.mkdir(delta_out_dir)
        
        with open(os.path.abspath(delta_out_dir) + '/Permuted_Centrality_Weighted_Score.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in subscore]))


def output_Myresults(result,output_dir):
    
    for a,b,c,d,e,count,g,h,i,j,k,l,m,n_ss,o_ss, p_ss,q_ss,n_dwss, o_dwss, p_dwss, q_dwss, n_cwss, o_cwss, p_cwss, q_cwss,n_dwsd, o_dwsd, p_dwsd, q_dwsd, n_cwsd, o_cwsd, p_cwsd, q_cwsd, delta in result:
        
        delta_out_dir1 = os.path.abspath(output_dir + "/delta_" + str(delta) + "/New_Results" )
        if not os.path.isdir(delta_out_dir1): os.mkdir(delta_out_dir1)
        

        delta_out_dir2 = os.path.abspath(output_dir + "/delta_" + str(delta) + "/New_Results" + "/Permuted_Results")
        if not os.path.isdir(delta_out_dir2): os.mkdir(delta_out_dir2)
        
        
        delta_out_dir3 = os.path.abspath(output_dir + "/delta_" + str(delta) + "/New_Results" + "/Scroe_functions")
        if not os.path.isdir(delta_out_dir3): os.mkdir(delta_out_dir3)
        
        delta_out_dir4 = os.path.abspath(output_dir + "/delta_" + str(delta) + "/New_Results" + "/Sum_Size")
        if not os.path.isdir(delta_out_dir4): os.mkdir(delta_out_dir4)
        
        delta_out_dir5 = os.path.abspath(output_dir + "/delta_" + str(delta) + "/New_Results" + "/Sum_Degree")
        if not os.path.isdir(delta_out_dir5): os.mkdir(delta_out_dir5)
        
        delta_out_dir6 = os.path.abspath(output_dir + "/delta_" + str(delta) + "/New_Results" + "/DegreeWSum_Size")
        if not os.path.isdir(delta_out_dir6): os.mkdir(delta_out_dir6)
        
        delta_out_dir7 = os.path.abspath(output_dir + "/delta_" + str(delta) + "/New_Results" + "/CentWSum_Size")
        if not os.path.isdir(delta_out_dir7): os.mkdir(delta_out_dir7)
        
        delta_out_dir8 = os.path.abspath(output_dir + "/delta_" + str(delta) + "/New_Results" + "/DegreeWSum_Degree")
        if not os.path.isdir(delta_out_dir8): os.mkdir(delta_out_dir8)
        
        delta_out_dir9 = os.path.abspath(output_dir + "/delta_" + str(delta) + "/New_Results" + "/CentWSum_Degree")
        if not os.path.isdir(delta_out_dir9): os.mkdir(delta_out_dir9)
        



        with open(os.path.abspath(delta_out_dir1) + '/Alpha_Significant.tsv', 'w') as f:f.write(str(c))
        
        with open(os.path.abspath(delta_out_dir1) + '/Gene_degree.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in h]))
        
        
        with open(os.path.abspath(delta_out_dir2) + '/Permuted_Subnet_Score.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in a]))
        
        with open(os.path.abspath(delta_out_dir2) + '/Permuted_Gene_degree.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in i]))

        with open(os.path.abspath(delta_out_dir2) + '/Permuted_Centrality_Weighted_Score.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in l]))
    
    
        with open(os.path.abspath(delta_out_dir2) + '/Permuted_Degree_Weighted_Score.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in m]))


        with open(os.path.abspath(delta_out_dir3) + '/Centrality_Weighted_Score.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in j]))
        
        with open(os.path.abspath(delta_out_dir3) + '/Degree_Weighted_Score.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in k]))
        
        with open(os.path.abspath(delta_out_dir3) + '/Subnet_Score.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in b]))
        
        
        #####sum and total number of gene in a subnetwork
        with open(os.path.abspath(delta_out_dir4) + '/Subnetscore_with_Significant.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in d]))
        with open(os.path.abspath(delta_out_dir4) + '/Components_with_Significant.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in e]))
        with open(os.path.abspath(delta_out_dir4) + '/Significant_Counts.tsv', 'w') as f:f.write(str(count))
        with open(os.path.abspath(delta_out_dir4) + '/Significant_Components.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in g]))
        
        ####sum and sum of degree
        with open(os.path.abspath(delta_out_dir5) + '/Components_with_Significant.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in o_ss]))
        with open(os.path.abspath(delta_out_dir5) + '/Significant_Counts.tsv', 'w') as f:f.write(str(p_ss))
        with open(os.path.abspath(delta_out_dir5) + '/Significant_Components.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in q_ss]))

        ######degree weigthed sum and total number of subnetwork
        with open(os.path.abspath(delta_out_dir6) + '/Components_with_Significant.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in o_dwss]))
        with open(os.path.abspath(delta_out_dir6) + '/Significant_Counts.tsv', 'w') as f:f.write(str(p_dwss))
        with open(os.path.abspath(delta_out_dir6) + '/Significant_Components.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in q_dwss]))


       ######degree weigthed sum and total number of subnetwork
        with open(os.path.abspath(delta_out_dir7) + '/Components_with_Significant.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in o_cwss]))
        with open(os.path.abspath(delta_out_dir7) + '/Significant_Counts.tsv', 'w') as f:f.write(str(p_cwss))
        with open(os.path.abspath(delta_out_dir7) + '/Significant_Components.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in q_cwss]))

      ######degree weigthed sum and sum of degree
        with open(os.path.abspath(delta_out_dir8) + '/Components_with_Significant.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in o_dwsd]))
        with open(os.path.abspath(delta_out_dir8) + '/Significant_Counts.tsv', 'w') as f:f.write(str(p_dwsd))
        with open(os.path.abspath(delta_out_dir8) + '/Significant_Components.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in q_dwsd]))

     ######centrality weigthed sum and sum of degree
        with open(os.path.abspath(delta_out_dir9) + '/Components_with_Significant.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score)) for score in o_cwsd]))
        with open(os.path.abspath(delta_out_dir9) + '/Significant_Counts.tsv', 'w') as f:f.write(str(p_cwsd))
        with open(os.path.abspath(delta_out_dir9) + '/Significant_Components.tsv', 'w') as f:f.write('\n'.join([' '.join(map(str,score))  for score in q_cwsd]))



##############################################################################################
######score vect
def score_vec(infmat, index2gene, gene2heat):
    """Create and return a similarity matrix and index to gene mapping for the given influence
        matrix and heat. Only genes with heat that are in the network will be included in the returned
        similarity matrix and index to gene mapping.
        
        Arguments:
        infmat -- 2D ndarray representing the full influence matrix
        index2gene -- dict mapping an index in the matrix to the name of the gene represented at that
        index in the influence matrix
        gene2heat -- dict mapping a gene name to the heat score for that gene
        directed -- if True, sim[i][j] = inf(i,j)*heat[i] and sim[i][j] != sim[j][i]
        if False, sim[i][j] = min(inf(i,j), inf(j,i))*max(heat(i), heat(j))
        
        """
    start_index = min(index2gene.keys())
    gene2index = dict((gene, index) for index, gene in index2gene.iteritems())
    
    # Identify genes in the given list that are also in the network
    genelist = sorted(set(gene2heat.keys()).intersection(gene2index.keys()))
    index2gene = dict(enumerate(genelist))
    infmat = np.asarray(infmat, dtype=np.float64)
    h = np.array([gene2heat[g] for g in genelist], dtype=np.float64)
    indices = np.array([gene2index[g]-start_index for g in genelist], dtype=np.int)
    ####find influence matrix
    M = infmat[np.ix_(indices, indices)]
    ####final state of genes socre
    scoreVec = np.dot(h,M)
    return genelist, scoreVec

#############################################################################################
def Matchscore_fun(genelist, scoreVec,consensus):
    subnetwork = [sorted(d['core'])+ sorted(d['expansion']) for d in consensus]
    score = np.column_stack((genelist, scoreVec))
    subnet_score= []
    size = np.shape(subnetwork)[0]
    for i in range(size):
        size2 = np.shape(subnetwork[i])[0]
        temp = [1.0]*size2
        for j in range(size2):
            ind = [k for k, x in enumerate(score[:,0] == subnetwork[i][j]) if x][0]
            temp[j]= float(score[ind][1])
        subnet_score.append(temp)
            #print subnet_score
    return subnet_score



def Matchscore_Permuted_fun(genelist, scoreVec, consensus):
    ####converte set to list
    subnetwork = [list(i) for i in consensus]
    score = np.column_stack((genelist, scoreVec))
    subnet_score= []
    size = np.shape(subnetwork)[0]
    for i in range(size):
        size2 = np.shape(subnetwork[i])[0]
        temp = [1.0]*size2
        for j in range(size2):
            ind = [k for k, x in enumerate(score[:,0] == subnetwork[i][j]) if x][0]
            temp[j]= float(score[ind][1])
        subnet_score.append(temp)
    return subnet_score


##############################################################################################
"""Functions to calculate permuted subnetwork score """
    
"""This is the functions to calculate single permuted subnetowrk score"""
def permuted_subnetScore((infmat, index2gene, heat_permutation, delta, verbose, min_cc_size, G)):
        
        
    permuted_sim, permuted_index2gene = hn.similarity_matrix(infmat, index2gene, heat_permutation, True, verbose=verbose)
    permuted_G = hn.weighted_graph(permuted_sim, permuted_index2gene, delta, directed=True)
        
    
    ###permuted_ccs is the permutated subnetwork
    permuted_ccs = hn.connected_components(permuted_G, min_cc_size)
    
    ####permuted gene degree in each subnetwork
    permuted_gene_degree = Gene_Degree_func(G,permuted_ccs)
    
    ####find centrality of gene in each subnetwork
    permuted_gene_centrality = Gene_Centrality_func(G,permuted_ccs)
    
    #### find final gene score for each gene
    permuted_genelist, permuted_scoreVec = score_vec(infmat,permuted_index2gene, heat_permutation)
    
    ####find subnetwork with score only
    permuted_subnetScore = Matchscore_Permuted_fun(permuted_genelist, permuted_scoreVec,permuted_ccs)
    
    #####find gene score that weighted by their own degree
    permuted_degree_weighted_score = Product_func(permuted_subnetScore,permuted_gene_degree)
    
    #####weighted gene by centrality
    permuted_centality_weighted_score  = Product_func(permuted_subnetScore,permuted_gene_centrality)
   
    
    return (permuted_subnetScore,permuted_gene_degree,permuted_degree_weighted_score, permuted_centality_weighted_score)
        
        
        
def calculate_permuted_subnetScore(infmat, index2gene, heat_permutations, delta, verbose, min_cc_size, G, num_cores=1):
    """Return a list of permuted subnetwork score for each permutation.

         Arguments:
        infmat -- 2D ndarray representing an influence matrix
        index2gene -- dict mapping an index in the matrix to the name of the gene represented at that
        index in the influence matrix
        heat_permutations -- iterable of dicts mapping a gene name to the permuted heat score for that gene
        delta -- threshold for edge weight removal
        sizes -- list of sizes for which the number of connected components of that sizes should be
        calculated in each permutation
        num_cores -- number of cores to use for running in parallel
        
        s"""
    
    if num_cores != 1:
          pool = mp.Pool(None if num_cores == -1 else num_cores)
          map_fn = pool.map
    else:
          map_fn = map
            
    args = [(infmat, index2gene, heat_permutation, delta, verbose,min_cc_size, G) for heat_permutation in heat_permutations]
            
    all_permuted_subnetsocre = map_fn(permuted_subnetScore, args)
            
    if num_cores != 1:
          pool.close()
          pool.join()

#print all_permuted_subnetsocre
            # Parse the results into a map of k -> counts
    permuted_sub_score = []
    permuted_gene_degree = []
    permuted_degree_weighted_score = []
    permuted_centality_weighted_score = []

    for (score, degree,weight_degree,weight_cen) in all_permuted_subnetsocre:
           permuted_sub_score = permuted_sub_score + score
           permuted_gene_degree = permuted_gene_degree + degree
           permuted_degree_weighted_score = permuted_degree_weighted_score + weight_degree
           permuted_centality_weighted_score = permuted_centality_weighted_score + weight_cen

# print permuted_centality_weighted_score
    return permuted_sub_score, permuted_gene_degree, permuted_degree_weighted_score, permuted_centality_weighted_score


#################################Find subnetwork significant
def MyStatistics_func(observed_data,observed_data2,data,data2,con,opt):
    ####opt =1 is sum of subnet, opt=2 is median
    if opt == 1:
        q= []
        s =[]
        
        for i in observed_data:
            q.append(sum(i))
            s.append(len(i))
        
        sum_subnetscore = []
        length_subnet = []
        
        for row in data:
            sum_subnetscore.append(sum(row))
            length_subnet.append(len(row))

    elif opt == 2:
         q= []
         s =[]
            
         for i in observed_data:
             q.append(sum(i))
        
         for j in observed_data2:
             s.append(sum(j))
                        
         sum_subnetscore = []
         length_subnet = []
                                
         for row in data:
             sum_subnetscore.append(sum(row))

         for row2 in data2:
             length_subnet.append(sum(row2))

    C = len(q)
    alpha = 1.0 - (1.0 - 0.05) ** (1.0/C)
    
    sum_subnetscore = np.array(sum_subnetscore, dtype = np.float)
    length_subnet = np.array(length_subnet, dtype = np.float)
    #print sum_subnetscore
    ####mean value for emperial data
    #mean_t_q = np.mean(sum_subnetscore)
    #mean_t_s = np.mean(length_subnet)
    std_t_q =  np.std(sum_subnetscore, ddof = 1)
    std_t_s =  np.std(length_subnet, ddof = 1)
    gamma = np.corrcoef(sum_subnetscore,length_subnet)[0, 1]
    h = float(len(sum_subnetscore)) ** (- 1.0 / 6.0)
    ####P vlaue vector
    results = []
    conponent = []
    sig_conponet = []
    count = 0.0
    ####calculate p_value for each input
    for i in range(C):
        w = np.exp(- ( (s[i] - length_subnet) / (np.sqrt(2.0) * h * std_t_s) ) ** 2)
        cd = norm.cdf( ( (q[i] - sum_subnetscore) / (h * std_t_q) - gamma * (s[i] - length_subnet) / (h * std_t_s) ) / np.sqrt(1.0 - gamma * gamma) )
        ####prevent outside the bound
        w[w==0.0] = 1e-323
        cd[cd==0.0] = 1e-323
        denom = sum(w)
        ###calculate p-value
        p_values = 1.0 - (sum( w * cd ) / denom)
        ####store significant component
        if p_values <= alpha:
            sig_conponet.append(con[i])
            count = count + 1.0
        
        significant = p_values <= alpha
        results.append(observed_data[i]+ [significant] + [p_values])
        conponent.append(con[i]+ [p_values] + [significant])

    return alpha,results,conponent,count,sig_conponet


#################Find the degree of gene within each subnetwork function
def Gene_Degree_func(G,ccs):
    """Find the degree of gene within each subnetwork"""
    gene_degree = []
    for subnetwork in ccs:
        H = G.subgraph(subnetwork)
        row = []
        for gene in subnetwork:
            row = row + [H.degree(gene)]
            
        gene_degree.append(row)

    return gene_degree

def Gene_Centrality_func(G,ccs):
    """Find centrality of each gene within subnetwork"""
    gene_centrality = []
    for subnetwork in ccs:
        H = G.subgraph(subnetwork)
        temp = nx.betweenness_centrality(H,normalized=True)
        row = []
        for gene in subnetwork:
            row = row + [temp[gene]]

        gene_centrality.append(row)

    return gene_centrality


def Product_func(s,d):
    conbind = []
    d = Add_one_func(d)
    for i, j in zip(s,d):
        conbind.append(np.array(i)*np.array(j))

    return conbind

def Add_one_func(d):
    New_d = []
    for i in d:
        New_d.append([x+1.0 for x in i])

    return New_d





