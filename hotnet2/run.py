#!/usr/bin/env python

# Load required modules
import json
import sys, os, shutil
import numpy as np

# Load local modules
import heat as hnheat, hotnet2 as hn, hnio, stats, permutations as p, My_Function as my, networkx as nx
from delta import get_deltas_for_network, get_deltas_for_heat
from constants import *

def run_helper(args, infmat, full_index2gene, G, nname, pnp, heat, hname, addtl_genes, get_deltas_fn, infmat_name="PPR", max_cc_sizes=[5, 10, 15, 20], verbose=0):
    """Helper shared by runHotNet2 and runClassicHotNet."""
    
    # Perform delta selection (if necessary)
    if args.deltas:
        deltas = args.deltas
    else:
        deltas = get_deltas_fn(full_index2gene, heat, args.network_permutations, args.num_cores, infmat, addtl_genes, pnp, infmat_name, max_cc_sizes, verbose)

    sim, index2gene = hn.similarity_matrix(infmat, full_index2gene, heat, True, verbose=verbose)

    results = []
    """Change from here """
    #permuted_sub_score_delta = []
    #subnet_geneScore_delta = []
    #Alpha_sig_delta = []
    #Subnetscore_sig_delta = []
    #Conponet_sig_delta = []
    #sig_Conp_delta = []
    #sig_Count_delta = []
    #degree_delta = []
    #permuted_degree_delta = []
    #degree_weighted_score_delta = []
    #cent_weighted_score_delta  = []
    My_results = []
    """end here"""
    
    for delta in deltas:
        
        """Change from here"""
        simG = hn.weighted_graph(sim, index2gene, delta, directed=True)
        ccs = hn.connected_components(simG, args.min_cc_size)

        """end here"""
        # calculate significance (using all genes with heat scores)
        if verbose > 4:
            print "* Performing permuted heat statistical significance..."
            print "\t- Using no. of components >= k (k \\in",
            print "[%s, %s]) as statistic" % (min(HN2_STATS_SIZES), max(HN2_STATS_SIZES))

        heat_permutations = p.permute_heat(heat, full_index2gene.values(),
                                           args.heat_permutations, addtl_genes,
                                           args.num_cores)

#print heat_permutations
                
#print heat_permutations[0]

        """change from here"""
        permuted_sub_score, permuted_gene_degree, permuted_degree_weighted_score, permuted_centality_weighted_score = my.calculate_permuted_subnetScore(infmat, full_index2gene, heat_permutations, delta, verbose, args.min_cc_size, G, args.num_cores)
        
        #permuted_sub_score_delta.append((permuted_sub_score, delta))
        #permuted_degree_delta.append((permuted_gene_degree,delta))
        
        #print permuted_sub_score_delta
        ####print permuted_degree_weighted_score
        #print permuted_degree_delta
        ####calcaule the degree of gene in each permuted subnetwork
        """ End here """



        sizes2counts = stats.calculate_permuted_cc_counts(infmat, full_index2gene,
                                                          heat_permutations, delta, HN2_STATS_SIZES, True,
                                                          args.num_cores)
                
        real_counts = stats.num_components_min_size(simG, HN2_STATS_SIZES)
        size2real_counts = dict(zip(HN2_STATS_SIZES, real_counts))
        sizes2stats = stats.compute_statistics(size2real_counts, sizes2counts,
                                               args.heat_permutations)
                                               #print ccs
        # sort ccs list such that genes within components are sorted alphanumerically, and components
        # are sorted first by length, then alphanumerically by name of the first gene in the component
        ccs = [sorted(cc) for cc in ccs]
        ccs.sort(key=lambda comp: comp[0])
        ccs.sort(key=len, reverse=True)
        
        #print ccs
        """changing here"""
        #####calculate all the gene score for gene list
        genelist, scoreVec = my.score_vec(infmat,full_index2gene, heat)
        #gene_socre = np.column_stack((genelist, scoreVec))
        subnet_geneScore = my.Matchscore_Permuted_fun(genelist,scoreVec,ccs)
        #print subnet_geneScore
        ####output the subnetwork with score only
        #print ccs
        #subnet_geneScore_delta.append((subnet_geneScore,delta))
        
        ####calculate the gene degree within each subnetwork
        gene_degree = my.Gene_Degree_func(G,ccs)
        #degree_delta.append((gene_degree,delta))
        
        #####find score weighted by degree
        degree_weighted_score = my.Product_func(subnet_geneScore,gene_degree)
        #degree_weighted_score_delta.append((degree_weighted_score,delta))
        
        #####find score of gene weighted by centrality
        gene_centrality = my.Gene_Centrality_func(G,ccs)
        cent_weighted_score = my.Product_func(subnet_geneScore,gene_centrality )
        
        #cent_weighted_score_delta.append((cent_weighted_score,delta))
        
        #####calculate significant test
        """Sum and Size of subnetwork"""
        Alpha_sig, Subnetscore_sig, Conponet_sig,count,sig_conponet = my.MyStatistics_func(subnet_geneScore,gene_degree, permuted_sub_score, permuted_gene_degree,ccs,1)
        
        """Degree weighted sum and size of subnetwork"""
        Alpha_sig2, Subnetscore_sig_dwss, Conponet_sig_dwss,count_dwss,sig_conponet_dwss = my.MyStatistics_func(degree_weighted_score, gene_degree, permuted_degree_weighted_score, permuted_gene_degree,ccs,1)
        
        """Centrality weighted sum and size of subnewtwork"""
        Alpha_sig3, Subnetscore_sig_cwss, Conponet_sig_cwss,count_cwss,sig_conponet_cwss = my.MyStatistics_func(cent_weighted_score,gene_degree, permuted_centality_weighted_score, permuted_gene_degree,ccs,1)
        
        
        
        """Sum and Sum of degree of gene within subnetwork"""
        Alpha_sig1, Subnetscore_sig_sd, Conponet_sig_sd,count_sd,sig_conponet_sd = my.MyStatistics_func(subnet_geneScore,gene_degree, permuted_sub_score, permuted_gene_degree,ccs,2)
        
        """Degree weighted sum and Sum of degree of gene within subnetwork"""
        Alpha_sig4, Subnetscore_sig_dwsd, Conponet_sig_dwsd,count_dwsd,sig_conponet_dwsd = my.MyStatistics_func(degree_weighted_score,gene_degree, permuted_degree_weighted_score, permuted_gene_degree,ccs,2)
        
        """Centrality weighted sum and Sum of degree of gene within subnetwork"""
        Alpha_sig5, Subnetscore_sig_cwsd, Conponet_sig_cwsd,count_cwsd,sig_conponet_cwsd = my.MyStatistics_func(cent_weighted_score,gene_degree, permuted_centality_weighted_score, permuted_gene_degree,ccs,2)
        
        
        #####store all the result
        ####this is My results
        My_results.append( (permuted_sub_score, subnet_geneScore, Alpha_sig, Subnetscore_sig, Conponet_sig,count, sig_conponet, gene_degree, permuted_gene_degree, cent_weighted_score, degree_weighted_score, permuted_centality_weighted_score, permuted_degree_weighted_score, Subnetscore_sig_sd, Conponet_sig_sd,count_sd,sig_conponet_sd,Subnetscore_sig_dwss, Conponet_sig_dwss,count_dwss,sig_conponet_dwss, Subnetscore_sig_cwss, Conponet_sig_cwss, count_cwss, sig_conponet_cwss, Subnetscore_sig_dwsd, Conponet_sig_dwsd,count_dwsd,sig_conponet_dwsd, Subnetscore_sig_cwsd, Conponet_sig_cwsd,count_cwsd,sig_conponet_cwsd, delta) )
        
        #print Subnetscore_sig
        """end here """
      
        results.append( (ccs, sizes2stats, delta) )
        ####this is My results
        
    return results,My_results

def get_deltas_hotnet2(full_index2gene, heat, num_perms, num_cores, _infmat, _addtl_genes,
                       permuted_networks_path, infmat_name, max_cc_sizes, verbose):
    # find smallest delta
    deltas = get_deltas_for_network(permuted_networks_path, heat, infmat_name, full_index2gene,
                                       MAX_CC_SIZE, max_cc_sizes, False, num_perms, num_cores, verbose)

    # and run HotNet with the median delta for each size
    return [np.median(deltas[size]) for size in deltas]

def get_deltas_classic(full_index2gene, heat, num_perms, num_cores, infmat, addtl_genes, min_cc_size, max_cc_size, verbose):
    # find delta that maximizes # CCs of size >= min_cc_size for each permuted data set
    deltas = get_deltas_for_heat(infmat, full_index2gene, heat, addtl_genes, num_perms, NUM_CCS,
                                    [min_cc_size], True, num_cores, verbose)

    # find the multiple of the median delta s.t. the size of the largest CC in the real data
    # is <= MAX_CC_SIZE
    medianDelta = np.median(deltas[min_cc_size])

    sim, index2gene = hn.similarity_matrix(infmat, full_index2gene, heat, False)

    for i in range(1, 11):
        G = hn.weighted_graph(sim, index2gene, i*medianDelta)
        largest_cc_size = max([len(cc) for cc in hn.connected_components(G)])
        if largest_cc_size <= max_cc_size:
            break

    # and run HotNet with that multiple and the next 4 multiples
    return [i*medianDelta for i in range(i, i+5)]
