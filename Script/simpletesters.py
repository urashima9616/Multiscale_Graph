'''
Modified Multiscale Network Generator 

The algorithms are based on the work from Alexander Gutfraind and Ilya Safro.


Test Scripts

'''

import os
import time
import numpy as np
import numpy.random as npr
import random, sys
import networkx as nx
#import matplotlib
#matplotlib.use('PDF')
#import matplotlib.pylab as pylab
#import pylab
import pdb
import cPickle
import graphutils

np.seterr(all='raise')

timeNow = lambda : time.strftime('%Y_%m_%d__%H_%M_%S', time.localtime())
       
def integrity_test():
    import algorithms
    print 'Integrity testing ...'
    graphs = {'karate': nx.generators.karate_club_graph(),
              'er200_025': nx.erdos_renyi_graph(n=200, p=0.25, seed=17),
              'er200_0001': nx.erdos_renyi_graph(n=200, p=0.001, seed=42)}

    params = {'verbose':True,
              'node_edit_rate': [0],
              'edge_edit_rate': [0],
              'node_growth_rate': [0],
              'verbose':False}
    for name,G in graphs.items():
        print name
        replica = algorithms.generate_graph(original=G, params=params)

        diff    = graphutils.graph_graph_delta(G, replica)
        assert diff['new_nodes'] == []
        assert diff['del_edges'] == []
        assert diff['new_nodes'] == []
        assert diff['del_edges'] == []
        
    print 'Integrity test: PASSED'

def sanity_test(G, params=None):
    ok = True
    if G.selfloop_edges() != []:
        print 'Warning: self-loops detected'
        ok = False
    if G.number_of_edges() == 0:
        print 'Warning: no edges'
        ok = False
    if G.is_directed():
        print 'Warning: the algorithm DOES NOT support directed graphs for now'
        ok = False
    if G.has_node(None):
        print 'Node with label None is not allowed in the graph'
        ok = False

    return ok

def smoke_test():
    import algorithms
    print 'Smoke testing ...'
    graphs = {'karate': nx.generators.karate_club_graph(),
              'er200_025': nx.erdos_renyi_graph(n=200, p=0.25, seed=42),
              'er200_0001': nx.erdos_renyi_graph(n=200, p=0.001, seed=42)}

    params = {'verbose':False,
              'node_edit_rate': [0.1/(1.+i) for i in xrange(100)],
              'edge_edit_rate': [0.1/(1.+i) for i in xrange(100)],
              'node_growth_rate': [0.1/(1.+i) for i in xrange(100)]}
    for name,G in graphs.items():
        print name
        #print '  nn=%d,ne=%d'%(G.number_of_nodes(), G.number_of_edges())
        replica = algorithms.generate_graph(original=G, params=params)
        #print '  nn=%d,ne=%d'%(replica.number_of_nodes(), replica.number_of_edges())
        assert G.selfloop_edges() == []

    print 'Smoke test: PASSED'
    print


def validate_params(params):
    valid_params = {'accept_chance_edges':None,
                    'algorithm':None,
                    'coarsening_density_limit':None, #TODO: document
                    'component_is_edited':None, #TODO: document
                    'deep_copying':None,
                    'deferential_detachment_factor':None,
                    'do_coarsen_tester':None, #TODO: document
                    'do_uncoarsen_tester':None, #TODO: document
                    'dont_cutoff_leafs':None,
                    'edge_edit_rate':None,
                    'edge_growth_rate':None,
                    'edge_welfare_fraction':None,
                    'edit_edges_tester':None, #TODO: document
                    'edit_nodes_tester':None, #TODO: document
                    'enforce_connected':None,
                    'fine_clustering':None,
                    'locality_algorithm':None,
                    'locality_bias_correction':None,
                    'long_bridging':None, #TODO: document
                    'maintain_edge_attributes':None,
                    'maintain_node_attributes':None,
                    'memoriless_interpolation':None,    
                    'minorizing_node_deletion':None, #TODO: document
                    'new_edge_horizon':None,
                    'node_edit_rate':None,
                    'node_growth_rate':None,
                    'num_deletion_trials':None, #TODO: document
                    'num_insertion_trials':None, #TODO: document
                    'num_pairs_to_sample':None, #TODO: document
                    'num_snapshots':None,
                    'num_v_cycles':None,
                    'post_processor':None,
                    'preserve_clustering_on_deletion':None,
                    'revise_graph_tester':None, #TODO: document
                    'retain_intermediates':None,
                    'seed_threshold_1':None, #TODO: document
                    'seed_threshold_2':None, #TODO: document
                    'stats_report_on_all_levels':None, #presumes retain_intermediates
                    'triangle_distribution_limit':None, #TODO: document
                    'suppress_warnings':None,
                    'verbose':None,
                    }
    bad_params = False
    for k in params:
        if k not in valid_params and (not str(k).startswith('_')):
            print 'Unknown parameter: %s'%k
            print
            bad_params = True
    
    if params.get('memoriless_interpolation', False) and params.get('deep_copying', True):
        print 'Memoriless_interpolation=True requires deep_copying=False'
        bad_params = True

    if params.get('stats_report_on_all_levels', False) and not params.get('retain_intermediates'):
        print 'Making retain_intermediates=True'
        params['retain_intermediates'] = True


    if bad_params:
        print 'Existing ...'
        sys.exit(1)


if __name__ == '__main__': 
    smoke_test()
    integrity_test()
