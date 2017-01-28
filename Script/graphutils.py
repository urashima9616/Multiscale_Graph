# -*- coding: utf-8 -*-
'''
Modified Multiscale Network Generator 

The algorithms are based on the work from Alexander Gutfraind and Ilya Safro.

Code to assist the generator

'''

'''
Contains code derived from NetworkX - Copyright (c) by various authors.
'''

import os
import time
import numpy as np
import numpy.random as npr
import random, sys
import networkx as nx
import pdb
import cPickle
import algorithms
import community

np.seterr(all='raise')

timeNow = lambda : time.strftime('%Y_%m_%d__%H_%M_%S', time.localtime())

#those are short methods (like lambda, but supporting pickling)
def a_avg_degree(G):
    return np.average(nx.degree(G).values())
def a_degree_connectivity(G):
    return np.average(nx.average_degree_connectivity(G).values())
def a_s_metric(G):
    return nx.s_metric(G, normalized=False)
def a_eccentricity(G):
    return np.average(nx.eccentricity(G).values())
def a_avg_shortest(G):
    return average_all_pairs_shortest_path_estimate(G, max_num_sources=100)
def a_avg_harmonic(G):
    return average_all_pairs_inverse_shortest_path_estimate(G, max_num_sources=100)
def a_avg_between(G):
    return np.average(nx.betweenness_centrality(G, normalized=True).values())

METRIC_ERROR = -999999
#any negative value is suitable, but not None or -np.inf (might complicate code for statistical analysis)

def average_all_pairs_shortest_path(G):
  if nx.number_connected_components(G)>1:
    return METRIC_ERROR
  else:
    length=nx.all_pairs_shortest_path_length(G)
    sum = 0.0
    count = 0.0
    for key1 in length.keys():
      for key2 in length[key1].keys():
        sum = sum + length[key1][key2]
	count = count + 1
    return sum/count

def average_all_pairs_shortest_path_estimate(G, max_num_sources=100):
    if nx.number_connected_components(G)>1:
        return METRIC_ERROR

    num_sources = min(G.number_of_nodes(), max_num_sources)
    sum = 0.0
    for source in random.sample(G.nodes(), num_sources):
        lengths = nx.single_source_shortest_path_length(G, source)
        sum += np.average(lengths.values())
    sum /= num_sources
    return sum


def average_all_pairs_inverse_shortest_path_estimate(G, max_num_sources=100):
#estimates the ''efficiency'' of the graph: the harmonic mean of the distances
#this is well-defined even in disconnected graphs
    if G.number_of_edges()<1:
        return METRIC_ERROR

    num_sources = min(G.number_of_nodes(), max_num_sources)
    tally = 0.0
    for source in random.sample(G.nodes(), num_sources):
        lengths = nx.single_source_shortest_path_length(G, source)
        tally += sum([1.0/lengths[node] for node in lengths if node!=source])
    tally = num_sources*(1.0/tally)
    return tally

def average_flow_closeness(G):
  if nx.number_connected_components(G)>1:
    return METRIC_ERROR
  else:
    length=nx.algorithms.current_flow_closeness_centrality(G)
    sum = 0.0
    count = 0.0
    for key1 in length.keys():
      sum = sum + length[key1]
      count = count + 1
    return sum/count

def average_eigenvector_centrality(G):
  #warning: this algorithm might not be suitable for disconnected graphs, since it creates additional zero eigenvalues
  if nx.number_connected_components(G)>1:
    return METRIC_ERROR
  else:
    length=nx.algorithms.eigenvector_centrality(G, 500, 0.0001)
    sum = 0.0
    count = 0.0
    for key1 in length.keys():
      sum = sum + length[key1]
      count = count + 1
    return sum/count

def algebraic_distance(G, params={}):
    '''
    takes: graph G, computational parameters

    returns:
    a distance dictionary, d[node1][node2]  giving the distance between the nodes

    ref:
            RELAXATION-BASED COARSENING AND
            MULTISCALE GRAPH ORGANIZATION
            DORIT RON, ILYA SAFRO, AND ACHI BRANDT

    this code is not currently used in the main generator algorithm

    wishlist: use sparse matrices
    '''

    H = nx.convert_node_labels_to_integers(G, discard_old_labels=False)
    metric              = params.get('metric', 'Linfinity')
    num_relaxations_r   = params.get('num_relaxations_r', 10)
    num_test_vectors_K  = params.get('num_test_vectors_K', 20)
    lazy_walk_param_w   = params.get('lazy_walk_param_w', 0.5)

    if metric != 'Linfinity':
        raise Exception('Metric other than Linifinity not implemented')

    distance = {}
    for node1 in H:
        distance[node1] = {}
        for node2 in H:
            if node1 < node2: #save time
               distance[node1][node2] = -np.inf

    LAP      = nx.laplacian(H)
    diag_vec = np.diag(LAP)
    DIAG     = np.diag(diag_vec)
    w_times_Dinv_times_D_minus_LAP = lazy_walk_param_w * np.dot(np.diag([1./el for el in diag_vec]),DIAG-LAP)

    for t in xrange(num_test_vectors_K):
        x = npr.rand(H.number_of_nodes(),1)

        for iteration in xrange(num_relaxations_r):
            x = (1-lazy_walk_param_w)*x + np.dot(w_times_Dinv_times_D_minus_LAP, x)

        for node1 in H:
            for node2 in H:
                dis = abs(x[node1]-x[node2])[0]
                if node1 < node2 and dis > distance[node1][node2]: #to save time, compute just the upper triangle of the matrix
                    distance[node1][node2] = dis

    #generate the distance dictionary in the original node labels, and including the diagonal and lower triangle
    ret = {}
    for u in G:
        ret[u] = {}
        for v in G:
            node1 = H.node_labels[u]
            node2 = H.node_labels[v]
            if node1 < node2:
                ret[u][v] = distance[node1][node2]
            elif node1 > node2:
                ret[u][v] = distance[node2][node1]
            else:
                ret[u][v] = 0.
    return ret

def bfs_distance_with_horizon(G, source, horizon=4, blocked_node=None):
#computes distance from every node to every neighbor at distance at most <horizon> hops
#    nodes further away are considered infinitely away
#no path is allowed through blocked_node
    G_adj = G.adj
    G_neighbors = lambda u: G_adj[u].keys()

    fringe = set(G_neighbors(source))
    distance_source  = {source:0}
    for d in xrange(1, horizon+1):
        new_fringe = []
        for v in fringe:
            if v not in distance_source and v!=blocked_node:
                distance_source[v] = d
                new_fringe += G_neighbors(v)
        fringe = set(new_fringe)

    return distance_source


def color_by_3d_distances(G, verbose):
    import matplotlib.pylab as pylab
    #cm=pylab.get_cmap('Paired')
    #cm=pylab.get_cmap('gist_rainbow')
    cm=pylab.get_cmap('RdBu')  #UFL
    
    if verbose:
        print 'Computing edge colors ...'
    max_dis = 0
    positions = {}
    for u,v,data in G.edges_iter(data=True):
        try:
            u_pos = positions[u]
        except:
            #u_pos = np.array([float(p) for p in G.node[u]['pos'][1:-1].split(',')])
            u_pos = np.array([float(p) for p in G.node[u]['pos'].split(',')])
            positions[u] = u_pos
        try:
            v_pos = positions[v]
        except:
            #v_pos = np.array([float(p) for p in G.node[v]['pos'][1:-1].split(',')])
            v_pos = np.array([float(p) for p in G.node[v]['pos'].split(',')])
            positions[v] = v_pos

        dis = np.sqrt(np.sum(np.power(u_pos-v_pos,2)))
        max_dis = max(max_dis, dis)

        data['dis'] = dis

    for u,v,data in G.edges_iter(data=True):
        dis = data.pop('dis')
        #data['color'] = '"%.3f %.3f %.3f"'%tuple(cm(dis/max_dis)[:3])
        data['color'] = '%.3f %.3f %.3f'%tuple(cm(dis/max_dis)[:3])
        #data['weight'] = 1.0

    return G


def color_new_nodes_and_edges(G, original, params=None):
#add red color to new components.  
    for node in G:
        G.node[node]['label'] = ''
        #d['style'] = 'filled'
        if node in original:
            G.node[node]['color']='black'
        else:
            G.node[node]['color']='blue'
    for u,v,d in G.edges_iter(data=True):
        if original.has_edge(u,v):
            d['color']='black'
        else:
            d['color']='blue'

    return G

def compare_nets(old_G, new_G, metrics=None, params={}):
    '''
    Report on the differences between two networks
    '''
    if metrics == None:
        metrics = default_metrics
    verbose   = params.get('verbose', True)
    precision   = params.get('reporting_precision', 2)
    formatstring = '\t%.'+str(precision)+'f\t%.'+str(precision)+'f\t%.'+str(precision)+'f%%'

    #TODO: at the moment, graph_graph_delta cannot find edges which were deleted then inserted back: it changes the edge attribute data
    errors = {}
    if verbose: 
        delta = graph_graph_delta(old_G, new_G)
        num_changed_nodes = len(delta['new_nodes']) + len(delta['del_nodes'])
        num_changed_edges = len(delta['new_edges']) + len(delta['del_edges'])
        if old_G.number_of_nodes() > 0:
            print 'New or deleted Nodes: %d (%.1f%%)'%(num_changed_nodes, 100*float(num_changed_nodes)/old_G.number_of_nodes())
            print 'New or deleted Edges: %d (%.1f%%)'%(num_changed_edges, 100*float(num_changed_edges)/old_G.number_of_edges())
            print
        print 'Name\t\t\tOld G\tNew G\tRelative Error'
        print 'statistics start ------------------------------------------------------------'
    for met_info in metrics:
        met_name = met_info['name']
        met_func = met_info['function']
        met_wt = met_info['weight']
        if met_info['optional'] > 0 or met_info['practical_node_limit'] < old_G.number_of_nodes():
            continue
        try:
            if verbose: 
                sys.stdout.write(met_name.center(20))
                sys.stdout.flush()
            old_value = met_func(old_G)
            new_value = met_func(new_G)
            if old_value != 0. and abs(old_value-METRIC_ERROR) > 1 and abs(new_value-METRIC_ERROR) > 1:
                error = met_wt*float(new_value-old_value)/old_value
            else:
                error = np.NaN
            if verbose: 
                print formatstring%(old_value,new_value,100*error)
            errors[met_name] = (old_value, new_value, error)
        except Exception,inst:
            print
            print 'Warning: could not compute '+met_name + ': '+str(inst)
    mean_error = np.average([v[2] for v in errors.values() if (v[2]!=np.NaN) and abs(v[2]-METRIC_ERROR) > 1])
    if verbose: 
        print 'statistics end ------------------------------------------------------------'
        print 'Mean difference: %.2f%%'%100*mean_error
        
    return mean_error, errors

def degree_assortativity(G):
#this wrapper helps avoid error due to change in interface name
    if hasattr(nx, 'degree_assortativity_coefficient'):
        return nx.degree_assortativity_coefficient(G)
    elif hasattr(nx, 'degree_assortativity'):
        return nx.degree_assortativity(G)
    else:
        raise ValueError, 'Cannot compute degree assortativity: method not available'

def graph_graph_delta(G, new_G, **kwargs):
#lists the changes in the two graphs, and reports the Jaccard similarity coefficient for nodes and for edges
    new_nodes = []
    del_nodes = []
    new_edges = []
    del_edges = []

    for node in G:
        if node not in new_G:
            del_nodes.append(node)
    for edge in G.edges():
        if not new_G.has_edge(*edge):
            del_edges.append(edge)

    for node in new_G:
        if node not in G:
            new_nodes.append(node)
    for edge in new_G.edges():
        if not G.has_edge(*edge):
            new_edges.append(edge)

    ret = {'new_nodes':new_nodes, 'del_nodes':del_nodes, 'new_edges':new_edges, 'del_edges':del_edges}

    num_nodes_original = G.number_of_nodes()
    num_edges_original = G.number_of_edges()
    if num_nodes_original + len(new_nodes) > 0:
        jaccard_nodes = float(num_nodes_original-len(del_nodes))/(num_nodes_original + len(new_nodes))
    else:
        jaccard_nodes = 0.
    if num_edges_original + len(new_edges) > 0:
        jaccard_edges = float(num_edges_original-len(del_edges))/(num_edges_original + len(new_edges))
    else:
        jaccard_edges = 0.
    ret['jaccard_nodes'] = jaccard_nodes
    ret['jaccard_edges'] = jaccard_edges

    return ret
    
     
def load_graph(path, params={}, list_types_and_exit=False):
    '''reads graph from path, using automatic detection of graph type
       to attempt AUTODETECTION use params['graph_type'] = AUTODETECT
    '''

    loaders = {
            'adjlist':nx.read_adjlist, 
            'adjlist_implicit_prefix':read_adjlist_implicit_prefix, 
            'graph6':nx.read_graph6, 
            'shp':nx.read_shp, 
            'dot':nx.read_dot, 
            'xdot':nx.read_dot, 
            'graph6_list':nx.read_graph6_list,  
            'sparse6':nx.read_sparse6, 
            'edges':nx.read_edgelist, 
            'elist':nx.read_edgelist, 
            'edgelist':nx.read_edgelist, 
            'graphml':nx.read_graphml,        
            'sparse6_list':nx.read_sparse6_list, 
            'gexf':nx.read_gexf,    
            'leda':nx.read_leda,            
            'weighted_edgelist':nx.read_weighted_edgelist,
            'gml':nx.read_gml, 
            'multiline_adjlist':nx.read_multiline_adjlist,
            'yaml':nx.read_yaml,
            'gpickle':nx.read_gpickle,           
            'pajek':nx.read_pajek,}

    raw_loaders = {
            'adjlist':nx.parse_adjlist, 
            'elist':nx.parse_edgelist, 
            'edgelist':nx.parse_edgelist, 
            'gml':nx.parse_gml, 
            'leda':nx.parse_leda,            
            'multiline_adjlist':nx.parse_multiline_adjlist,
            'pajek':nx.parse_pajek,
            }

    known_extensions = {
            'gml':'gml',
            'dot':'dot',
            'xdot':'dot',
            'edges':'edgelist',
            'alist':'adjlist',
            'adjlist':'adjlist',
            'edgelist':'edgelist',
            'elist':'edgelist',
            'pajek':'pajek',}

    if list_types_and_exit:
        return loaders.keys()

    def sane_graph(G, params={}):
        if G.number_of_nodes() > 0:
            return True
        else:
            return False

    G = None
    graph_type = params.get('graph_type', 'AUTODETECT')
    read_params= params.get('read_params', {})
    skip_sanity= params.get('skip_sanity', False)

    if not os.path.exists(path):
        raise ValueError, 'Path does not exist: %s'%path

    if graph_type in loaders:
        if graph_type in ['edges', 'elist', 'edgelist']:
            print "Default weight is 1. To indicate weight, each line should use the format: node1 node2 {'weight':positive_wt}"
        try:
            G = loaders[graph_type](path=path, **read_params, create_using=nx.DiGraph());
            if not sane_graph(G) and not skip_sanity:
                print 'Warning: Sanity test failed!'
                print
                graph_type = None
        except:
            print 'Graph read error.'
            raise

    if G == None and graph_type != 'AUTODETECT':
        raise Exception,'Unable to load graphs of type '+str(graph_type)

    extension_guess = os.path.splitext(path)[1][1:]
    if G == None and extension_guess in known_extensions:
        print 'Attempting auto-detection of graph type.'

        if params.get('verbose', True):
            print 'Warning: Trying to auto-detect graph type by extension'
        graph_type = known_extensions[extension_guess]
        if params.get('verbose', True):
            print 'Guessing type: '+str(graph_type)
        try:
            G = loaders[graph_type](path=path)
            assert sane_graph(G) or skip_sanity
        except Exception, inst:
            print 'Graph read error.  This might be caused by malformed edge data or unicode errors.'
            print inst

    if G == None and graph_type in raw_loaders:
        if params.get('verbose', True):
            print 'Trying raw read...'
        try:
            f = open(path, 'rb')
            lines = f.readlines()
            G = raw_loaders[graph_type](lines=lines)
            del lines
            assert sane_graph(G) or skip_sanity
        except Exception, inst:
            print 'Graph read error:'
            print inst
        finally:
            try:
                f.close()
            except:
                pass

    if G == None:
        if params.get('verbose', True):
            print 'Warning: Trying to guess graph type iteratively: this often FAILS'
        for graph_type in loaders:
            try:
                if params.get('verbose', True):
                    sys.stdout.write(graph_type + '? ')
                G = loaders[graph_type](path=path)
                if sane_graph(G) or skip_sanity:
                    if params.get('verbose', True):
                        print(' Yes!')
                        print 'Successfully detected type: '+str(graph_type)
                    break
                else:
                    if params.get('verbose', True):
                        print(' No.')
                    G = None
                #wishlist: attempt edgelist before adjlist
            except:
                if params.get('verbose', True):
                    print(' No.')

    if G == None:
        raise Exception, 'Could not load graph.  None of the available loaders succeeded.'
    
    #postprocessing
    if graph_type == 'dot':
        G.name = os.path.split(path)[1]  #otherwise the output is terrible

    return G

def powerlaw_mle(G, xmin=6.):
    #estimate the power law exponent based on Clauset et al., http://arxiv.org/abs/0706.1062, 
    #for simplicity, we avoid the MLE calculation of Eq. (3.5) and instead use the approximation of Eq. (3.7)
    #the power law is only applied for nodes of degree > xmin, so it's not suitable for others
    degseq = G.degree().values()

    #print np.array(degseq).transpose()

    if xmin < 6:
        print 'Warning: the estimator uses an approximation which is not suitable for xmin < 6'

    degseqLn = [np.log(xi/(xmin-0.5)) for xi in degseq if xi >= xmin]
    degseqLn.sort() #to reduce underflow. 

    #print degseqLn
    return 1. + len(degseqLn) / sum(degseqLn)


def read_adjlist_implicit_prefix(path, create_using=None):
    '''
    reads network files formatted as:

    "
    15606 45878
    2 3 6 7 
    1 4 6 9 
    "
    and so on.
    first line: num_nodes num_edges
    second lines (and the rest of the lines in the file):
    [implicit node = line number - 1] neighbor1 neighbor2 ... 
    empty lines are degree=0 nodes
    '''
    
    if create_using == None:
        G = nx.Graph()
    else:
        G = create_using()

    try:
        with open(path, 'r') as file_handle:
            header_data = file_handle.next().split(' ')
            node_num = 1
            for line in file_handle:
                line = line.strip()
                if line == '':
                    G.add_node(node_num)
                else:
                    G.add_edges_from([(node_num,int(v)) for v in line.split(' ')])
                node_num += 1
    except Exception,inst:
        if 'node_num' not in locals():
            raise
        raise IOError, 'Parse error on line %d'%(node_num+1)
    
    expected_num_nodes = int(header_data[0])
    expected_num_edges = int(header_data[1])

    if G.number_of_nodes() != expected_num_nodes or G.number_of_edges() != expected_num_edges:
        raise IOError, 'Failed to correctly parse input file. Expected nn=%d,ne=%d; Read %d,%d'%(expected_num_nodes,expected_num_edges,G.number_of_nodes(),G.number_of_edges())

    return G

def safe_pickle(path, data, params=None):
    with open(path, 'wb') as f:
        cPickle.dump(data, f)
        if type(params) != type({}) or params.get('verbose', True): 
            print 'pickled to: '+str(path)


def test_algebraic_distance():
    print 'Testing Algebraic distance'
    #test1: nodes nearby on the path graph should land nearby
    print 'test1 ...'
    G1 = nx.path_graph(10)
    distance1 = algebraic_distance(G1)
    true_distance1 = []
    alg_distance1 = []
    for node1 in G1:
        for node2 in G1:
            if node1 > node2:
                continue
            true_distance1.append(abs(node1-node2))
            alg_distance1.append(distance1[node1][node2])

    val1 = np.corrcoef(true_distance1, alg_distance1)[0,1]
    print 'correlation: %.2f'%val1
    assert val1 > 0.8
    print 'passed.'

    #test2: same for grid graph
    G2=nx.grid_graph(dim=[10,10])
    distance2 = algebraic_distance(G2)
    true_distance2 = []
    alg_distance2 = []
    for node1 in G2:
        for node2 in G2:
            if node1 > node2:
                continue
            true_distance2.append(abs(node1[0]-node2[0]) + abs(node1[1]-node2[1]))
            alg_distance2.append(distance2[node1][node2])

    val2 = np.corrcoef(true_distance2, alg_distance2)[0,1]
    print 'correlation: %.2f'%val2
    assert val2 > 0.5
    print 'passed.'

def test_average_path_length():
    print 'Testing avg path length estimator'
    G = nx.barabasi_albert_graph(300, 5)
    #G = nx.cycle_graph(300)

    estimated_avg = average_all_pairs_shortest_path_estimate(G, max_num_sources=200)

    true_lengths = nx.all_pairs_shortest_path_length(G)
    true_avg = np.average([np.average(true_lengths[node].values()) for node in G])

    print 'Estimate: %f'%estimated_avg
    print 'True:     %f'%true_avg

    assert abs(estimated_avg-true_avg)/true_avg < 0.03
    print 'PASSED'

def test_bfs():
    print 'Testing BFS'
    G = nx.path_graph(5)
    distances_path0 = bfs_distance_with_horizon(G, source=0, horizon=2)
    assert distances_path0[0] == 0
    assert distances_path0[1] == 1
    assert distances_path0[2] == 2
    assert 3 not in distances_path0
    assert 4 not in distances_path0
    distances_path1 = bfs_distance_with_horizon(G, source=1, horizon=2)
    assert distances_path1[0] == 1     
    assert distances_path1[1] == 0     
    assert distances_path1[2] == 1     
    assert distances_path1[3] == 2     
    assert 4 not in distances_path1

    ER100 = nx.erdos_renyi_graph(100, 0.02)
    true_d       = nx.all_pairs_shortest_path_length(ER100)
    cc1 = nx.connected_components(ER100)[0]
    for node1 in cc1:
        horizon_d_node1    = bfs_distance_with_horizon(ER100, source=node1, horizon=4)
        horizon_dinf_node1 = bfs_distance_with_horizon(ER100, source=node1, horizon=1000)
        for node2 in cc1:
            if node2 in horizon_d_node1:
                assert true_d[node1][node2] == horizon_d_node1[node2]
            assert true_d[node1][node2] == horizon_dinf_node1[node2]
    print 'PASSED'

    s='''
    import networkx as nx
    import graphutils
    #G = nx.grid_2d_graph(10, 10)
    G = nx.erdos_renyi_graph(200, 0.2)
    #graphutils.bfs_distance_with_horizon(G, source=(1,5), horizon=10)
    graphutils.bfs_distance_with_horizon(G, source=15, horizon=10)
    '''
    import timeit
    t=timeit.Timer(stmt=s)
    num_trials = 100
    print '%f usec/pass'%(t.timeit(number=num_trials)/num_trials)



def test_inverse_mean_path_length():
    print 'Testing BFS'
    G = nx.erdos_renyi_graph(100, 0.02)
    eff_est = average_all_pairs_inverse_shortest_path_estimate(G, max_num_sources=100)
    print 'Estimate: '+str(eff_est)
    eff_tru = average_all_pairs_inverse_shortest_path_estimate(G, max_num_sources=G.number_of_nodes())
    print 'True:     '+str(eff_tru)
    assert abs(eff_est-eff_tru)/eff_tru < 0.05
    print 'PASSED'

def test_powerlaw_mle():
    print 'Testing Power law MLE estimator'
    G = nx.barabasi_albert_graph(100, 5)
    print 'nn: %d, alpha: %f'%(G.number_of_nodes(),powerlaw_mle(G))
    G = nx.barabasi_albert_graph(1000, 5)
    print 'nn: %d, alpha: %f'%(G.number_of_nodes(),powerlaw_mle(G))
    G = nx.barabasi_albert_graph(10000, 5)
    print 'nn: %d, alpha: %f'%(G.number_of_nodes(),powerlaw_mle(G))
    G = nx.barabasi_albert_graph(100000, 5)
    print 'nn: %d, alpha: %f'%(G.number_of_nodes(),powerlaw_mle(G))
    print 'Expected: 2.9 (or thereabout)'


def write_dot_helper(G, path, encoding='utf-8'):
    #a simplified implementation of dot writer
    #needed in the Windows platform where pygraphviz is not available
    #loses label information
    with open(path, mode='wb') as f:
        header = 'strict graph ' + getattr(G, 'name', 'replica') + ' {\n'.encode(encoding)
        f.write(header) 
        for line in nx.generate_edgelist(G, ' -- ', False):
            line =' %s;\n'%line
            f.write(line.encode(encoding))
        f.write('}\n'.encode(encoding))

def write_graph(G, path, params={}, list_types_and_exit=False):
    '''reads graph from path, using automatic detection of graph type
    '''

    writers = {
            'adjlist':nx.write_adjlist, 
            'dot':nx.write_dot, 
            'xdot':nx.write_dot, 
            'edges':nx.write_edgelist, 
            'elist':nx.write_edgelist, 
            'edgelist':nx.write_edgelist, 
            'weighted_edgelist':nx.write_weighted_edgelist, 
            'graphml':nx.write_graphml,        
            'gml':nx.write_gml, 
            'gpickle':nx.write_gpickle,           
            'pajek':nx.write_pajek,
            'yaml':nx.write_yaml}
    if os.name == 'nt':
        writers['dot'] = write_dot_helper
        writers['xdot'] = write_dot_helper

    if list_types_and_exit:
        return writers.keys()

    write_params = params.get('write_params', {})
    skip_sanity  = params.get('skip_sanity', False)

    graph_type = os.path.splitext(path)[1][1:]

    if graph_type in writers:
        try:
            writers[graph_type](G=G, path=path, **write_params)
        except Exception, inst:
            print 'Graph write error:'
            print inst

            print 'Attempting to write to DOT format'
            nx.write_dot(G, path)
            print 'Done.'
    else:
        raise Exception,'Unable to write graphs of type: '+str(graph_type)


default_metrics = []
default_metrics += [{'name':'num nodes',          'weight':1, 'optional':0, 'practical_node_limit': np.inf, 'function':nx.number_of_nodes}]
default_metrics += [{'name':'density',            'weight':1, 'optional':2, 'practical_node_limit': np.inf, 'function':nx.density}]
default_metrics += [{'name':'num edges',          'weight':1, 'optional':0, 'practical_node_limit': np.inf, 'function':nx.number_of_edges}]
default_metrics += [{'name':'num comps',          'weight':1, 'optional':0, 'practical_node_limit': np.inf, 'function':nx.number_connected_components}]
default_metrics += [{'name':'clustering',         'weight':1, 'optional':0, 'practical_node_limit': np.inf, 'function':nx.average_clustering}]
default_metrics += [{'name':'average degree',     'weight':1, 'optional':0, 'practical_node_limit': np.inf, 'function':a_avg_degree}]
default_metrics += [{'name':'degree assortativity', 'weight':1, 'optional':2, 'practical_node_limit': np.inf, 'function':degree_assortativity}]
default_metrics += [{'name':'degree connectivity', 'weight':1, 'optional':2, 'practical_node_limit': np.inf, 'function':a_degree_connectivity}]

default_metrics += [{'name':'total deg*deg',       'weight':1, 'optional':0, 'practical_node_limit': np.inf, 'function':a_s_metric}]
default_metrics += [{'name':'mean ecc',           'weight':1, 'optional':0, 'practical_node_limit': np.inf, 'function':a_eccentricity}]
#default_metrics += [{'name':'L eigenvalue sum',   'weight':1, 'optional':0, 'practical_node_limit': np.inf, 'function':lambda G: sum(nx.spectrum.laplacian_spectrum(G)).real}]
default_metrics += [{'name':'average shortest path',   'weight':1, 'optional':0, 'practical_node_limit': np.inf, 'function':a_avg_shortest}]
default_metrics += [{'name':'harmonic mean path',   'weight':1, 'optional':0, 'practical_node_limit': np.inf, 'function':a_avg_harmonic}]
#flow_closeness appears to be broken in NX 1.6
default_metrics += [{'name':'avg flow closeness',   'weight':1, 'optional':1, 'practical_node_limit': np.inf, 'function':average_flow_closeness}]
default_metrics += [{'name':'avg eigvec centrality',   'weight':1, 'optional':1, 'practical_node_limit': np.inf, 'function':average_eigenvector_centrality}]
default_metrics += [{'name':'avg between. central.',   'weight':1, 'optional':0, 'practical_node_limit': np.inf, 'function':a_avg_between}]
default_metrics += [{'name':'modularity',          'weight':1, 'optional':0, 'practical_node_limit': np.inf, 'function':community.louvain_modularity}]
default_metrics += [{'name':'powerlaw exp',          'weight':1, 'optional':0, 'practical_node_limit': np.inf, 'function':powerlaw_mle}]
#'optional' runs from 0 (always used) to 5 (never)


if __name__ == '__main__': 
    test_algebraic_distance()
    test_bfs()
    test_average_path_length()
    test_inverse_mean_path_length()
    test_powerlaw_mle()

