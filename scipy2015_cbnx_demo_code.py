import numpy as np
import networkx as nx
from itertools import chain, combinations, tee
from graph_enumerator import powerset


def completeDiGraph(nodes):
    """
    returns a directed graph with all possible edges
    
    Variables:
    nodes are a list of strings that specify the node names
    """
    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    edgelist = list(combinations(nodes,2))
    edgelist.extend([(y,x) for x,y in edgelist])
    edgelist.extend([(x,x) for x in nodes])
    G.add_edges_from(edgelist)
    return G

def conditionalSubgraphs(G,condition_list):
    """
    Returns a graph generator/iterator such that any conditions specified in condition_list 
    are met by some subgraph of G.
    This is intended to be used in conjunction with completeDiGraph or any graph which subgraphs 
    are expected to be taken.
    
    Variables: 
    G is a graph from which subgraphs will be taken.
    condition_list is a list of first order functions that define conditions subgraphs of G need to meet to be output.
    Functions in condition_list should return a single boolean value for every graph passed into them.
    """
    try: 
        condition_list[0]
    except TypeError:
        raise TypeError("""
        Subsampling from a graph requires passing in a list of conditions encoded
        as first-class functions that accept networkX graphs as an input and return boolean values.""")
    
    for edges in powerset(G.edges()):
        G_test = G.copy()
        G_test.remove_edges_from(edges)
        if all([c(G_test) for c in condition_list]):
            yield G_test

def create_no_self_loops_condition():
    """
    This factory allows us to specify that there are no valid self-loops
    This returns a function that takes an graph argument (G). 
    
    NB: This is a common assumption of causal graphs, because they are not considered to be extended through time.
    """

    def no_self_loops_condition(G):
        return not(any([(y,y) in G.edges() for y in G.nodes()]))
    return no_self_loops_condition
            
def create_path_complete_condition(transmit_node_pairs):
    """
    This factory allows us to specify that there are valid directed paths between pairs of nodes.
    This returns a function that takes an graph argument (G) 
    and verifies that for the list of node pairs the graph meets those dependency conditions. 
    
    NB: This is useful for making known indirect dependencies explicit.
    
    Variables:
    node_list is a list of 2-tuples of nodes that will have valid direct paths 
    from the first of the nodes to the second.
    """

    def path_complete_condition(G):
        return all([nx.has_path(G,x,y) for x,y in transmit_node_pairs])
    return path_complete_condition

def create_no_input_node_condition(node_list):
    """
    This factory allows us to specify that no directed can be directed into a set of nodes.
    This returns a function that takes an graph argument (G) and verifies that 
    none of the nodes in node_list are child nodes. 
    
    NB: This is useful for making interventions explicit over a set of graphs.
    
    Variables:
    node_list is a list of nodes that will have no parents
    """
    
    def no_input_node_condition(G):
        return all([G.in_degree(y)==0 for y in node_list])
    return no_input_node_condition

def new_conditional_graph_set(graph_set,condition_list):
    """
    This returns a copy of the old graph_set and a new graph generator which has 
    the conditions in condition_list applied to it.
    
    Warning: This function will devour the iterator that you include as the graph_set input, 
    you need to redeclare the variable as one of the return values of the function.
    
    Thus a correct use would be:    
    a,b = new_conditional_graph_set(a,c)
    
    The following would not be a correct use:
    x,y = new_conditional_graph_set(a,c)
    
    Variables: 
    graph_set is a graph-set generator
    condition_list is a list of first order functions returning boolean values when passed a graph.
    """
    
    try: 
        condition_list[0]
    except TypeError:
        raise TypeError("""
        Subsampling from a graph requires passing in a list of conditions encoded
        as first-class functions that accept networkX graphs as an input and return boolean values.""")
    graph_set_newer, graph_set_test = tee(graph_set,2)
    def gen():
        for G in graph_set_test:
            G_test = G.copy()
            if all([c(G_test) for c in condition_list]):
                yield G_test
    return graph_set_newer, gen()

def extract_remove_self_loops():
    def remove_self_loops(G):
        graph = G.copy()
        graph.remove_edges_from(graph.selfloop_edges())
        return graph
    return remove_self_loops

def filter_Graph(G,filter_set):
    graph = G.copy()
    for f in filter_set:
        graph = f(graph)
    return graph


def sample_from_graph(G,func_dictionary=None,k = 1):
    """
    This is the function that samples from the rich networkX Bayes Net graph using the parameterization specified in the node attributes.
    
    Variables:
    G is the graph being sampled from.
    k is the number of samples.
    """
    if func_dictionary == None:
        func_dictionary = {"choice": np.random.choice}

    nodes_dict = G.nodes(data = True)
    node_ids = np.array(G.nodes())
    state_spaces = [(node[0],node[1]["state_space"]) for node in nodes_dict]
    orphans = [node for node in nodes_dict if node[1]["parents"]==[]]
    sample_values = np.empty([len(state_spaces),k],dtype='U20')
    sampled_nodes = []

    for node in orphans:
        ## sample k values for all orphan nodes
        samp_func = string_to_sample_function(node[1]["sample_function"],func_dictionary)
        samp_states = node[1]["state_space"]
        samp_distribution = node[1]["distribution"]
        samp_index = G.nodes().index(node[0])
        sample_values[samp_index,:]  = samp_func(samp_states,size=[1,k],p=samp_distribution)
        sampled_nodes.append(node[0])
        
    while set(sampled_nodes) < set(G.nodes()):
        nodes_to_sample = check_if_parents_filled(G,sampled_nodes)
        #nodes_to_sample returns a list of node names that need to be sampled
        
        for n in nodes_to_sample:
            #extracts the indices of the parents of the node to be sampled and their values
            parent_indices = [(parent,G.nodes().index(parent)) for parent in G.node[n]["parents"]]
            parent_vals = [(parent[0],sample_values[parent[1],:]) for parent in parent_indices]
            
            #extracts sample index
            samp_index = G.nodes().index(n)
            sample_values[samp_index,:] = conditional_sampling(G,n,parent_vals,func_dictionary,k)
            sampled_nodes.append(n)
        
    return {node:sample_values[G.nodes().index(node)] for node in sampled_nodes}
       
    
    
def check_if_parents_filled(G,sampled_nodes):
    """
    This function will return those nodes who have not yet been sampled, whose parents have been sampled.
    Variables:
    G is a networkX graph
    sampled_nodes are a list of node names
    """
    check_nodes = [x for x in G.nodes() if x not in sampled_nodes]
    nodes_to_be_sampled = []
    for node in G.nodes(data = True):
        if (node[0] in check_nodes) & (node[1]["parents"]<=sampled_nodes):
            nodes_to_be_sampled.append(node[0])
        
    if len(nodes_to_be_sampled)==0: 
        raise RuntimeError("You should never be running this when no values are returned")
    return nodes_to_be_sampled

def nodeset_query(G,node_set,attrib=[]):
    """
    This is a helper function for querying particular attributes from a node  
    Variables:
    G is a networkX style graph
    node_set is a list of node names that are in G
    attrib are a list of attributes associated with the nodes in G
    """
    if len(attrib)==0:
        return [node for node in G.nodes(data = True) if node[0] in node_set]
    else:
        return_val = []
        for node in G.nodes(data=True):
            if node[0] in node_set:
                return_val.append((node[0],{attr:node[1][attr] for attr in attrib}))
        return return_val
    
    
def conditional_sampling(G,node,parent_vals,func_dictionary, k = 1):
    """
    This function takes a graph as input, a node to sample from in that graph and a set of values for the parents of that node.
    This function should not be consulted for variables without any parents.
    Variables: 
    G is a networkX style graph
    node is a node in G
    parent_vals are the values of the parents of node realized k times
    returns an array of values 
    """
    
    try: node in G
    except KeyError:
        print("{} is not in graph".format(n))
    
    output = np.empty(k,dtype="U20")
    for i in np.arange(k):
        par_val_list = []
        for parent in parent_vals:
            par_val_list.append(tuple([parent[0],parent[1][i]]))
        samp_distribution = G.node[node]["distribution"][tuple(par_val_list)]

    
        samp_func = string_to_sample_function(G.node[node]["sample_function"],func_dictionary)
        samp_states = G.node[node]["state_space"]
#         output.append(samp_func(samp_states,size=[1],p=samp_distribution))
        temp_output = samp_func(samp_states,size=1,p=samp_distribution)
        output[i] = temp_output[0]
    return output

def string_to_sample_function(func_name, func_dictionary=None):
    """
    This allows the function to be passed in as a string that is mapped to a first-class function to other methods.
    sample_function is a string that maps onto a function in the dictionary defined below.
    This takes two arguments a func_dictionary 
    """
    if func_dictionary == None:
        func_dictionary = {"choice": np.random.choice}
        
    try: func_dictionary[func_name]
    except KeyError:
        print("{} is not a function defined in the dictionary you passed.".format(func_name))
    
    return func_dictionary[func_name]

def print_prob_est(test):
    for key,value in test.items():
        for unique_element in set(value):
            prob_est = sum(sum([value==unique_element]))/len(value)
            print("p̂({}={}) = {} ± {:.2e}".format(key,unique_element,prob_est,np.sqrt(prob_est/len(value))))
        print("\n")
