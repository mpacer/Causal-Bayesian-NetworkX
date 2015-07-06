"""
This is the original version of the code that is written about in the section `Example: |cbnx| implementation for sprinkler graph` in the paper "Causal Bayesian NetworkX" by M.D. Pacer, which can be found at https://github.com/michaelpacer/scipy_proceedings/tree/2015/papers/mike_pacer.
"""


node_prop_list = [
    ("rain",
     {"state_space":("raining","dry"), "sample_function": "choice",
      "parents":[], 
      "distribution":[.2,.8]
     }),
    ("sprinkler",
     {"state_space":("on","off"),"sample_function": "choice",
      "parents":["rain"], 
      "distribution":{(("rain","raining"),):[.01,.99],
                      (("rain","dry"),):[.4,.6]}
     }),
    ("grass_wet",
     {"state_space":("wet","notWet"),"sample_function": "choice",
      "parents":["rain","sprinkler"],
      "distribution":{(("rain","raining"),("sprinkler","on")):[.99,.01],
                      (("rain","raining"),("sprinkler","off")):[.8,.2],
                      (("rain","dry"),("sprinkler","on")):[.9,.1],
                      (("rain","dry"),("sprinkler","off")):[0,1]}

     })
]

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
        
    return sample_values
   


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
