"""
This is the original version of the code that is written about in the section `Causal Bayesian NetworkX: Graphs` in the paper "Causal Bayesian NetworkX" by M.D. Pacer, which can be found at https://github.com/michaelpacer/scipy_proceedings/tree/2015/papers/mike_pacer.
"""

def completeDiGraph(nodes):
    """
    Building a max-graph from a set of n nodes.
    This graph has :math:`n^2` edges.
    Variables:
    nodes are a list of strings comprising node names
    """

    G = nx.DiGraph() # Creates new graph
    G.add_nodes_from(nodes) # adds nodes to graph
    edgelist = list(combinations(nodes,2)) 
    # list of directed edges
    edgelist.extend([(y,x) for x,y in edgelist)
    #add symmetric edges
    edgelist.extend([(x,x) for x in nodes]) 
    # add self-loops
    G.add_edges_from(edgelist) # add edges to graph
    return G

def filter_Graph(G,filter_set):
    """
    This allows us to apply a set of filters encoded 
    as closures that take a graph as input
    and return a graph as output.
    """
    graph = G.copy()
    for f in filter_set:
        graph = f(graph)
    return graph

def conditionalSubgraphs(G,condition_list):
    """
    Returns a graph iterator of subgraphs of G 
    meeting conditions in condition_list.

    Variables: 
    G: a graph from which subgraphs will be taken.
    condition_list: a list of condition functions.
    
    Functions in condition_list have i/o defined as
    input: graph, generated as a subgraph of G
    output: Bool, whether graph passes condition
    """

    for edges in powerset(G.edges()):
        G_test = G.copy()
        G_test.remove_edges_from(edges)
        if all([c(G_test) for c in condition_list]):
            
            yield G_test

def new_conditional_graph_set(graph_set,cond_list):
    """
    Returns graph_set & a new iterator which has 
    conditions in cond_list applied to it.
    
    Warning: This function will devour the iterator 
    you include as the `graph_set` input, 
    you need to redeclare the variable as 
    one of the return values of the function.
    
    Thus a correct use would be:    
    a,b = new_conditional_graph_set(a,c)
    
    The following would not be a correct use:
    x,y = new_conditional_graph_set(a,c)
    
    Variables: 
    graph_set: graph-set iterator generator
    cond_list: list conditions
        input: a graph.
        output: boolean value
    """
    
    graph_set_newer, graph_set_test = tee(graph_set,2)
    def gen():
        for G in graph_set_test:
            G_test = G.copy()
            if all([c(G_test) for c in condition_list]):
                yield G_test
    return graph_set_newer, gen()


def create_path_complete_condition(node_pairs):
    """ Cretaes a graph condition that returns true only if a subgraph has a directed path from the first node to the second for every pair of nodes int he node_pairs list.
    """
    def path_complete_condition(G):
        return all([nx.has_path(G,x,y) for x,y in node_pairs])
    return path_complete_condition

def extract_remove_self_loops_filter():
    def remove_self_loops_filter(G):
        graph = G.copy()
        graph.remove_edges_from(graph.selfloop_edges())
        return graph
    return remove_self_loops_filter
