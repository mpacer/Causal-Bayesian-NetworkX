import networkx as nx
import json
from networkx.readwrite import json_graph
from itertools import chain, combinations
# from earthquake_loglikelihood import ll_per_graph

def powerset(iterable):
#    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
#
    len_powerset = 0
    powerset_vals = chain.from_iterable(combinations(s, r) for r in range(len(s)+1))
    return powerset_vals


def clean_json_adj_load(file_name):
    with open(file_name) as d:
        json_data = json.load(d)
    H = json_graph.adjacency_graph(json_data)
    for edge_here in H.edges():
        del(H[edge_here[0]][edge_here[1]]["id"])
    return H

def clean_json_adj_loads(json_str):
    json_data = json.loads(json_str)
    H = json_graph.adjacency_graph(json_data)
    for edge_here in H.edges():
        del(H[edge_here[0]][edge_here[1]]["id"])
    return H

def intervention_effects(graph):    
    f = lambda x: x[0].endswith("int")
    return  [x for x in graph.edges() if f(x)]            

def cause_observation_pairings(graph):    
    f = lambda x: x[0].endswith("★") and x[1].endswith("out")
    return  [x for x in graph.edges() if f(x)]

def hidden_cause_pairs(graph):
    f = lambda x: x[0].endswith("★") and x[1].endswith("★")
    return [x for x in graph.edges() if f(x)]
    

            
def completeDiGraph(nodes):
    """
    returns a directed graph with all possible edges
    
    Variables:
    nodes are a list of strings that specify the node names
    """
    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    edgelist = list(combinations(nodes,2))
    edgelist.extend([(y,x) for x,y in list(combinations(nodes,2))])
    edgelist.extend([(x,x) for x in nodes])
    G.add_edges_from(edgelist)
    return G

def filter_maxGraph(G,filter_set):
    graph = G.copy()
    for f in filter_set:
        graph = f(graph)
    return graph

def partialConditionalSubgraphs(G,edge_set,condition_list):
    try: 
        condition_list[0]
    except TypeError:
        raise TypeError("""
        Subsampling from a graph requires passing in a list of conditions encoded
        as first-class functions that accept networkX graphs as an input and return boolean values.""")
    edge_powerset = powerset(edge_set)
 
    for edges in powerset(edge_set):
        G_test = G.copy()
        G_test.remove_edges_from(edges)
        if all([c(G_test) for c in condition_list]):
            yield G_test

def conditionalSubgraphs(G,condition_list):
    try: 
        condition_list[0]
    except TypeError:
        raise TypeError("""
        Subsampling from a graph requires passing in a list of conditions encoded
        as first-class functions that accept networkX graphs as an input and return boolean values.""")
    edge_powerset = powerset(G.edges())
 
    for edges in powerset(G.edges()):
        G_test = G.copy()
        G_test.remove_edges_from(edges)
        if all([c(G_test) for c in condition_list]):
            
            yield G_test

def create_path_complete_condition(transmit_node_pairs):
    def path_complete_condition(G):
        return all([nx.has_path(G,x,y) for x,y in transmit_node_pairs])
    return path_complete_condition

def create_no_input_node_condition(node_list):
    def no_input_node_condition(G):
        return all([G.in_degree(y)==0 for y in node_list])
    return no_input_node_condition

def create_no_self_loops():
    def no_self_loops(G):
        return not(any([(y,y) in G.edges() for y in G.nodes()]))
    return no_self_loops
    
def create_explicit_parent_condition(parentage_tuple_list):
    """
    This states for a child node, what its explicit parents are.
    """
    def explicit_parent_condition(G):
        return all(
            [sorted(G.in_edges(y[0])) == sorted([(x,y[0]) for x in y[1]]) 
             for y in parentage_tuple_list])
    return explicit_parent_condition

def create_explicit_child_condition(parentage_tuple_list):
    """ 
    This states for a parent node, what its explicit children are.
    """
    def explicit_child_condition(G):
        return all(
            [sorted(G.out_edges(y[0])) == sorted([(y[0],x) for x in y[1]]) 
             for y in parentage_tuple_list])
    return explicit_child_condition

def create_no_direct_arrows_condition(node_pair_list):
    def no_direct_arrows_condition(G):
        return not(any([y in G.edges() for y in node_pair_list]))
    return no_direct_arrows_condition

def create_no_output_node_condition(node_list):
    def no_output_node_condition(G):
        return all([G.out_degree(y)==0 for y in node_list])
    return no_output_node_condition


def extract_remove_self_loops():
    def remove_self_loops(G):
        graph = G.copy()
        graph.remove_edges_from(graph.selfloop_edges())
        return graph
    return remove_self_loops

def extract_remove_inward_edges(exceptions_from_removal):
    """
    This covers both no input and explicit_child_parentage.
    """
    
    def remove_inward_edges(G):
        graph = G.copy()
        list_of_children = [x[0] for x in exceptions_from_removal if len(x[1]) > 0]
        list_of_orphans = [x[0] for x in exceptions_from_removal if len(x[1]) == 0]
        
        for orphan in list_of_orphans:
            graph.remove_edges_from([edge for edge in graph.edges() if edge[1] == orphan])
        
        for child in list_of_children:
            current_edges = graph.in_edges(child)
            valid_edges = [(y,x[0]) for x in exceptions_from_removal if x[0] == child  for y in x[1]]
            graph.remove_edges_from([edge for edge in current_edges if edge not in valid_edges])
        
        return graph
    return remove_inward_edges

def extract_remove_outward_edges(exceptions_from_removal):
    """
    This covers both no_output and explicit_parent_offspring.
    """
    
    def remove_outward_edges(G):
        graph = G.copy()
        list_of_parents = [x[0] for x in exceptions_from_removal if len(x[1]) > 0]
        list_of_barrens = [x[0] for x in exceptions_from_removal if len(x[1]) == 0]

        for barren in list_of_barrens:
            graph.remove_edges_from([edge for edge in graph.edges() if edge[0] == barren])
            
        for parent in list_of_parents:
            current_edges = graph.out_edges(parent)
            valid_edges = [(x[0],y) for x in exceptions_from_removal if x[0] == parent for y in x[1]]
            graph.remove_edges_from([edge for edge in current_edges if edge not in valid_edges])
            
        return graph
    return remove_outward_edges

# def add_edge_attribute(graph,edge,attribute_name,attribute_value):
#     graph[edge[0]][edge[1]][attribute_name]=attribute_value
#     pass
    
# def add_multiple_edge_attributes(graph,edge_list,attribute_name,attribute_value):
#     for edge in edge_list:
#         add_edge_attribute(graph,edge,attribute_name,attribute_value)
#     pass

# def add_gamma_attribute_values(graph,edge_list,base_rate,scale):
#     pass
