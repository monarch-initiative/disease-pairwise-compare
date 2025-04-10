# Imports
import os
import argparse
import scipy
import pickle
import gzip
import copy
import multiprocessing as mp
import numpy as np
import pandas as pd
import networkx as nx
from collections import Counter
from typing import List, Optional
from pydantic import BaseModel

# Data vis
import matplotlib.pyplot as plt
import xarray

# Pheval phenopacket utils, Ontology tools
from pheval.utils.phenopacket_utils import phenopacket_reader
from pheval.utils.phenopacket_utils import PhenopacketUtil
from semsimian import Semsimian


# ?
class Patient(BaseModel):
    sample_name: str
    phenopacket_path: str
    
    # Pulled from phenopacket
    phenotype_ids: list
    phenotype_count: int
    
    disease_name: str
    disease_id: str
    
    gene_symbol: str
    gene_id: str
    
    # For sssom mappings
    disease_mapped: Optional[str] = None
    gene_mapped: Optional[str] = None


def read_sssom_to_lookup(fpath, exact_match=True):
    """
    Assumes monarch initiative header style of # and then first line is header
    """
    
    sssom_map = {}
    with open(fpath, 'r') as infile:
        header = 0
        cols = {}
        for line in infile:

            
            line = line.strip('\r').strip('\n')
            
            # Weird thing happening with the monarch sssom files (new line character most likely as the final line)
            if len(line) == 0:
                break
                
            if line[0] == "#":
                continue
            
            header += 1
            cols = line.split('\t')
            if header == 1:
                col_inds = {v:i for i,v in enumerate(cols)}
                continue
            
            # Our actual data here
            map_key = cols[col_inds["object_id"]]
            map_val = cols[col_inds["subject_id"]]
            map_type = cols[col_inds["predicate_id"]]
            
            # Only deal with exact matches
            if map_type != "skos:exactMatch":
                continue
            
            sssom_map.update({map_key:map_val})
    
    val_count = len(set(list(sssom_map.values())))
    print("- {} unique terms mapping to {} unique terms...".format(format(len(sssom_map), ','), format(val_count, ',')))
    return sssom_map


def divide_workload(data_list, num_proc: int=1) -> list:
    """
    Meant to divide up the elements in data_list into num_proc equal portions
    by iteratively adding each element to a basket never repeating the same basket until all baskets have an equal amount
    If num_proc == 1 then the original input list will be returned nested in a top layer list i.e. [data_list]
    """

    # Deal with our edge case at the very begginning which then is used as input into the second potential edge case
    ndata_elements = len(data_list)
    if ndata_elements < num_proc:
        num_proc = ndata_elements

    # Edge case
    if num_proc <= 1:
        return [data_list]
    else:
        baskets = [[] for i in range(0, num_proc)]
        index_count = 0
        for d in data_list:
            baskets[index_count].append(d)
            if index_count == (num_proc-1):
                index_count = 0
            else:
                index_count += 1

        #print("- Workload divided into {} portions with each portion recieving {} elements respectively...".format(num_proc, [format(len(b), ',') for b in baskets]))
        return baskets



####################################
### Monarch kg data extract code ###
def get_entities_from_file(infile_path: str, select_data: set, match_column: str ='category') -> list[dict]:
    """
    Expects a tabular file with the header as the first line. Will pull rows from the file where
    where the column whos name matches match_column are found within the select_data 
    (ie we only want rows where column-xyz matches a name within the select_data argument)
    
    Allows one to pull specific node and edge types from the monarch knowledge graph
    """
    
    # Instead of using pandas to load into memory, we will be memory conscious and load in only select nodes
    return_data = []
    with open(infile_path) as infile:
        header_info = {v:i for i,v in enumerate(infile.readline().strip('\r').strip('\n').split('\t'))}
        
        # Ensure relevant columns exist in header before we proceed
        if (match_column not in header_info) or ("id" not in header_info):
            print("- ERROR Header line is expected to be first line of file with columns 'category' and 'id'")
            return None
        
        # Loop through file and pull nodes
        for line in infile:
            cols = line.strip('\r').strip('\n').split('\t')
            entity = cols[header_info[match_column]]
            if entity not in select_data:
                continue
            
            # Convert row into dictionary and append to list
            return_data.append({k:cols[v] for k,v in header_info.items()})
            
    print("- Entities returned {} of types...".format(format(len(return_data), ',')))
    for k, v in Counter([c[match_column] for c in return_data]).items():
        print("- {} {}".format(k, format(v, ',')))
                                                   
    return return_data


def make_graph(nodes, edges):
    """
    Creates a networkx graph from a list of nodes [{'id':id, attribute1:'xyz',...}] and edges of the same formatting.
    The 'id' for each node is used as the node id within the graph and the rest is added in as attributes.
    Edges are only added if both nodes exist in the graph and all other information is added as an edge attribute
    """
    
    # Record keeping
    node_attr_count = 0
    edge_attr_count = 0
    dangling_edges = 0
    repeat_edges = 0
    
    # Initiate graph and add nodes
    kg_graph = nx.Graph()
    for node in nodes:
        
        # Add node
        kg_graph.add_node(node["id"])
        
        # Add attributes
        attributes = {node["id"]:{k:v for k,v in node.items() if k != "id"}}
        
        nx.set_node_attributes(kg_graph, attributes)
        node_attr_count += len(attributes[node["id"]])
    

    # Add edges
    for edge in edges:
        
        # Ensure both nodes of edge exist in graph
        e1, e2 = (edge["subject"], edge["object"])
        eid, eid_rev = (e1, e2), (e2, e1)
        if (not kg_graph.has_node(e1)) or (not kg_graph.has_node(e2)):
            dangling_edges += 1
            continue
        
        # Check for potential repeat edge
        elif kg_graph.has_edge(e1, e2):
            repeat_edges += 1
            
            #repeat_edges[eid].append(edge)
            #repeat_edges[eid_rev].append(edge)
            
            # The "repeat_edges" (subject, object) demonstrate the consolodation of different ontologies --> mondo
            #for r in repeat_edges[eid]:
            #    print(r["object"], r["subject"], r["predicate"], r["category"])
            #print("################")
            #print("################")
        
        else:
            ##repeat_edges.update({eid:[edge]})
            ##repeat_edges.update({eid_rev:[edge]})
            #print(repeat_edges[eid])
            
            # Add edge along with attributes
            kg_graph.add_edge(e1, e2)
            attributes = {eid:{k:v for k,v in edge.items()}}
            nx.set_edge_attributes(kg_graph, attributes)
            edge_attr_count += len(attributes[eid])
    
    
    # Remove any nodes with zero neighbors
    keep_nodes = set([node_id for node_id in kg_graph.nodes() if len(set(kg_graph.neighbors(node_id))) > 0])
    nsubgraph = kg_graph.subgraph(keep_nodes)
    
    print("- Graph created with {} nodes, and {} edges".format(format(nsubgraph.number_of_nodes(), ','),
                                                               format(nsubgraph.number_of_edges(), ',')))
    print("- Dangling edges found {}".format(format(dangling_edges, ',')))
    
    return nsubgraph



################################
### Pairwise comparison code ###
def generate_pairwise_comparison_args(sample_data):
    
    keys = list(sample_data.keys())
    comps = []
    ind = 0
    for i,k in enumerate(keys[0:-1]):
        p_ids = sample_data[k]
        for k2 in keys[i+1:]:
            p_ids2 = sample_data[k2]
            comps.append([k, k2, p_ids, p_ids2, ind])
            ind += 1
    
    print("- {} Pairwise comparison terms generated...".format(format(len(comps), ',')))
    print("- {} Labels generated...".format(format(len(keys), ',')))
    return comps, keys


def semsim_termset_compare(hp_db_path, input_args):
    
    # Load in data
    semsim_obj = Semsimian(spo=None, resource_path=hp_db_path)

    # Create output datastructure
    sim_metrics = ["jaccard_similarity", "ancestor_information_content", "phenodigm_score"]
    output_data = {"comp_name":[],
                   "comp_index":[],
                   "jaccard_similarity":[],
                   "ancestor_information_content":[],
                   "phenodigm_score":[]}
    
    # For each set of input arguments, compute jaccar, ic, and phenodigm sim scores and add to our output dict
    for i, argset in enumerate(input_args):
        
        termset1, termset2 = argset[2], argset[3]
        sname1, sname2 = argset[0], argset[1]
        comp_index = argset[4]
        output_data["comp_name"].append(tuple([sname1, sname2]))
        output_data["comp_index"].append(comp_index)
        
        for sim in sim_metrics:
            score = semsim_obj.termset_comparison(set(termset1), set(termset2), score_metric=sim)
            output_data[sim].append(score)

        if i % 1000 == 0:
            print("- processed {}/{}".format(i+1, format(len(input_args), ',')))

    return output_data


def merge_and_sort_parallel_output(output, labels, outpath=False):
    
    # Merge results from previous step into single data structure
    # Each element in the output list is a dictionary {sival:[pval,pval,...], }
    merged_data = {}
    for p in output:
        for k,v in p.items():
            if k not in merged_data:
                merged_data.update({k:[]})

            if type(v) != type([]):
                v = list(v)
            merged_data[k] += v
    
    for k,v in merged_data.items():
        print(k, len(v))
    # Sort all data based on the "comp_index"... 
    # This allows us to recapitulate the input "collapsed" pairwise comparison matrix
    for mkey in list(merged_data.keys()):
        print(mkey, "merging...")

        if mkey == "comp_index":
            continue
        
        # Must avoid combining strings with integers for array notation
        if mkey == "comp_name":
            cmat = sorted([[str(v),str(ind)] for v,ind in zip(merged_data[mkey], merged_data["comp_index"])], key=lambda x: int(x[1]))
        
        else:
            cmat = sorted([[v,ind] for v,ind in zip(merged_data[mkey], merged_data["comp_index"])], key=lambda x: x[1])
        
        cmat = np.asarray(cmat).T[0]
        merged_data[mkey] = cmat
    
    merged_data["comp_index"] = sorted(merged_data["comp_index"])
    merged_data.update({"labels":labels})
    if outpath != False:
        pickle.dump(merged_data, gzip.open(outpath, "wb"))
    
    return merged_data




if __name__ == '__main__':
    ################
	## ARG PARSE ###
    def parse_input_command():
        parser = argparse.ArgumentParser(description='Phenopacket Pairwise Similarity...')
        parser.add_argument("-p","--project_dir", help="Top most project directory", required=True, type=str, default=None)
        parser.add_argument("-c","--num_proc", help="number of cpu cores to use", required=False, type=int, default=2)
        return parser.parse_args()

    args = parse_input_command()
    ############################
    
    ###############
    ### PROGRAM ###

    # Data and results variables / paths
    nodes_path = os.path.join(args.project_dir, "monarch_kg", "monarch-kg_nodes.tsv")
    edges_path = os.path.join(args.project_dir, "monarch_kg", "monarch-kg_edges.tsv")
    hp_db_path = os.path.join(args.project_dir, "monarch_kg", "hp.db")
    res_path = os.path.join(args.project_dir, "results", "disease_pairwise_similarity_results.pkl.gz")

    # Initiate the biolink model categories we wish to pull for our nodes and edges
    valid_nodes = set(["biolink:Gene", 
                       "biolink:Disease", 
                       "biolink:PhenotypicFeature"])

    valid_edges = set(["biolink:DiseaseToPhenotypicFeatureAssociation", 
                       "biolink:CausalGeneToDiseaseAssociation"])

    # Initiate our graph from monarch nodes and edges .tsv files
    nodes = get_entities_from_file(nodes_path, select_data=valid_nodes, match_column='category')
    edges = get_entities_from_file(edges_path, select_data=valid_edges, match_column='category')
    graph = make_graph(nodes, edges)
    
    # Read graph data to get hypergeom testing params and datastructures
    dis_nodes = [n for n in graph.nodes() if graph.nodes[n]['category'] == 'biolink:Disease']
    phen_nodes = [n for n in graph.nodes() if graph.nodes[n]['category'] == 'biolink:PhenotypicFeature']
    
    # Sanity check to ensure only HP (human phenotype) terms have been pulled into the graph
    if len([ppp for ppp in phen_nodes if not ppp.startswith("HP:")]) > 0:
        sys.exit("- ERROR, Non HP terms found in graph associated with diseases. Exiting...")
    
    # Generate disease --> phenotype mapping 
    # Note, graph building process could be skipped by just selecting 
    dis_to_hp = {d:list(set([nn for nn in list(graph.neighbors(d))
                        if graph.nodes[nn]['category'] == "biolink:PhenotypicFeature"])) for d in dis_nodes}


    # Progress statement
    print("- Disease nodes present {}".format(format(len(dis_nodes), ',')))
    print("- Phenotype nodes present {}".format(format(len(phen_nodes), ',')))
    print("- Preprocessing complete...")


    # Subsample data (for testing)
    subsample_data = {k:dis_to_hp[k] for k in list(dis_to_hp.keys())[0:25]}

    comp_args, labs = generate_pairwise_comparison_args(subsample_data)
    #comp_args, labs = generate_pairwise_comparison_args(dis_to_hp)

    # Divide workload for parallel processing
    div_args = divide_workload(comp_args, num_proc=args.num_proc)
    div_semsim_objs = [copy.copy(hp_db_path) for i in range(0, args.num_proc)]

    # Setup parallel processing overhead, kick off jobs via asynchronous processing, and retrieve results
    #semsim_termset_compare(div_semsim_objs[0], div_args[0]+div_args[1])

    output = mp.Queue()
    pool = mp.Pool(processes=args.num_proc)
    results = [pool.apply_async(semsim_termset_compare, args=(sem, inargs,)) for sem, inargs in zip(div_semsim_objs, div_args)]
    output = [p.get() for p in results]
    pool.close()
    pool.join()

    print("- Semantic similarity comparisons completed...")
    print("- Combining and writing results to {}...".format(res_path))
    comp_data = merge_and_sort_parallel_output(output=output, labels=labs, outpath=res_path)