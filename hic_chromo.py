import maxflow
from cooler_extended_arms import CoolerExtendedArms
import math
import pandas as pd
import argparse

def check_dict(d, item):
    if item not in d:
        d[item] = len(d)
    return d[item]


def load_network(c, chr_number, patient, sv_filepath):
    chr_number = str(chr_number)
    sv_master_table = pd.read_csv(sv_filepath)
    sv = sv_master_table[(sv_master_table["chrom1"] == chr_number) &
                                         (sv_master_table["chrom2"] == chr_number) &
                                        (sv_master_table["unique_id"] == patient)]
    sv = sv.drop(columns=['unique_id', 'chrom2', 'chrom1'])
    network_init = {}
    vertex_dict = {}
    break_set = set(sv.start1.unique())
    break_set.update(set(sv.start2.unique()))
    for breakpoint_1 in break_set:
        for breakpoint_2 in break_set:
            if breakpoint_1 == breakpoint_2:
                continue
            value = c.get_single_hic_score(breakpoint_1, breakpoint_2, chr_number)
            if math.isnan(value):
                continue
            v1 = check_dict(vertex_dict, breakpoint_1)
            v2 = check_dict(vertex_dict, breakpoint_2)
            network_init[(str(v1), str(v2))] = value*value*value
    return network_init, len(vertex_dict), vertex_dict


def WFind_Densest_Subgraph(cluster_threshold_size, c, chr_number, patient, sv_filepath):
    ''' This function performs the binary search of the density of subgraph and finds the densest subgraph.'''
    min_degree = 0

    sum_w = 0.0
    sum_ = 0
    network, number_of_nodes, vertex_dict = load_network(c, chr_number, patient, sv_filepath)
    for from_node, to_node in network:
        w = network[(from_node, to_node)]
        sum_w += w
        sum_ += 1
    number_of_edges = number_of_nodes * number_of_nodes
    max_degree = sum_w

    subgraph = []
    if sum_w == 0:
        return subgraph
    difference = (1.0 / (number_of_nodes * (number_of_nodes + 1))) * (sum_w / sum_)
    while (max_degree - min_degree >= difference):
        least_density = (max_degree + min_degree) / 2.0
        source_segment = Wmake_graph(number_of_nodes, number_of_edges, least_density, network)
        if len(source_segment) <= cluster_threshold_size:
            max_degree = least_density
        else:
            min_degree = least_density
            subgraph = source_segment
    break_points = list()
    for node in subgraph:
        break_points.append(list(vertex_dict.keys())[list(vertex_dict.values()).index(node)])
    return break_points


def Wmake_graph(number_of_nodes, number_of_edges, least_density, network):
    ''' Constructs the network as per the specifications given by Goldberg'''
    graph = maxflow.Graph[float](number_of_nodes, number_of_edges)
    nodes = graph.add_nodes(number_of_nodes)
    degrees = {}
    sum_w = 0.0
    for from_node, to_node in network:
        w = network[(from_node, to_node)]
        sum_w += w
        graph.add_edge(nodes[int(from_node)], nodes[int(to_node)], w, w)

        if from_node in degrees:
            degrees[from_node] += w
        else:
            degrees[from_node] = w
        if to_node in degrees:
            degrees[to_node] += w
        else:
            degrees[to_node] = w
    for i in range(number_of_nodes):
        if str(i) not in degrees:
            degrees[str(i)] = 0
        graph.add_tedge(nodes[i], sum_w, sum_w + 2 * least_density - degrees[str(i)])
    source_segment = []
    '''Computes the max-flow in the graph'''
    max_flow = graph.maxflow()
    '''The following section of code finds which node belongs to which cutset.'''
    for i in nodes:
        if (graph.get_segment(nodes[i]) == 0):
            source_segment.append(nodes[i])
    return source_segment


def detect(resolution, hic_filepath, sv_filepath, out_filepath):
    c = CoolerExtendedArms(hic_filepath)
    df = pd.read_csv(sv_filepath)
    patients_set = set(df.unique_id.unique())
    result = list()
    for patient in patients_set:
        patient_dict = {}
        chromosomes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X"]
        for cr_number in chromosomes:
            answer = WFind_Densest_Subgraph(10, c, cr_number, patient, sv_filepath)
            patient_dict['patient_id'] = patient
            patient_dict[f"{cr_number}chr"] = answer
        result.append(pd.Series(patient_dict))
    pd.DataFrame(result,
                 columns=['patient_id', '1chr', '2chr', '3chr', '4chr', '5chr', '6chr', '7chr', '8chr', '9chr', '10chr',
                          '11chr', '12chr', '13chr', '14chr', '15chr', '16chr', '17chr', '18chr', '19chr', '20chr',
                          '21chr', '22chr', 'Xchr']).to_csv(out_filepath)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='This utilite performs chromothripsis structural variations detection by Hi-C data')
    parser.add_argument("--sv", required=True, type=str, help="Path to SVs file, .csv format")
    parser.add_argument("--hic", required=True, type=str, help="Path to HiC file - .cool or .mcool")
    parser.add_argument("--r", default=400000, type=int, help="Resolution of Hi-C file, required if Hi-C is .mcool file")
    parser.add_argument("--o", required=True, type=str, help="Output file name, .csv format")

    args = parser.parse_args()
    hic = args.hic
    if ".mcool" in hic:
        hic = f"{hic}::/resolutions/{args.r}"

    detect(args.r, hic, args.sv, args.o)
