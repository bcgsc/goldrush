#!/usr/bin/env python3
'''
Detect hairpin structure artifacts in oxford nanopore reads
'''
import argparse
from collections import defaultdict, namedtuple, Counter
import igraph as ig
import ntlink_utils
import sys

MinimizerInfo = namedtuple("MinimizerInfo", ["pos", "strand", "time_seen"])


def tally_sequence_lengths(index_filename):
    "Tally the sequence lengths from the faidx index file"
    lengths = {}
    with open(index_filename, 'r') as fin:
        for line in fin:
            name, length = line.strip().split("\t")[:2]
            lengths[name] = int(length)
    return lengths

def is_valid_mx(list_mx_info):
    "Return true if the given minimizer has a multiplicity of 2 and those are found from different strands"
    if len(list_mx_info) != 2:
        return False
    info1, info2 = list_mx_info
    if info1.strand == info2.strand:
        return False
    return True


def filter_ordered_sketch(mx_line):
    "Given a line of a indexlr file, parse and only keep minimizers with multiplicity of 2, strictly on different strands"

    mx_info = defaultdict()  # mx -> [Minimizer info]
    line = mx_line.strip().split("\t")
    if len(line) > 1:
        name, mxs_all = line
        mx_pos_strands = mxs_all.split(" ")
        for mx_pos_strand in mx_pos_strands:
            mx, pos, strand = mx_pos_strand.split(":")
            if mx not in mx_info:
                mx_info[mx] = [MinimizerInfo(int(pos), strand, 1)] #!! TODO: change numbers?
            else:
                mx_info[mx].append(MinimizerInfo(int(pos), strand, 2))

        mx_info = {mx: mx_info[mx] for mx in mx_info if is_valid_mx(mx_info[mx])}
        mxs = [mx_pos_strand.split(":")[0] for mx_pos_strand in mx_pos_strands if mx_pos_strand.split(":")[0] in mx_info]
    #print(mx_info)
    #print(mxs)
    return mx_info, mxs

def set_edge_attributes(graph, edge_attributes): #!! TODO: from ntJoin code
    "Sets the edge attributes for a python-igraph graph"
    graph.es()["weight"] = [edge_attributes[e]['weight'] for e in sorted(edge_attributes.keys())]

def build_graph(list_mxs, mx_info):
    "Builds an undirected graph: nodes=minimizers; edges=between adjacent minimizers"
    graph = ig.Graph()

    vertices = set()
    edges = defaultdict(dict)  # source -> target -> weight

    for i, j in zip(range(0, len(list_mxs)),
                    range(1, len(list_mxs))):
        if list_mxs[i] in edges and \
                list_mxs[j] in edges[list_mxs[i]]:
            edges[list_mxs[i]][list_mxs[j]] += 1
        elif list_mxs[j] in edges and \
                list_mxs[i] in edges[list_mxs[j]]:
            edges[list_mxs[j]][list_mxs[i]] += 1
        else:
            edges[list_mxs[i]][list_mxs[j]] = 1
        vertices.add(list_mxs[i])
    if list_mxs:
        vertices.add(list_mxs[-1])

    formatted_edges = [(s, t) for s in edges for t in edges[s]]

    graph.add_vertices(list(vertices))
    graph.add_edges(formatted_edges)

    edge_attributes = {ntlink_utils.edge_index(graph, s, t): {"weight": edges[s][t]}
                       for s in edges for t in edges[s]}
    set_edge_attributes(graph, edge_attributes)

    return graph

def filter_graph_global(graph): #!! TODO: From ntJoin
    "Filter the graph globally based on minimum edge weight"
    to_remove_edges = [edge.index for edge in graph.es()
                       if edge['weight'] < 2]
    new_graph = graph.copy()
    new_graph.delete_edges(to_remove_edges)
    return new_graph


def is_graph_linear(graph): #!! TODO: From ntJoin
    "Given a graph, return True if all the components are linear"
    for component in graph.components():
        component_graph = graph.subgraph(component)
        if not all(u.degree() < 3 for u in component_graph.vs()):
            return False
    return True

def determine_source_vertex(source_nodes, graph, mx_info):
    min_pos_sources = [(node, min(mx_info_entry.pos for mx_info_entry in mx_info[ntlink_utils.vertex_name(graph, node)]))
                       for node in source_nodes]
    source_target = sorted(min_pos_sources, key=lambda x:x[1])
    return source_target[0][0], source_target[1][0]

def get_positions(path, mx_info, graph, time_seen):
    "Return positions in ordered path for time_seen specified"
    pos_path = []
    for mx in path:
        for info in mx_info[ntlink_utils.vertex_name(graph, mx)]:
            if info.time_seen == time_seen:
                pos_path.append(info.pos)
    return pos_path

def paths_largely_increasing_or_decreasing(paths):
    "Return true if the positions in the paths are largely increasing or decreasing"
    consistent_positions = True
    print(paths)
    for time_seen in paths:
        path = paths[time_seen]
        if all(x < y for x, y in zip(path, path[1:])):
            continue
        if all(x > y for x, y in zip(path, path[1:])):
            continue
        print([x < y for x, y in zip(path, path[1:])])
        tally = Counter([x < y for x, y in zip(path, path[1:])])
        positive_perc = tally[True] / float(len(path) - 1) * 100
        negative_perc = 100 - positive_perc
        print(path, positive_perc)
        if positive_perc >= 95 or negative_perc >= 95: #!! TODO: magic number
            continue
        consistent_positions = False  # If make it here, none of criteria were met

    return consistent_positions


def parse_mx_path(path, mx_info, graph):
    "Parse a path of mx to find the mapped regions per arm of putative hairpin"
    pos_paths = {1: get_positions(path, mx_info, graph, time_seen=1),
                2: get_positions(path, mx_info, graph, time_seen=2)}
    if paths_largely_increasing_or_decreasing(pos_paths):
        mapped_extents = [(min(pos_paths[time_seen][0], pos_paths[time_seen][-1]),
                           max(pos_paths[time_seen][0], pos_paths[time_seen][-1]))
                          for time_seen in pos_paths]
        return mapped_extents
    return None



def find_paths(graph, mx_info): #!! TODO Snippet adapted from ntJoin
    "Find paths/mapped segments through the graph"
    mapped_regions = []
    for subcomponent in graph.components():
        subcomponent_graph = graph.subgraph(subcomponent)
        source_nodes = [node.index for node in subcomponent_graph.vs() if node.degree() == 1]
        if len(source_nodes) == 2:
            #print(source_nodes)
            source, target = determine_source_vertex(source_nodes, subcomponent_graph, mx_info)
            #print(source, target)
            path = subcomponent_graph.get_shortest_paths(source, target)[0]
            num_edges = len(path) - 1
            if len(path) == len(subcomponent_graph.vs()) and \
                    num_edges == len(subcomponent_graph.es()) and len(path) == len(set(path)):
                # All the nodes/edges from the graph are in the simple path, no repeated nodes
                mx_regions = parse_mx_path(path, mx_info, subcomponent_graph)
                if mx_regions is not None:
                    mapped_regions.append(mx_regions)
    #print(mapped_regions)
    return mapped_regions

def find_max_covered_mapped_regions(mapped_regions):
    "Given mapped regions, return the maximum bases covered by a pair"
    max_bp = 0
    print(mapped_regions)
    total_bases_covered = sum(region[1] - region[0] for regions in mapped_regions for region in regions)
    print(total_bases_covered)
    for region_sets in mapped_regions:
        bp_covered = sum([region[1] - region[0] for region in region_sets])
        if bp_covered > max_bp:
            max_bp = bp_covered
    print(max_bp)
    return max_bp

def print_graph(graph, list_mx_info, prefix): ## !! TODO: for troubleshooting
    "Prints the minimizer graph in dot format"
    out_graph = prefix + ".mx.dot"
    outfile = open(out_graph, 'a')

    outfile.write("graph G {\n")

    colours = ["red", "green", "blue", "purple", "orange",
               "turquoise", "pink", "yellow", "orchid", "salmon"]
    for node in graph.vs():
        mx_ctg_pos_labels = str(list_mx_info[node['name']])
        node_label = "\"%s\" [label=\"%s\n%s\"]" % (node['name'], node['name'], mx_ctg_pos_labels)
        outfile.write("%s\n" % node_label)

    for edge in graph.es():
        outfile.write("\"%s\" -- \"%s\"" %
                      (ntlink_utils.vertex_name(graph, edge.source),
                       ntlink_utils.vertex_name(graph, edge.target)))
        weight = edge['weight']
        if weight == 1:
            colour = colours[0]
        elif weight == 2:
            colour = "lightgrey"
        else:
            colour = "black"
        outfile.write(" [weight=%s color=%s]\n" % (weight, colour))

    outfile.write("}\n")

def get_mx_pos(mx_data, time_seen=1):
    "Return position of mx with desired time_seen"
    for mx_info_el in mx_data:
        if mx_info_el.time_seen == time_seen:
            return mx_info_el.pos

def filter_branches(graph, mx_info):
    "For each branch node, only keep two edges most consistent in terms of positions"
    branch_nodes = [node.index for node in graph.vs() if node.degree() > 2]
    to_remove_edges = []
    for node in branch_nodes:
        node_edges = []
        for edge in graph.incident(node):
            source, target = graph.es()[edge].source, graph.es()[edge].target
            source, target = ntlink_utils.vertex_name(graph, source), \
                             ntlink_utils.vertex_name(graph, target)
            if source != target:
                total_dist = abs(get_mx_pos(mx_info[source], time_seen=1) - get_mx_pos(mx_info[target], time_seen=1))
                total_dist += abs(get_mx_pos(mx_info[source], time_seen=2) - get_mx_pos(mx_info[target], time_seen=2))
                node_edges.append((edge, total_dist))
            else:
                # source == target, if self edge
                to_remove_edges.append(edge)
        if len(node_edges) > 2:
            add_remove_edges = sorted(node_edges, key=lambda x:x[1], reverse=True)[:-2]
            to_remove_edges.extend([rm_edge[0] for rm_edge in add_remove_edges])

    new_graph = graph.copy()
    new_graph.delete_edges(to_remove_edges)
    return new_graph

def detect_hairpins(args, seq_lengths):
    "Read through minimizers for each read, and determine whether it is a putative hairpin artifact"
    hairpins = 0
    total_reads = 0

    print("Name", "Status", sep="\t", file=sys.stderr)   

    with open(args.MX, 'r') as mx_in:
        for mx_line in mx_in:
            name, _ = mx_line.strip().split("\t")
            mx_info, mxs = filter_ordered_sketch(mx_line)
            graph = build_graph(mxs, mx_info)
            print_graph(graph, mx_info, "test_before")
            graph = filter_branches(graph, mx_info)
            print_graph(graph, mx_info, "test_branch")
            #graph = filter_graph_global(graph)
            #print_graph(graph, mx_info, "test")
            if is_graph_linear(graph): #!! TODO: add more sophisticated filter
                #print("HERE - linear")
                mapped_regions = find_paths(graph, mx_info)
                max_covered = find_max_covered_mapped_regions(mapped_regions)
                #print(max_covered)
                if max_covered/seq_lengths[name]*100 >= args.perc:
                    print(name, "Hairpin", sep="\t", file=sys.stderr)
                    hairpins += 1
                else:
                    print(name, "Non-hairpin", sep="\t", file=sys.stderr)
            total_reads += 1
    return hairpins, total_reads


def main():
    "Detect hairpin structures in input nanopore reads from minimizer sketches"
    parser = argparse.ArgumentParser(description="Detect hairpin artifacts in nanopore reads")
    parser.add_argument("MX", help="Input minimizers TSV file, or '-' if piping to standard in")
    parser.add_argument("-i", "--index", help="samtools faidx index for input reads", required=True, type=str)
    parser.add_argument("--perc", help="Percentage of read", type=float, default=90)
    args = parser.parse_args()

    args.MX = "/dev/stdin" if args.MX == "-" else args.MX

    seq_lengths = tally_sequence_lengths(args.index)
    hairpins, total_reads = detect_hairpins(args, seq_lengths)
    print("Total reads analyzed:", total_reads)
    print("Total hairpins identified:", hairpins)

if __name__ == "__main__":
    main()
