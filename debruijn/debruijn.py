#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

from itertools import combinations
import matplotlib.pyplot as plt
import textwrap
import statistics
from random import randint
import argparse
from ast import Delete
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
matplotlib.use("Agg")

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage="{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file (default contigs.fasta)")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as an image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file, 'r') as f:
        for i in f:
            yield next(f).strip()
            next(f)
            next(f)


def cut_kmer(read, kmer_size):
    for i in range(0, len(read)-kmer_size+1):
        yield str(read[i:i+kmer_size])


def build_kmer_dict(fastq_file, kmer_size):
    kmer_dict = {}
    for i in read_fastq(fastq_file):
        for kmer in cut_kmer(i, kmer_size):
            if kmer not in kmer_dict:
                kmer_dict[kmer] = 0
            kmer_dict[kmer] += 1
    return kmer_dict


def build_graph(kmer_dict):
    g = nx.DiGraph()
    for k, n in kmer_dict.items():
        g.add_edge(k[:-1], k[1:], weight=n)
    return g


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    graph.remove_edges_from(path_list)
    if delete_entry_node and delete_sink_node:
        l = list(set([item for t in path_list for item in t]))
        graph.remove_nodes_from(l)
    if delete_entry_node and not delete_sink_node:
        for i in path_list:
            graph.remove_nodes_from(i[:len(i)-1])
    if delete_sink_node and not delete_entry_node:
        for i in path_list:
            graph.remove_nodes_from(i[1:len(i)])
    if not (delete_entry_node and delete_sink_node):
        for i in path_list:
            graph.remove_nodes_from(i[1:len(i)-1])
    return graph


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    std_weight = statistics.stdev(weight_avg_list)
    std_length = statistics.stdev(path_length)
    if(std_weight > 0):
        path_list.pop(weight_avg_list.index(max(weight_avg_list)))
        remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    if(std_weight == 0):
        if(std_length > 0):
            path_list.pop(path_length.index(max(path_length)))
            remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
        if(std_length == 0):
            path_list.pop(randint(0, len(path_length)))
            remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    return graph


def path_average_weight(graph, path):
    """Compute the weight of a path"""
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])


def solve_bubble(graph, ancestor_node, descendant_node):
    path_list = []
    path_length = []
    weight_avg_list = []
    for path in nx.all_simple_paths(graph, ancestor_node, descendant_node):
        path_list.append(path)
        path_length.append(len(path))
        weight_avg_list.append(path_average_weight(graph, path))
    graph = select_best_path(graph, path_list, path_length, weight_avg_list)
    return graph


def simplify_bubbles(graph):
    bubble = False
    for n in graph.nodes():
        pred = [e for e in graph.predecessors(n)]
        if len(list(pred)) > 1:
            for comb in combinations(pred, 2):
                ancestor = nx.lowest_common_ancestor(graph, comb[0], comb[1])
                if ancestor != None:
                    bubble = True
                    break
        if bubble == True:
            break
    if bubble:
        graph = simplify_bubbles(solve_bubble(graph, ancestor, n))
    return graph


def solve_entry_tips(graph, starting_nodes):
    tip = False
    for n in graph.nodes():
        pred = [e for e in graph.predecessors(n)]
        if len(pred) > 1:
            if(len(list(set(pred).intersection(starting_nodes))) > 1):
                tip = True
                break
            if tip == True:
                break
    if tip == True:
        graph = solve_entry_tips(select_best_path())

    return graph


def solve_out_tips(graph, ending_nodes):
    pass


def get_starting_nodes(graph):
    return [n for n, d in graph.in_degree() if d == 0]


def get_sink_nodes(graph):
    return [u for u, d in graph.out_degree() if not d]


def get_contigs(graph, starting_nodes, ending_nodes):
    l = []
    s = ""
    for i in starting_nodes:
        for j in ending_nodes:
            for path in nx.all_simple_edge_paths(graph, i, j):
                for k in range(0, len(path), 2):
                    s += path[k][0]
            s += j
            l.append((s, len(s)))
            s = ""
    return l


def save_contigs(contigs_list, output_file):
    with open(output_file, "w") as f:
        count = 0
        for i in contigs_list:
            f.write(f">contig_{count} len={i[1]}\n")
            f.write(f"{i[0]}")
            count += 1
    f.close()


def draw_graph(graph, graphimg_file):
    """Draw the graph
    """
    fig, ax = plt.subplots()
    elarge = [(u, v)
              for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    # print(elarge)
    esmall = [(u, v)
              for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    # print(elarge)
    # Draw the graph with networkx
    # pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


# ==============================================================
# Main program
# ==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    #args = get_arguments()
    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)


if __name__ == '__main__':
    main()
