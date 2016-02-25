import argparse
from collections import defaultdict
import csv
import networkx as nx
import operator
import pybedtools
import sys


def get_node_size(node):
    return int(node[4])


def find_consensus_calls(graph):
    """
    Find all connected subgraphs such that calls with only self-self overlaps
    will exist as singleton graphs while all nodes that overlap each other
    directly or transitively will be clustered in the same graph.
    """
    for subgraph in nx.connected_component_subgraphs(graph):
        # Collect all nodes in this group by their start positions.
        nodes_by_start = defaultdict(list)
        for node in subgraph.nodes():
            try:
                nodes_by_start[int(node[1])].append(node)
            except:
                sys.stderr.write("Exception with node (%s) and graph (%s)\n and " % (node, graph))
                raise

        # Find the start position with the most nodes (i.e., the consensus start
        # position) and return one of the nodes from this set as the
        # representative node.
        consensus_nodes = max(nodes_by_start.values(), key=lambda nodes: len(nodes))

        # Report the representative node along with the number of nodes in the
        # subgraph corresponding to the support for the representative event.
        print "\t".join(consensus_nodes[0] + (str(len(subgraph.nodes())),))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("calls", help="BED file of calls from various callsets")
    parser.add_argument("--reciprocal_overlap", type=float, default=0.5, help="proportion of reciprocal overlap required to consider two calls the same")
    args = parser.parse_args()

    # Load all calls from the given BED file. Loading directly into a BedTool
    # from the filename is not guaranteed to work if the input BED file has the
    # same number of columns as a standard bedN format (e.g., bed12) and any of
    # those columns contain values of a type that differs from the standard.
    with open(args.calls, "r") as fh:
        reader = csv.reader(fh, delimiter="\t")
        calls = pybedtools.BedTool([row for row in reader])
        columns_per_call = len(list(calls[0]))

    # Intersect the given calls with themselves. Self-self overlaps will be
    # reported along with any other overlaps matching the given reciprocal
    # overlap proportion.
    intersected_calls = calls.intersect(b=calls, f=args.reciprocal_overlap, r=True, wao=True)

    # Create a graph connecting all calls that share a reciprocal overlap.
    current_contig = None
    for call in intersected_calls:
        if current_contig != call[0]:
            # If this isn't the first pass through the calls and we've found a
            # new contig, find the consensus calls and print them.
            if not current_contig is None:
                find_consensus_calls(graph)

            # If we've switched to a different contig, create a new graph for
            # that contig.
            current_contig = call[0]
            graph = nx.Graph()

        left = tuple(call[:columns_per_call])
        right = tuple(call[columns_per_call:-1])
        graph.add_edge(left, right)
    else:
        # If we've finished processing the last call, get consensus calls for
        # the final contig's graph.
        find_consensus_calls(graph)
