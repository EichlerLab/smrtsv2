import argparse
from collections import defaultdict
import csv
import intervaltree
import networkx as nx
import operator
import pybedtools
import sys

# Constants related to zero-based indexing of fields from the SV call format.
CHROMOSOME = 0
START = 1
END = 2
QUERY_NAME = 7
QUERY_START = 10
QUERY_END = 11
QUERY_LENGTH = 13

def get_node_size(node):
    return int(node[4])


def find_consensus_calls(graph, tiling_path):
    """
    Find all connected subgraphs such that calls with only self-self overlaps
    will exist as singleton graphs while all nodes that overlap each other
    directly or transitively will be clustered in the same graph.

    The consensus node has to be from the corresponding contig in the given
    chromosome's tiling path IntervalTree instance.
    """
    for subgraph in nx.connected_component_subgraphs(graph):
        # Find a node whose variant originates from the corresponding tiling
        # path contig.
        consensus_node = None
        for node in subgraph.nodes():
            # Search for a variant's start position in the tiling path.
            overlaps = tiling_path.search(int(node[START]))

            # Output all overlaps for the given start position where the contig
            # name is the same in the variant as the tiling path.
            for overlap in overlaps:
                if overlap[2] == node[QUERY_NAME]:
                    consensus_node = node
                    break

            if consensus_node is not None:
                break

        # Report the representative node along with the number of nodes in the
        # subgraph corresponding to the support for the representative event.
        if consensus_node is not None:
            print "\t".join(consensus_node + (str(len(subgraph.nodes())),))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("calls", help="BED file of calls from various callsets")
    parser.add_argument("tiling_path", help="BED file of tiling path through local assembly contigs. The fourth column is the contig name.")
    parser.add_argument("action", choices=("intersect", "window"), help="type of comparison to use to identify overlapping calls")
    parser.add_argument("--reciprocal_overlap", type=float, default=0.5, help="proportion of reciprocal overlap required to consider two calls the same")
    parser.add_argument("--window", type=int, default=100, help="inspect a window on either side of each call position to determine overlap")
    args = parser.parse_args()

    # Load all calls from the given BED file. Loading directly into a BedTool
    # from the filename is not guaranteed to work if the input BED file has the
    # same number of columns as a standard bedN format (e.g., bed12) and any of
    # those columns contain values of a type that differs from the standard.
    with open(args.calls, "r") as fh:
        reader = csv.reader(fh, delimiter="\t")
        calls = pybedtools.BedTool([row for row in reader])
        columns_per_call = len(list(calls[0]))

    if args.action == "intersect":
        # Intersect the given calls with themselves. Self-self overlaps will be
        # reported along with any other overlaps matching the given reciprocal
        # overlap proportion.
        intersected_calls = calls.intersect(b=calls, f=args.reciprocal_overlap, r=True, wao=True)
    else:
        # Inspect a window on either side of each call to determine
        # "overlap". This approach is useful for insertions which have single
        # base pair lengths, but may be effectively overlapping if close enough
        # together.
        intersected_calls = calls.window(b=calls, w=args.window)

    # Load the tiling path.
    chrom_tiling = {}
    with open(args.tiling_path, "r") as fh:
        for line in fh:
            vals = line.strip().split()
            if vals[0] not in chrom_tiling:
                chrom_tiling[vals[0]] = intervaltree.IntervalTree()

            # Store the interval of the current tile by chromosome along with
            # the name of the contig composing the tile.
            chrom_tiling[vals[0]].addi(int(vals[1]), int(vals[2]), vals[3])

    # Create a graph connecting all calls that share a reciprocal overlap.
    current_contig = None
    for call in intersected_calls:
        # Skip calls from chromosomes that are not in the tiling path.
        if call[0] not in chrom_tiling:
            continue

        if current_contig != call[0]:
            # If this isn't the first pass through the calls and we've found a
            # new contig, find the consensus calls and print them.
            if not current_contig is None:
                find_consensus_calls(graph, chrom_tiling[current_contig])

            # If we've switched to a different contig, create a new graph for
            # that contig.
            current_contig = call[0]
            graph = nx.Graph()

        left = tuple(call[:columns_per_call])

        # If we're intersecting, omit the final column that contains the total
        # bases overlapping between inputs. Otherwise, we're windowing and need
        # to keep the entire output.
        if args.action == "intersect":
            right = tuple(call[columns_per_call:-1])
        else:
            right = tuple(call[columns_per_call:])

        graph.add_edge(left, right)
    else:
        # If we've finished processing the last call, get consensus calls for
        # the final contig's graph.
        find_consensus_calls(graph, chrom_tiling[current_contig])
