import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import operator
import re


def make_windows(length, window, slide):
    """
    For a given length, return an iterator for intervals of length `window` with
    a slide of `slide`.

    >>> list(make_windows(8, 4, 0))
    [(0, 4), (4, 8)]
    >>> list(make_windows(8, 5, 0))
    [(0, 5), (5, 8)]
    >>> list(make_windows(8, 8, 0))
    [(0, 8)]
    >>> list(make_windows(8, 4, 2))
    [(0, 4), (2, 6), (4, 8)]
    >>> list(make_windows(8, 5, 2))
    [(0, 5), (2, 7), (4, 8)]
    >>> list(make_windows(7, 8, 0))
    [(0, 7)]
    """

    if slide == 0:
        windows = range(0, length, window)
    else:
        windows = range(0, length, slide)

    for start in windows:
        yield (start, min(start + window, length))

        # At most, only output one window at the end of the sequence.
        if length <= start + window:
            break


def fragment_sequence(sequence, window, slide=0):
    """Fragment a given sequence to the requested window length without a slide.

    >>> fragment_sequence("ACTGACTG", 4, 0)
    ['ACTG', 'ACTG']
    >>> fragment_sequence("ACTGACTG", 5, 0)
    ['ACTGA', 'CTG']
    >>> fragment_sequence("ACTGACTG", 8, 0)
    ['ACTGACTG']

    Fragment a given sequence to the requested window length with a slide.

    >>> fragment_sequence("ACTGACTG", 4, 2)
    ['ACTG', 'TGAC', 'ACTG']
    >>> fragment_sequence("ACTGACTG", 5, 2)
    ['ACTGA', 'TGACT', 'ACTG']

    Remove gap bases from input sequence and return the longest non-gap
    fragment. Don't return any sequence if the entire input is gap bases.

    >>> fragment_sequence("NNNNNNNN", 4, 2)
    []
    >>> fragment_sequence("ACTGNNNN", 4, 2)
    ['ACTG']
    >>> fragment_sequence("ACTGNNTA", 4, 2)
    ['ACTG']
    >>> fragment_sequence("ACNNACTA", 4, 2)
    ['ACTA']
    """

    # Check sequence for gap bases and keep the longest of the non-gap pieces in
    # the sequence.
    sequences = []
    sequence_pieces = [(piece, len(piece)) for piece in re.split("N+", sequence) if len(piece) > 0]

    if len(sequence_pieces) > 0:
        sorted_sequence_pieces = sorted(sequence_pieces, key=operator.itemgetter(1), reverse=True)
        sequence = sorted_sequence_pieces[0][0]
        sequence_length = sorted_sequence_pieces[0][1]

        if sequence_length > window:
            # Split sequence into two or more reads of the given length.
            window_ranges = make_windows(sequence_length, window, slide)
            sequences = [sequence[start:end] for start, end in window_ranges]
        else:
            # Output the sequence as is.
            sequences = [sequence]

    return sequences


if __name__ == "__main__":

    # Parse arguments
    parser = argparse.ArgumentParser(
        description='Create sequence fragments with a sliding-window over input sequences.')

    parser.add_argument('input', help='FASTA sequences to fragment.')
    parser.add_argument('output', help='Fragmented FASTA sequences.')
    parser.add_argument('window', type=int, help='Length to fragment each output fragment.')
    parser.add_argument('--slide', type=int, default=0, help='Number of bases to move when sliding fragment windows.')

    args = parser.parse_args()

    with open(args.output, 'w') as oh:
        for seq_record in SeqIO.parse(args.input, 'fasta'):
            sequences = fragment_sequence(str(seq_record.seq), args.window, args.slide)

            records = []

            for i in range(len(sequences)):
                records.append(SeqRecord(
                    Seq(sequences[i]),
                    id='{:s}_{:d}'.format(seq_record.id, i),
                    description=''
                ))

            SeqIO.write(records, oh, 'fasta')
