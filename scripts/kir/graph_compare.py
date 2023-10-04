#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import argparse
from pathlib import Path
from pafpy import PafFile
from gfapy import Gfa
import networkx as nx
import sys
from collections import defaultdict
import re
import matplotlib.colors as mcolors
COLORS = [c for c in mcolors.TABLEAU_COLORS.values()]


def rect(x, y, width, height, opacity=1, fill="black"):
    return '<rect x="{}" y="{}" width="{}" height="{}" fill-opacity="{}" fill="{}" />'.format(x, y, width, height, opacity, fill)


def text(x, y, size, text):
    return '<text x="{}" y="{}" font-family="monospace" font-size="{}">{}</text>'.format(x, y, size, text)


def line(x1, y1, x2, y2, color="black", width=1, opacity=1):
    return '<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" stroke-width="{}" stroke-opacity="{}" />'.format(x1, y1, x2, y2, color, width, opacity)


def path(x0, y0, x1, y1, x2, y2, x3, y3, color='#777', width=5, opacity=1):
    return '<path d="M {} {} C {} {}, {} {}, {} {}" stroke="{}" fill="transparent" stroke-width="{}" stroke-opacity="{}"/>'.format(x0, y0, x1, y1, x2, y2, x3, y3, color, width, opacity)


def poly(x0, y0, x1, y1, x2, y2, x3, y3, opacity=0.5, color="black"):
    return '<path d="M {} {} L {} {} L {} {} L {} {} z" fill="{}" fill-opacity="{}"/>'.format(x0, y0, x1, y1, x2, y2, x3, y3, color, opacity)


def join(x1, y1, x2, y2, start_to_right, end_to_right, dx=100, opacity=1.0):
    x1c = x1 + dx if start_to_right else x1 - dx
    y1c = y1
    x2c = x1 + dx if end_to_right else x2 - dx
    y2c = y2
    return path(x1, y1, x1c, y1c, x2c, y2c, x2, y2, opacity=opacity)


def svg(elements, width=800, height=800):
    head = '<svg width="{}" height="{}" xmlns="http://www.w3.org/2000/svg">'.format(
        width, height)
    tail = '</svg>'
    return head + '\n'.join(elements) + tail


def aligned_pairs(cigar, ignore_match=False):
    """
    given `cigar=":733*ga:12536*ga:21*cg:1013*ga:2138+t:6250+t:142*ct:1058+a:283*ct:2736+a:5472+t:4382+t:2454"`
    supported: :,*,+,-
    ignored: =,~

    >>> cigar = "-ac:2+tt:1*ac"
    >>> aligned_pairs(cigar)
    [(0, 0), (1, 0), (2, 0), (3, 1), (4, 2), (4, 3), (4, 4), (5, 5)]
    >>> aligned_pairs(cigar, ignore_match=True)
    [(0, 0), (1, 0), (4, 2), (4, 3), (5, 5)]
    """
    pairs = []
    # ref query
    i, j = 0, 0
    for segment in re.findall('\*[acgt][acgt]|:\d+|\+[acgt]+|-[acgt]+', cigar):
        type = segment[0]
        if type == ':':
            # :n
            n = int(segment[1:])
            for _ in range(n):
                if not ignore_match:
                    pairs.append((i, j))
                i += 1
                j += 1
        elif type == '*':
            # *xy
            x = segment[1]
            y = segment[2]
            pairs.append((i, j))
            i += 1
            j += 1
        elif type == '+':
            # +s (insertion)
            s = segment[1:]
            n = len(s)
            for _ in range(n):
                pairs.append((i, j))
                j += 1
        elif type == '-':
            # -s (deletion)
            s = segment[1:]
            n = len(s)
            for _ in range(n):
                pairs.append((i, j))
                i += 1
    return pairs


def reverse(strand):
    if strand == '+':
        return '-'
    elif strand == '-':
        return '+'
    else:
        raise Exception()


def parse_gfa(filename):
    graph = nx.DiGraph()
    seqs = dict()
    with open(filename) as f:
        for line in f:
            segments = line.split('\t')
            if segments[0] == 'S':
                name = segments[1]
                seq = segments[2]

                graph.add_node((name, '+'))
                graph.add_node((name, '-'))
                seqs[name] = len(seq)
            elif segments[0] == 'L':
                name_a = segments[1]
                strand_a = segments[2]
                name_b = segments[3]
                strand_b = segments[4]
                # a -> b
                graph.add_edge((name_a, strand_a), (name_b, strand_b))
                # b -> a
                graph.add_edge((name_b, reverse(strand_b)),
                               (name_a, reverse(strand_a)))
    return graph, seqs


def parse_paf(filename):
    match = defaultdict(list)
    haps = dict()
    seqs = dict()
    with PafFile(filename) as paf:
        for record in paf:
            # nm = record.get_tag("NM").value
            match[record.qname].append(record)
            haps[record.tname] = record.tlen
            seqs[record.qname] = record.qlen
    return match, haps, seqs


def eprint(*args):
    print(*args, file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description='generate SVG')
    parser.add_argument('gfa', type=Path, help='GFA of asm')
    parser.add_argument('paf', type=Path, help='PAF of genome vs asm')
    parser.add_argument('--width', type=int,
                        default=1000, help='svg width')
    parser.add_argument('--margin', type=int, default=30, help='')
    parser.add_argument('--box_height', type=int, default=20, help='')
    parser.add_argument('--min_length', type=int, default=1000, help='')
    parser.add_argument('--order', type=str, nargs='+')
    parser.add_argument('--draw_mismatch_threshold', type=int, default=100,
                        help='if the number of mismatches is below this threshold, visualize mismatch position in red line. To disable mismatch visualization, set threshold to zero.')
    args = parser.parse_args()

    graph, seqs = parse_gfa(args.gfa)
    eprint(graph.nodes)

    match, haps, seqs_from_paf = parse_paf(args.paf)
    eprint(match)
    eprint(haps)
    n_haps = len(haps)
    assert n_haps == 2, "genome is not diploid"

    for seqname, seqlen in seqs_from_paf.items():
        if seqs[seqname] != seqlen:
            print('seq "{}" length not match gfa={} paf={} using paf'.format(
                seqname, seqs[seqname], seqlen), file=sys.stderr)
            seqs[seqname] = seqlen

    # debug output
    hapnames = sorted(haps.keys())
    for i, hapname in enumerate(hapnames):
        print('HAP', i, hapname, haps[hapname], file=sys.stderr)

    def bestmatch(records):
        """convert List[paf records] into best hit (hap, pos, strand) (in terms of match length)"""
        if len(records) == 0:
            return (hapnames[0], 0, '+')
        record = max(records, key=lambda record: record.mlen)
        return (record.tname, record.tstart, str(record.strand))

    seqpositions = {seqname: bestmatch(match[seqname])
                    for seqname in seqs.keys()}
    seqnames = sorted([seqname for seqname in seqs.keys()],
                      key=lambda seqname: seqpositions[seqname])
    for i, seqname in enumerate(seqnames):
        print('SEQ', i, seqname, seqpositions[seqname], file=sys.stderr)

    if args.order:
        print('ORDER', args.order, file=sys.stderr)
        seqnames = []
        for order in args.order:
            segs = order.split(',')
            seqname = segs[0]
            if len(segs) == 3:
                # "seqname,position,strand"
                position = int(segs[1])
                strand = segs[2]
                assert strand == '+' or strand == '-'
                seqpositions[seqname] = (hapnames[0], position, strand)
            seqnames.append(seqname)

    # bp vs pixel conversion
    width = args.width
    margin = args.margin
    box_height = args.box_height
    bp_per_px = max(haps.values()) / (width - 2 * margin)
    font_size = 20
    text_margin = 3
    px_per_row = margin + box_height
    # raw position function
    def x(bp): return bp / bp_per_px
    def y(n): return margin + n * px_per_row
    n_seqs = len(seqnames)
    height = y(n_haps + n_seqs)

    def hap_to_y(hapname):
        if hapname == hapnames[0]:
            return y(0)
        else:
            return y(n_seqs + 1)

    def seq_to_y(seqname):
        i = next((i for i, s in enumerate(seqnames) if s == seqname))
        return y(i + 1)

    def seq_to_x(seqname):
        hapname, start, strand = seqpositions[seqname]
        length = max(args.min_length, seqs[seqname])
        left = margin + x(start)
        right = margin + x(start) + x(length)
        return (left, right)

    print('<svg width="{}" height="{}" xmlns="http://www.w3.org/2000/svg">'.format(width, height))

    for hapname in hapnames:
        x_hap = margin
        y_hap = hap_to_y(hapname)
        print(text(x_hap, y_hap - text_margin, font_size,
              "{} {}bp".format(hapname, haps[hapname])))
        print(rect(x_hap, y_hap, x(haps[hapname]), box_height, fill="#555"))

    for i, seqname in enumerate(seqnames):
        color = COLORS[i % len(COLORS)]
        hapname, start, strand = seqpositions[seqname]
        y_seq = seq_to_y(seqname)
        x_seq_left, x_seq_right = seq_to_x(seqname)
        x_seq_width = x_seq_right - x_seq_left
        print(text(x_seq_left, y_seq - text_margin,
              font_size, "{}{} {}bp".format(seqname, strand, seqs[seqname])))
        print(rect(x_seq_left, y_seq, x_seq_width, box_height, fill=color))

        # show alignment
        for m in match[seqname]:
            hapname = m.tname
            mapq = m.mapq

            y_hap = hap_to_y(hapname)
            y_seq = seq_to_y(seqname)
            if hapname == hapnames[0]:
                y_hap += box_height
            else:
                y_seq += box_height

            x_seq_start = x_seq_left + x(m.qstart)
            x_hap_start = margin + x(m.tstart)
            print(line(x_seq_start, y_seq, x_hap_start, y_hap))

            x_seq_end = x_seq_left + x(m.qend)
            x_hap_end = margin + x(m.tend)
            print(line(x_seq_end, y_seq, x_hap_end, y_hap))

            print(poly(x_seq_start, y_seq, x_seq_end, y_seq,
                       x_hap_end, y_hap, x_hap_start, y_hap,
                       color=color if mapq > 0 else '#bbb',
                       opacity=0.2 if mapq > 0 else 0.2))

            nm = m.get_tag("NM").value
            cigar = m.get_tag("cs").value
            if nm <= args.draw_mismatch_threshold:
                for (tindex, qindex) in aligned_pairs(cigar, ignore_match=True):
                    x_hap = margin + x(m.tstart + tindex)
                    x_seq = x_seq_left + x(m.qstart + qindex)
                    print(line(x_seq, y_seq, x_hap, y_hap,
                          color="red", opacity=0.5))

        # show graph
        for edge in graph.out_edges((seqname, strand)):
            _, (seqname_child, strand_child) = edge
            print("EDGE", seqname, strand, seqname_child,
                  strand_child, file=sys.stderr)
            if seqname_child not in seqnames:
                continue
            # parent
            start_to_right = strand == '+'
            x_s = x_seq_right if start_to_right else x_seq_left
            y_s = seq_to_y(seqname) + box_height / 2
            # child
            _, _, strand_layout_child = seqpositions[seqname_child]
            x_seq_left_child, x_seq_right_child = seq_to_x(seqname_child)
            end_to_right = strand_child != strand_layout_child
            x_t = x_seq_right_child if end_to_right else x_seq_left_child
            y_t = seq_to_y(seqname_child) + box_height / 2
            # draw curve line
            print(join(x_s, y_s, x_t, y_t, start_to_right, end_to_right, opacity=0.5))
    print('</svg>')


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    main()
