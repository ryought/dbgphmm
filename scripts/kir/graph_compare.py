#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import argparse
from pathlib import Path
from pafpy import PafFile
from Bio import SeqIO
import networkx as nx
import sys
from collections import defaultdict
import re
import matplotlib.colors as mcolors
from scipy.optimize import minimize
COLORS = [c for c in mcolors.TABLEAU_COLORS.values()]


def rect(x, y, width, height, opacity=1, fill="black", stroke="black", stroke_width=0):
    return '<rect x="{}" y="{}" width="{}" height="{}" fill-opacity="{}" fill="{}" stroke="{}" stroke-width="{}" />'.format(x, y, width, height, opacity, fill, stroke, stroke_width)


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


def parse_gfa(filename, ignore_n=True, ignore_zero_copy=True):
    graph = nx.DiGraph()
    seqs = dict()
    with open(filename) as f:
        for line in f:
            segments = line.split('\t')
            if segments[0] == 'S':
                name = segments[1]
                seq = segments[2]
                copy_num = next(
                    (int(s[5:]) for s in segments[3:] if s.startswith('DP:f:')), None)
                if ignore_n:
                    seq = seq.replace('n', '')
                if ignore_zero_copy and copy_num == 0:
                    continue
                if seq:
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


def parse_fa(filename):
    graph = nx.DiGraph()
    seqs = dict()
    fasta = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    for name, record in fasta.items():
        graph.add_node((name, '+'))
        graph.add_node((name, '-'))
        seqs[name] = len(record.seq)
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


def hsl(h: float):
    assert 0 <= h <= 1
    return 'hsl({}, 100%, 50%)'.format(h * 360)


def spring_layout(seqpositions, seqs, seqnames, match, graph, min_identity):
    """
    seqpositions[seqname] = (hapname, start, strand)
    """
    x0 = [seqpositions[seqname][1] for seqname in seqnames]

    def fun(x):
        y = 0
        for i, seqname in enumerate(seqnames):
            for m in match[seqname]:
                identity = m.mlen / m.blen
                if identity >= min_identity:
                    y += (x[i] - m.tstart) ** 2

            hapname, start, strand = seqpositions[seqname]
            for edge in graph.out_edges((seqname, strand)):
                _, (seqname_child, strand_child) = edge
                j = next((i for i, s in enumerate(
                    seqnames) if s == seqname_child), None)
                if j is not None:
                    y += 0.01 * ((x[i] + seqs[seqname]) - x[j]) ** 2
        return y

    eprint('minimizing', x0)
    x = minimize(fun, x0).x
    eprint('x=', x)

    seqpositions_new = dict()
    for i, seqname in enumerate(seqnames):
        (hapname, start, strand) = seqpositions[seqname]
        seqpositions_new[seqname] = (hapname, x[i], strand)

    return seqpositions_new


def main():
    parser = argparse.ArgumentParser(description='generate SVG')
    parser.add_argument('gfa_or_fa', type=Path, help='GFA or FA of asm')
    parser.add_argument('paf', type=Path, help='PAF of genome vs asm')
    parser.add_argument('--width', type=int,
                        default=1000, help='svg width')
    parser.add_argument('--margin', type=int, default=30, help='')
    parser.add_argument('--box_height', type=int, default=20, help='')
    parser.add_argument('--min_length', type=int, default=0, help='')
    parser.add_argument('--layout', choices=['mapq', 'first'], default='mapq',
                        help=('how to determine (primary) unitig position?'
                              'mapq: alignment of highest mapq will be used'
                              'first: first mapping in PAF file will be used'
                              '(useful when with manually reordered PAF file'))
    parser.add_argument('--order', type=str, nargs='+')
    parser.add_argument('--no_spring_layout', action='store_true', help='')
    parser.add_argument('--shade_by_identity',
                        action='store_true', help='')
    parser.add_argument('--min_identity', type=float, default=0.999, help='')
    parser.add_argument('--hide_text', action='store_true', help='')
    parser.add_argument('--draw_mismatch_primary_only',
                        action='store_true', help='')
    parser.add_argument('--draw_mismatch_threshold', type=int, default=None,
                        help='if the number of mismatches is below this threshold, visualize mismatch position in red line. To disable mismatch visualization, set threshold to zero.')
    args = parser.parse_args()

    if args.gfa_or_fa.suffix == '.gfa':
        graph, seqs = parse_gfa(args.gfa_or_fa)
    elif args.gfa_or_fa.suffix == '.fa':
        graph, seqs = parse_fa(args.gfa_or_fa)
    else:
        raise Exception("args.gfa_or_fa is neither .gfa nor .fa")

    match, haps, seqs_from_paf = parse_paf(args.paf)
    n_haps = len(haps)
    assert n_haps == 2, "genome is not diploid"

    for seqname, seqlen in seqs_from_paf.items():
        if seqname in seqs and seqs[seqname] != seqlen:
            print('seq "{}" length not match gfa={} paf={} using paf'.format(
                seqname, seqs[seqname], seqlen), file=sys.stderr)
            seqs[seqname] = seqlen

    # debug output
    hapnames = sorted(haps.keys())
    for i, hapname in enumerate(hapnames):
        print('HAP', i, hapname, haps[hapname], file=sys.stderr)

    def hit(record):
        return (record.tname, record.tstart, str(record.strand))

    def bestmatch(records):
        """convert List[paf records] into best hit (hap, pos, strand) (in terms of match length)"""
        if len(records) == 0:
            return (hapnames[0], 0, '+')

        if args.layout == 'first':
            record = records[0]
        elif args.layout == 'mapq':
            record = max(records, key=lambda record: record.mlen)
        else:
            raise ValueError('invalid layout mode')
        return hit(record)

    seqpositions = {seqname: bestmatch(match[seqname])
                    for seqname in seqs.keys()}
    seqnames = sorted([seqname for seqname in seqs.keys()],
                      key=lambda seqname: seqpositions[seqname])
    for i, seqname in enumerate(seqnames):
        print('SEQ', i, seqname, seqpositions[seqname], file=sys.stderr)

    if not args.no_spring_layout:
        # adjust horizontal layout by spring layout of alignment and graph connection
        seqpositions = spring_layout(
            seqpositions, seqs, seqnames, match, graph, args.min_identity
        )
        # sort by starting positions to determine vertical ordering
        seqnames = sorted([seqname for seqname in seqs.keys()],
                          key=lambda seqname: seqpositions[seqname])

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
    text_margin = 3
    margin = args.margin - 15 if args.hide_text else args.margin
    box_height = args.box_height
    bp_per_px = max(haps.values()) / (width - 2 * margin)
    font_size = 20
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

    elements_top = []
    elements = []

    # (1) haplotypes
    for hapname in hapnames:
        x_hap = margin
        y_hap = hap_to_y(hapname)
        if not args.hide_text:
            elements_top.append(text(x_hap, y_hap - text_margin, font_size,
                                     "{} {}bp".format(hapname, haps[hapname])))
        elements_top.append(rect(x_hap, y_hap,
                                 x(haps[hapname]), box_height,
                                 fill="#555", stroke="#555", stroke_width=2))

    # (2) sequences (contigs/unitigs)
    for i, seqname in enumerate(seqnames):
        color = COLORS[i % len(COLORS)]
        hapname, start, strand = seqpositions[seqname]
        y_seq = seq_to_y(seqname)
        x_seq_left, x_seq_right = seq_to_x(seqname)
        x_seq_width = x_seq_right - x_seq_left
        if not args.hide_text:
            elements_top.append(text(x_seq_left, y_seq - text_margin,
                                     font_size, "{}{} {}bp".format(seqname, strand, seqs[seqname])))
        elements_top.append(rect(x_seq_left, y_seq, x_seq_width, box_height,
                                 fill=color, stroke="#555", stroke_width=2))

        # (2a) show alignment
        for m in match[seqname]:
            hapname = m.tname
            mapq = m.mapq

            y_hap = hap_to_y(hapname)
            y_seq = seq_to_y(seqname)
            if hapname == hapnames[0]:
                y_hap += box_height
            else:
                y_seq += box_height

            # alignment start
            x_seq_start = x_seq_left + x(m.qstart)
            x_hap_start = margin + x(m.tstart)
            # alignment end
            x_seq_end = x_seq_left + x(m.qend)
            x_hap_end = margin + x(m.tend)

            # shade the best match
            if args.shade_by_identity:
                identity = m.mlen / m.blen
                # opacity(1.00) = 1
                # opacity(0.99) = 0
                min_identity = args.min_identity

                if identity >= min_identity:
                    ratio = max(0, (identity - min_identity) /
                                (1 - min_identity)) if min_identity < 1 else 1
                    # opacity = ratio * 0.3
                    opacity = 0.3
                    # eprint(m.qname, m.tname, opacity)
                    # lines
                    elements.append(
                        line(x_seq_start, y_seq, x_hap_start, y_hap, opacity=opacity))
                    elements.append(
                        line(x_seq_end, y_seq, x_hap_end, y_hap, opacity=opacity))
                    # shade
                    elements.append(poly(x_seq_start, y_seq, x_seq_end, y_seq,
                                         x_hap_end, y_hap, x_hap_start, y_hap,
                                         color=color, opacity=opacity))
                    # mismatch sign
                    cigar = m.get_tag("cs").value
                    for (tindex, qindex) in aligned_pairs(cigar, ignore_match=True):
                        x_hap = margin + x(m.tstart + tindex)
                        x_seq = x_seq_left + x(m.qstart + qindex)
                        # mismatch
                        elements.append(line(x_seq, y_seq, x_hap, y_hap,
                                             color=color, opacity=0.5))
            else:
                # lines
                elements.append(line(x_seq_start, y_seq, x_hap_start, y_hap))
                elements.append(line(x_seq_end, y_seq, x_hap_end, y_hap))

                if bestmatch(match[seqname]) == hit(m):
                    elements.append(poly(x_seq_start, y_seq, x_seq_end, y_seq,
                                         x_hap_end, y_hap, x_hap_start, y_hap,
                                         color=color, opacity=0.2))
                else:
                    elements.append(poly(x_seq_start, y_seq, x_seq_end, y_seq,
                                         x_hap_end, y_hap, x_hap_start, y_hap,
                                         color='#bbb', opacity=0.1))
                # mismatch
                nm = m.get_tag("NM").value
                cigar = m.get_tag("cs").value
                if args.draw_mismatch_threshold is None or nm <= args.draw_mismatch_threshold:
                    if args.draw_mismatch_primary_only and not m.is_primary():
                        continue
                    for (tindex, qindex) in aligned_pairs(cigar, ignore_match=True):
                        x_hap = margin + x(m.tstart + tindex)
                        x_seq = x_seq_left + x(m.qstart + qindex)
                        # mismatch
                        elements.append(line(x_seq, y_seq, x_hap, y_hap,
                                             color=color, opacity=0.5))

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
            elements_top.append(join(x_s, y_s, x_t, y_t,
                                     start_to_right, end_to_right, opacity=0.5))
    s = svg(elements + elements_top, width=width, height=height)
    print(s)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    main()
