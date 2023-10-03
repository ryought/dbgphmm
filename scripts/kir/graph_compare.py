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


def rect(x, y, width, height, fill="red"):
    return '<rect x="{}" y="{}" width="{}" height="{}" fill-opacity="0.5" fill="{}" />'.format(x, y, width, height, fill)


def text(x, y, size, text):
    return '<text x="{}" y="{}" font-size="{}">{}</text>'.format(x, y, size, text)


def line(x1, y1, x2, y2):
    return '<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="black" stroke-width="1" />'.format(x1, y1, x2, y2)


def path(x0, y0, x1, y1, x2, y2, x3, y3):
    return '<path d="M {} {} C {} {}, {} {}, {} {}" stroke="green" fill="transparent" stroke-width="5"/>'.format(x0, y0, x1, y1, x2, y2, x3, y3)


def join(x1, y1, x2, y2, start_to_right, end_to_right, dx=100):
    x1c = x1 + dx if start_to_right else x1 - dx
    y1c = y1
    x2c = x1 + dx if end_to_right else x2 - dx
    y2c = y2
    return path(x1, y1, x1c, y1c, x2c, y2c, x2, y2)


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
                seqs[name] = seq
            elif segments[0] == 'L':
                name_a = segments[1]
                strand_a = segments[2]
                name_b = segments[3]
                strand_b = segments[4]
                graph.add_edge((name_a, strand_a), (name_b, strand_b))
    return graph, seqs


def parse_paf(filename):
    match = defaultdict(list)
    haps = dict()
    with PafFile(filename) as paf:
        for record in paf:
            # nm = record.get_tag("NM").value
            match[record.qname].append(record)
            haps[record.tname] = record.tlen
    return match, haps


def eprint(*args):
    print(*args, file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description='generate SVG')
    parser.add_argument('gfa', type=Path, help='GFA of asm')
    parser.add_argument('paf', type=Path, help='PAF of genome vs asm')
    parser.add_argument('--width', type=int, default=1000, help='svg width')
    parser.add_argument('--margin', type=int, default=30, help='')
    parser.add_argument('--box_height', type=int, default=20, help='')
    args = parser.parse_args()

    graph, seqs = parse_gfa(args.gfa)
    eprint(graph.nodes)
    n_seqs = len(seqs)

    match, haps = parse_paf(args.paf)
    eprint(match)
    eprint(haps)
    n_haps = len(haps)
    assert n_haps == 2, "genome is not diploid"

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

    width = args.width
    margin = args.margin
    box_height = args.box_height
    bp_per_px = max(haps.values()) / (width - 2 * margin)
    font_size = 20
    px_per_row = margin + box_height
    # raw position function
    def x(bp): return bp / bp_per_px
    def y(n): return margin + n * px_per_row
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
        length = len(seqs[seqname])
        left = margin + x(start)
        right = margin + x(start) + x(length)
        return (left, right)

    print('<svg width="{}" height="{}" xmlns="http://www.w3.org/2000/svg">'.format(width, height))

    for hapname in hapnames:
        x_hap = margin
        y_hap = hap_to_y(hapname)
        print(text(x_hap, y_hap, font_size, hapname))
        print(rect(x_hap, y_hap, x(haps[hapname]), box_height, fill="red"))

    for i, seqname in enumerate(seqnames):
        hapname, start, strand = seqpositions[seqname]
        y_seq = seq_to_y(seqname)
        x_seq_left, x_seq_right = seq_to_x(seqname)
        x_seq_width = x_seq_right - x_seq_left
        print(text(x_seq_left, y_seq, font_size, "{} {}".format(seqname, strand)))
        print(rect(x_seq_left, y_seq, x_seq_width, box_height, fill="blue"))

        # show alignment
        for m in match[seqname]:
            hapname = m.tname

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

        # show graph
        for edge in graph.out_edges((seqname, strand)):
            _, (seqname_child, strand_child) = edge
            print("EDGE", seqname, strand, seqname_child,
                  strand_child, file=sys.stderr)
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
            print(join(x_s, y_s, x_t, y_t, start_to_right, end_to_right))
    print('</svg>')


if __name__ == '__main__':
    main()
