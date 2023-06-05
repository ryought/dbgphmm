#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parsers for dbgphmm output files

# Supported

* INSPECT

# ToDo

* DBG
* GFA
* POST
* PATHS
* MAP
"""
import csv
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict
import argparse
from pathlib import Path
from dataclasses import dataclass
import re
import numpy as np
import glob
import networkx as nx
csv.field_size_limit(1000000000)

p0 = 0.8


def split_at(s: str, delimiter: str) -> List[str]:
    """
    >>> split_at('e21-e52+e3+', 'e')
    ['e21-', 'e52+', 'e3+']
    >>> split_at('', 'e')
    []
    """
    xs = [x for x in range(len(s)) if s[x] == delimiter]
    return [
        s[xs[i]:(xs[i+1] if i+1 < len(xs) else len(s))]
        for i in range(len(xs))
    ]


def option_int(s: str) -> int:
    return -1 if s == '?' else int(s)


def option_float(s: str) -> float:
    return -1.0 if s == '?' else float(s)


@dataclass(frozen=True)
class Update:
    edge: int
    is_up: bool

    def parse(s):
        """
        >>> Update.parse("e21+")
        Update(edge=21, is_up=True)
        >>> Update.parse("e5-")
        Update(edge=5, is_up=False)
        """
        # direction
        if s[-1] == '+':
            is_up = True
        elif s[-1] == '-':
            is_up = False
        else:
            raise Exception
        # edge
        edge = int(s.lstrip('e').rstrip('+-'))
        return Update(
            edge=edge,
            is_up=is_up,
        )

    def parse_list(s):
        """
        >>> Update.parse_list("e21+e5+e7-")
        [Update(edge=21, is_up=True), Update(edge=5, is_up=True), Update(edge=7, is_up=False)]
        """
        if s:
            return [Update.parse(u) for u in split_at(s, 'e')]
        else:
            return []


@dataclass(frozen=True)
class MapInfo:
    edge_compact: int
    base: int
    edge_full: int
    type: str
    prob: float

    def parse(s):
        """
        >>> MapInfo.parse('62-0:M(n1734)-7.45907')
        MapInfo(edge_compact=62, base=0, edge_full=1734, type='M', prob=-7.45907)
        >>> MapInfo.parse('99-11:I(n1)0.0000')
        MapInfo(edge_compact=99, base=11, edge_full=1, type='I', prob=0.0)
        """
        m = re.match('(\d+)-(\d+):([MID])\(n(\d+)\)(.*)', s)
        return MapInfo(
            edge_compact=int(m[1]),
            base=int(m[2]),
            type=m[3],
            edge_full=int(m[4]),
            prob=float(m[5]),
        )


@dataclass(frozen=True)
class MapInfoList:
    p_forward: float
    p_backward: float
    infos_forward: List[MapInfo]
    infos_backward: List[MapInfo]
    infos_merged: List[MapInfo]


def parse_bmap_file(filename: Path) -> List[List[MapInfoList]]:
    """
    bmap[read][pos] = MapInfoList
    """
    ret = []
    current_read = None
    with open(filename) as f:
        for line in f:
            if line[0] == '#':
                continue
            s = line.rstrip().split('\t')
            if s:
                read = int(s[0])
                pos = int(s[1])
                type = s[2]
                p = float(s[3])
                infos = [MapInfo.parse(t) for t in s[4].split(',')]

                if current_read != read:
                    print(current_read)
                    ret.append([])
                    current_read = read
                    assert len(ret) == current_read + 1

                if type == 'F':
                    p_forward = p
                    infos_forward = infos
                elif type == 'B':
                    p_backward = p
                    infos_backward = infos
                elif type == 'S':
                    infos_merged = infos
                    ret[current_read].append(MapInfoList(
                        p_forward=p_forward,
                        p_backward=p_backward,
                        infos_forward=infos_forward,
                        infos_backward=infos_backward,
                        infos_merged=infos_merged,
                    ))
                    assert len(ret[current_read]) == pos

    return ret


@dataclass(frozen=True)
class CopyNums:
    """
    C section in INSPECT
    """
    k: int
    id: int
    posterior: float
    log_likelihood: float
    log_prior: float
    log_euler: float
    genome_size: int
    dist_from_true: int
    infos: List[List[Update]]
    copy_nums: List[int]

    def parse(s):
        """
        >>> c = CopyNums.parse('20\tC\t5\t0.8\t-100\t-50\t92.2\t100\t2\t[e21-e52+,e2+e5-]\t[1,1,1,1,2]')
        >>> c == CopyNums(
        ...     k=20,
        ...     id=5,
        ...     posterior=0.8,
        ...     log_likelihood=-100.0,
        ...     log_prior=-50.0,
        ...     log_euler=92.2,
        ...     genome_size=100,
        ...     dist_from_true=2,
        ...     infos=[
        ...         [Update(edge=21, is_up=False), Update(edge=52, is_up=True)],
        ...         [Update(edge=2, is_up=True), Update(edge=5, is_up=False)],
        ...     ],
        ...     copy_nums=[1,1,1,1,2],
        ... )
        True
        >>> c = CopyNums.parse('20\tC\t5\t0.8\t-100\t-50\t30\t100\t2\t[]\t[1,1,1,1,2]')
        >>> c == CopyNums(
        ...     k=20,
        ...     id=5,
        ...     posterior=0.8,
        ...     log_likelihood=-100.0,
        ...     log_prior=-50.0,
        ...     log_euler=30.0,
        ...     genome_size=100,
        ...     dist_from_true=2,
        ...     infos=[],
        ...     copy_nums=[1,1,1,1,2],
        ... )
        True
        """
        t = s.split()
        assert(t[1] == 'C')
        infos = [
            Update.parse_list(info)
            for info in t[9].lstrip('[').rstrip(']').split(',') if len(info) != 0
        ]
        copy_nums = list(map(int, t[10].lstrip('[').rstrip(']').split(',')))
        return CopyNums(
            k=int(t[0]),
            id=int(t[2]),
            posterior=float(t[3]),
            log_likelihood=float(t[4]),
            log_prior=float(t[5]),
            log_euler=float(t[6]),
            genome_size=int(t[7]),
            dist_from_true=option_int(t[8]),
            infos=infos,
            copy_nums=copy_nums,
        )


@dataclass(frozen=True)
class Edge:
    """
    E section in INSPECT
    """
    k: int
    id: int
    copy_num_true: int
    copy_num_posterior: float
    p_true: float
    p_zero: float
    distribution: str

    def parse(s):
        """
        >>> Edge.parse('20\tE\te50\t1\t2.00\t1.00\t0.00\tp(2)=1.0,p(3)=0.1')
        Edge(k=20, id=50, copy_num_true=1, copy_num_posterior=2.0, p_true=1.0, p_zero=0.0, distribution='p(2)=1.0,p(3)=0.1')
        """
        t = s.split()
        assert(t[1] == 'E')
        return Edge(
            k=int(t[0]),
            id=int(t[2].lstrip('e')),
            copy_num_true=option_int(t[3]),
            copy_num_posterior=float(t[4]),
            p_true=option_float(t[5]),
            p_zero=float(t[6]),
            distribution=t[7],
        )


@dataclass(frozen=True)
class Inspect:
    k: int
    copy_nums: List[CopyNums]
    edges: List[Edge]
    props: Dict[str, str]


def parse_inspect_file(filename: Path, k: int) -> Inspect:
    """
    Parse INSPECT tsv file
    """
    cs = []
    es = []
    gs = {}
    with open(filename) as f:
        for line in f:
            s = line.rstrip().split('\t')
            if s[0][0] == '#':
                continue
            if s:
                if s[1] == 'C':
                    cs.append(CopyNums.parse(line))
                elif s[1] == 'E':
                    es.append(Edge.parse(line))
                elif s[1] == 'G':
                    key = s[2]
                    value = s[3]
                    gs[key] = value
    return Inspect(
        k=k,
        copy_nums=cs,
        edges=es,
        props=gs,
    )


def parse_inspect_files_by_prefix(prefix: Path) -> List[Inspect]:
    ret = []
    for filename in glob.glob(str(prefix) + '*.inspect'):
        m = re.match('.*k(\d+).*', filename)
        if m:
            k = int(m.group(1))
            inspect = parse_inspect_file(filename, k)
            ret.append(inspect)
    ret.sort(key=lambda i: i.k)
    return ret


@dataclass(frozen=True)
class Emit:
    node: int
    prob: float


def parse_emits_list(s: str) -> List[Emit]:
    """
    >>> a = parse_emits_list('42:99.1,3:0.8,1:0.0,24:0.0')
    >>> a == [Emit(node=42, prob=99.1), Emit(node=3, prob=0.8), Emit(node=1, prob=0.0), Emit(node=24, prob=0.0)]
    True
    """
    return [Emit(node=int(t.split(':')[0]), prob=float(t.split(':')[1]))
            for t in s.split(',')]


def parse_map_file(filename: Path) -> List[List[List[Emit]]]:
    """
    """
    ret = []
    current_read = None
    currens_pos = None
    with open(filename) as f:
        for line in f:
            if line[0] == '#':
                continue
            s = line.rstrip().split('\t')
            if s:
                read = int(s[0])
                pos = int(s[1])
                base = s[2]
                emits = parse_emits_list(s[3])

                if current_read != read:
                    print(current_read)
                    ret.append([])
                ret[read].append(emits)
                current_read = read
                current_pos = pos
    return ret


def jaccard(a: List[Emit], b: List[Emit]) -> float:
    """
    >>> a = [Emit(node=42, prob=99.1), Emit(node=3, prob=0.8), Emit(node=1, prob=0.0)]
    >>> b = [Emit(node=42, prob=99.1), Emit(node=3, prob=0.8), Emit(node=0, prob=0.0)]
    >>> jaccard(a, b) == 2/4
    True
    """
    sa = set(e.node for e in a)
    sb = set(e.node for e in b)
    return len(sa & sb) / len(sa | sb)


def jaccards(map_a, map_b):
    ret = []
    n_reads = len(map_a)
    for i in range(n_reads):
        s = [jaccard(map_a[i][j], map_b[i][j]) for j in range(len(map_a[i]))]
        ret.append(s)
    return ret


def parse_dbg_file(filename: Path):
    """
    """
    g = nx.MultiDiGraph()
    with open(filename) as f:
        for line in f:
            if line[0] == '#':
                continue
            s = line.rstrip().split('\t')
            if s[0] == 'K':
                k = int(s[1])
            elif s[0] == 'N':
                node = int(s[1])
                kmer = s[2]
                g.add_node(node, id=node, kmer=kmer)
            elif s[0] == 'E':
                edge = int(s[1])
                source = int(s[2])
                target = int(s[3])
                kmer = s[4]
                copy_num = int(s[5])
                full_edges = [int(t) for t in s[6].split(',')]
                g.add_edge(source, target, id=edge, kmer=kmer,
                           copy_num=copy_num, full_edges=full_edges)
    return g


def plt_common_settings(plt, inspects, prob=False):
    ks = [inspect.k for inspect in inspects]
    plt.xscale('log')
    plt.xlabel('k')
    plt.xticks(ks, ks, rotation=90, fontsize=5)
    plt.grid(axis='both')
    if prob:
        plt.ylim(0-0.1, 1+0.1)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
