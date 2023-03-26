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
"""
from dataclasses import dataclass
from typing import List, Tuple, Dict
import argparse
from pathlib import Path
import re
import glob
import csv
csv.field_size_limit(1000000000)


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
    genome_size: int
    dist_from_true: int
    infos: List[List[Update]]
    copy_nums: List[int]

    def parse(s):
        """
        >>> c = CopyNums.parse('20\tC\t5\t0.8\t-100\t-50\t100\t2\t[e21-e52+,e2+e5-]\t[1,1,1,1,2]')
        >>> c == CopyNums(
        ...     k=20,
        ...     id=5,
        ...     posterior=0.8,
        ...     log_likelihood=-100.0,
        ...     log_prior=-50.0,
        ...     genome_size=100,
        ...     dist_from_true=2,
        ...     infos=[
        ...         [Update(edge=21, is_up=False), Update(edge=52, is_up=True)],
        ...         [Update(edge=2, is_up=True), Update(edge=5, is_up=False)],
        ...     ],
        ...     copy_nums=[1,1,1,1,2],
        ... )
        True
        """
        t = s.split()
        assert(t[1] == 'C')
        infos = [
            [Update.parse(u) for u in split_at(info, 'e')]
            for info in t[8].lstrip('[').rstrip(']').split(',')
        ]
        copy_nums = list(map(int, t[9].lstrip('[').rstrip(']').split(',')))
        return CopyNums(
            k=int(t[0]),
            id=int(t[2]),
            posterior=float(t[3]),
            log_likelihood=float(t[4]),
            log_prior=float(t[5]),
            genome_size=int(t[6]),
            dist_from_true=int(t[7]),
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
            copy_num_true=int(t[3]),
            copy_num_posterior=float(t[4]),
            p_true=float(t[5]),
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
        k = int(m.group(1))
        inspect = parse_inspect_file(filename, k)
        ret.append(inspect)
    ret.sort(key=lambda i: i.k)
    return ret


def draw_graph_props(inspects):
    """
    n_nodes and n_edges etc for all k
    """
    pass


def draw_posterior_of_copy_nums(inspects):
    """
    * posterior or likelihood of sampled copy numbers
    * n_samples
    """
    pass


def draw_posterior_of_edges(inspects):
    """
    p_zero or p_true of each edge by its true copy num and k
    (k-mer classification result)
    """
    pass


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('prefix', type=Path, help='prefix of output data')
    # parser.add_argument('--backward', action='store_true')
    # parser.add_argument('--show', action='store_true')
    args = parser.parse_args()
    inspects = parse_inspect_files_by_prefix(args.prefix)
    for inspect in inspects:
        print(inspect.k)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    main()
