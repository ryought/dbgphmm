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
import csv
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict
import argparse
from pathlib import Path
from dataclasses import dataclass
import re
import numpy as np
import glob
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
        >>> c = CopyNums.parse('20\tC\t5\t0.8\t-100\t-50\t100\t2\t[]\t[1,1,1,1,2]')
        >>> c == CopyNums(
        ...     k=20,
        ...     id=5,
        ...     posterior=0.8,
        ...     log_likelihood=-100.0,
        ...     log_prior=-50.0,
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
            for info in t[8].lstrip('[').rstrip(']').split(',') if len(info) != 0
        ]
        copy_nums = list(map(int, t[9].lstrip('[').rstrip(']').split(',')))
        return CopyNums(
            k=int(t[0]),
            id=int(t[2]),
            posterior=float(t[3]),
            log_likelihood=float(t[4]),
            log_prior=float(t[5]),
            genome_size=int(t[6]),
            dist_from_true=option_int(t[7]),
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


def plt_common_settings(plt, inspects, prob=False):
    ks = [inspect.k for inspect in inspects]
    plt.xscale('log')
    plt.xticks(ks, ks, rotation=90, fontsize=5)
    plt.grid(axis='both')
    if prob:
        plt.ylim(0-0.1, 1+0.1)


def draw_graph_props(inspects: List[Inspect], filename=None):
    """
    n_nodes and n_edges etc for all k

    * n_edges/n_nodes in compact/full
    * degree
    * average copy_nums
    """
    plt.figure(figsize=(10, 10))
    ks = [inspect.k for inspect in inspects]
    ns = [int(inspect.props['n_edges_full']) for inspect in inspects]
    ms = [int(inspect.props['n_edges_compact']) for inspect in inspects]
    plt.subplot(3, 1, 1)
    plt.plot(ks, ns, label='n_edges_full', marker='o')
    plt.ylabel('n_edges_full')
    plt_common_settings(plt, inspects)

    plt.subplot(3, 1, 2)
    plt.plot(ks, ms, label='n_edges_compact', marker='o')
    plt.ylabel('n_edges_compact')
    plt_common_settings(plt, inspects)

    plt.subplot(3, 1, 3)
    es = [int(inspect.props['n_emittable_edges']) for inspect in inspects]
    plt.plot(ks, es, marker='o')
    plt.ylabel('n_emittable_edges')
    plt_common_settings(plt, inspects)

    if filename:
        plt.savefig(filename)
        plt.clf()
        plt.close()
    else:
        plt.show()


def draw_posterior_of_copy_nums(inspects, filename=None):
    """
    * posterior or likelihood of sampled copy numbers
    * n_samples
    """

    # posterior
    plt.figure(figsize=(10, 10))
    plt.subplot(4, 1, 1)
    plt.ylabel('P(X|R)')
    plt_common_settings(plt, inspects, prob=True)
    for inspect in inspects:
        # false copy_nums
        copy_nums = inspect.copy_nums
        ks = [inspect.k for c in copy_nums if c.dist_from_true != 0]
        ps = [c.posterior for c in copy_nums if c.dist_from_true != 0]
        plt.scatter(ks, ps, c='blue', marker='o', alpha=0.5)

        # true copy_nums
        ks = [inspect.k for c in copy_nums if c.dist_from_true == 0]
        ps = [c.posterior for c in copy_nums if c.dist_from_true == 0]
        plt.scatter(ks, ps, c='red', marker='*', alpha=0.5)

    # likelihood and prior
    plt.subplot(4, 1, 2)
    plt.ylabel('log P(R|X)')
    plt_common_settings(plt, inspects)
    xs, ys, zs = [], [], []
    txs, tys = [], []
    for inspect in inspects:
        copy_nums = inspect.copy_nums
        xs += [inspect.k for c in copy_nums]
        ys += [c.log_likelihood for c in copy_nums]
        zs += [c.log_prior for c in copy_nums]
        # true copy_nums only
        txs += [inspect.k for c in copy_nums if c.dist_from_true == 0]
        tys += [c.log_likelihood for c in copy_nums if c.dist_from_true == 0]
    plt.scatter(xs, ys, c=zs, marker='o', alpha=0.5)
    plt.colorbar(orientation='horizontal', aspect=100, label='log P(X)')
    plt.scatter(txs, tys, c='red', marker='*', alpha=0.5)

    # genome sizes
    plt.subplot(4, 1, 3)
    plt.ylabel('genome size')
    plt_common_settings(plt, inspects)
    xs, ys, zs = [], [], []
    for inspect in inspects:
        copy_nums = inspect.copy_nums
        xs += [inspect.k for c in copy_nums]
        ys += [c.genome_size for c in copy_nums]
        zs += [c.log_likelihood for c in copy_nums]
    plt.scatter(xs, ys, c=zs, marker='o', alpha=0.5)
    plt.colorbar(orientation='horizontal', aspect=100, label='log P(R|X)')

    # samples
    plt.subplot(4, 1, 4)
    plt.ylabel('# samples')
    plt_common_settings(plt, inspects)
    xs = [inspect.k for inspect in inspects]
    ys = [len(inspect.copy_nums) for inspect in inspects]
    plt.scatter(xs, ys)

    plt.tight_layout()

    if filename:
        plt.savefig(filename)
        plt.clf()
        plt.close()
    else:
        plt.show()


def draw_posterior_of_edges(inspects, filename=None):
    """
    p_zero or p_true of each edge by its true copy num and k
    (k-mer classification result)
    """
    plt.figure(figsize=(10, 10))

    # P(X(e)=X_true(e))
    plt.subplot(4, 1, 1)
    plt.ylabel('P(X(e)=X_true(e))')
    plt_common_settings(plt, inspects, prob=True)
    xs, ys, zs = [], [], []
    for inspect in inspects:
        xs += [e.k for e in inspect.edges]
        ys += [e.p_true for e in inspect.edges]
        zs += [e.copy_num_true for e in inspect.edges]
    plt.scatter(xs, ys, c=zs, alpha=0.5)
    plt.colorbar(orientation='horizontal', aspect=100, label='X_true(e)')

    # P(X(e)=0) for e: X_true(e)=0
    plt.subplot(4, 1, 2)
    plt.ylabel('P(X(e)=0)')
    plt_common_settings(plt, inspects, prob=True)
    xs, ys, zs = [], [], []
    for inspect in inspects:
        xs += [e.k for e in inspect.edges if e.copy_num_true == 0]
        ys += [e.p_zero for e in inspect.edges if e.copy_num_true == 0]
    plt.scatter(xs, ys, alpha=0.5, marker='o', label='X_true(e)=0')

    # P(X(e)=0) for e: X_true(e)>0
    xs, ys, zs = [], [], []
    x1s, y1s, z1s = [], [], []
    for inspect in inspects:
        if any(e.copy_num_true != -1 for e in inspect.edges):
            xs += [e.k for e in inspect.edges if e.copy_num_true != 0]
            ys += [e.p_zero for e in inspect.edges if e.copy_num_true != 0]
        else:
            x1s += [e.k for e in inspect.edges if e.copy_num_true == -1]
            y1s += [e.p_zero for e in inspect.edges if e.copy_num_true == -1]
    plt.scatter(xs, ys, alpha=0.5, marker='x', label='X_true(e)>0')
    plt.scatter(x1s, y1s, alpha=0.5, marker='x', label='X_true(e)=?')
    plt.legend()

    #
    plt.subplot(4, 1, 3)
    plt.ylabel('# purged edges')
    plt_common_settings(plt, inspects)
    x1s, x2s, x3s, y1s, y2s, y3s = [], [], [], [], [], []
    for inspect in inspects:
        if any(e.copy_num_true != -1 for e in inspect.edges):
            # true positive
            x1s.append(inspect.k)
            y1s.append(
                len([e for e in inspect.edges if e.copy_num_true == 0 and e.p_zero > p0]))
            # false positive
            x2s.append(inspect.k)
            y2s.append(
                len([e for e in inspect.edges if e.copy_num_true > 0 and e.p_zero > p0]))
        else:
            # true copy num is unknown
            x3s.append(inspect.k)
            y3s.append(
                len([e for e in inspect.edges if e.copy_num_true == -1 and e.p_zero > p0]))
    plt.scatter(x1s, y1s, marker='o', label='X_true(e)=0')
    plt.scatter(x2s, y2s, marker='x', label='X_true(e)>0')
    plt.scatter(x3s, y3s, marker='x', label='X_true(e)=?')
    plt.legend()

    # E[X(e)] vs X_true(e)
    plt.subplot(4, 1, 4)
    plt.ylabel('E[X(e)] - X_true(e)')
    plt_common_settings(plt, inspects)
    xs, ys, zs = [], [], []
    for inspect in inspects:
        if any(e.copy_num_true != -1 for e in inspect.edges):
            xs += [e.k for e in inspect.edges if e.copy_num_true != -1]
            ys += [e.copy_num_posterior -
                   e.copy_num_true for e in inspect.edges if e.copy_num_true != -1]
    plt.scatter(xs, ys, alpha=0.5, label='')

    plt.tight_layout()

    if filename:
        plt.savefig(filename)
        plt.clf()
        plt.close()
    else:
        plt.show()


def inspect_false_positives(inspects):
    for inspect in inspects:
        for edge in inspect.edges:
            if edge.copy_num_true > 0 and edge.p_zero > p0:
                print('k={} edge={}'.format(inspect.k, edge))


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('prefix', type=Path, help='prefix of output data')
    parser.add_argument('--draw', action='store_true')
    parser.add_argument('--query', action='store_true')
    # parser.add_argument('--show', action='store_true')
    args = parser.parse_args()
    inspects = parse_inspect_files_by_prefix(args.prefix)
    if len(inspects) == 0:
        raise Exception('no inspect file found')

    if args.draw:
        draw_graph_props(
            inspects,
            filename=str(args.prefix) + '.graph_props.pdf',
        )
        draw_posterior_of_copy_nums(
            inspects,
            filename=str(args.prefix) + '.post_copy_nums.pdf',
        )
        draw_posterior_of_edges(
            inspects,
            filename=str(args.prefix) + '.post_edges.pdf',
        )
    else:
        inspect_false_positives(inspects)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    main()
