#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict

from dbgphmm import parse_inspect_files_by_prefix, Inspect, p0, plt_common_settings


def draw_graph_props(inspects: List[Inspect], filename=None):
    """
    n_nodes and n_edges etc for all k.

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
    fks, fps = [], []
    tks, tps = [], []
    for inspect in inspects:
        # false copy_nums
        copy_nums = inspect.copy_nums
        fks += [inspect.k for c in copy_nums if c.dist_from_true != 0]
        fps += [c.posterior for c in copy_nums if c.dist_from_true != 0]
        # true copy_nums
        tks += [inspect.k for c in copy_nums if c.dist_from_true == 0]
        tps += [c.posterior for c in copy_nums if c.dist_from_true == 0]
    plt.scatter(fks, fps, c='blue', marker='o', alpha=0.5, label='X')
    plt.scatter(tks, tps, c='red', marker='*', alpha=0.5, label='X_true')
    plt.legend()

    # likelihood and prior
    plt.subplot(4, 1, 2)
    plt.ylabel('log P(R|X)')
    plt.xlabel('k')
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
    main()
