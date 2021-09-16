#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
visualize training log

python3 plotter.py hoge.csv hoge.json
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

from parsers import *

def dist(state_a, state_b):
    return np.linalg.norm(np.array(state_a) - np.array(state_b))

def diff(state_a, state_b):
    return [b - a for (a, b) in zip(state_a, state_b)]

def diff_index(state_a, state_b):
    diffs = []
    for i, (a, b) in enumerate(zip(state_a, state_b)):
        if a != b:
            diffs.append((i, b-a))
    return diffs

def plot_temp(logs, stats, filename):
    x = np.array([log.id for log in logs])

    plt.subplot(4, 1, 1)
    plt.plot(x, np.array([log.temp for log in logs]))
    plt.ylabel('temp')

    plt.subplot(4, 1, 2)
    plt.plot(x, np.array([log.now_score for log in logs]))
    plt.plot(x, np.array([log.now_score_prior for log in logs]))
    plt.plot(x, np.array([log.now_score_forward for log in logs]))
    plt.ylabel('score')

    plt.subplot(4, 1, 3)
    plt.plot(x, np.array([log.now_size for log in logs]))
    plt.ylabel('size')

    plt.subplot(4, 1, 4)
    plt.plot(x, np.array([sum(log.now_state) for log in logs]))
    plt.plot(x, np.array([dist(log.now_state, logs[0].now_state) for log in logs]))
    plt.ylabel('dist')

    plt.xlabel('id')
    plt.savefig(filename)

def plot_projected(logs, stats, filename):
    pca = PCA(n_components=2)
    X = np.array([log.now_state for log in logs])
    C = np.array([log.now_score for log in logs])
    T = np.array([log.id for log in logs])
    pca.fit(X)
    Y = pca.transform(X)

    plt.subplot(2, 1, 1)
    plt.scatter(Y[:, 0], Y[:, 1], c=C)
    plt.scatter(Y[0, 0], Y[0, 1], marker='+', color='red')
    plt.colorbar()
    plt.xlabel('PC1')
    plt.ylabel('PC2')

    plt.subplot(2, 1, 2)
    plt.scatter(Y[:, 0], Y[:, 1], c=T)
    plt.scatter(Y[0, 0], Y[0, 1], marker='+', color='red')
    plt.colorbar()
    plt.xlabel('PC1')
    plt.ylabel('PC2')

    plt.savefig(filename)

def plot_grad_projected(logs, comments, stats, filename):
    # pca = PCA(n_components=2)
    tsne = TSNE(n_components=2, random_state = 0, perplexity = 30, n_iter = 400)
    X = np.array([log.next_state for log in logs])
    C = np.array([log.next_score for log in logs])
    # pca.fit(X)
    # Y = pca.transform(X)
    print('TSNE')
    Y = tsne.fit_transform(X)

    x = np.array([log.now_state for log in logs])
    c = np.array([log.now_score for log in logs])
    # y = pca.transform(x)

    plt.subplot(2, 1, 1)
    plt.scatter(Y[:, 0], Y[:, 1], c=C, alpha=0.2)
    # plt.scatter(y[:, 0], y[:, 1], marker='+', color='red')
    plt.colorbar()
    plt.xlabel('PC1')
    plt.ylabel('PC2')

    plt.savefig(filename)

def plot_grad_projected_with_true(logs, true_logs, stats, filename, method='pca'):
    X = np.array(
        [log.next_state for log in logs] + [log.next_state for log in true_logs]
    )
    C = np.array(
        [log.next_score for log in logs] + [log.next_score for log in true_logs]
    )
    N = len(logs)
    M = len(true_logs)

    # Y <- X
    if method == 'pca':
        pca = PCA(n_components=2)
        Y = pca.fit_transform(X)
    elif method == 'tsne':
        tsne = TSNE(n_components=2, random_state = 0, perplexity = 30, n_iter = 400)
        Y = tsne.fit_transform(X)
    else:
        return

    x = np.array([log.now_state for log in logs])
    c = np.array([log.now_score for log in logs])
    # y = pca.transform(x)

    plt.subplot(2, 1, 1)

    plt.scatter(Y[:, 0], Y[:, 1], c=C[:], alpha=0.2)

    # plt.scatter(Y[:N, 0], Y[:N, 1], c=C[:N], alpha=0.2)
    # plt.scatter(Y[N:N+M, 0], Y[N:N+M, 1], c=C[N:N+M], alpha=0.2, marker='+')

    # plt.scatter(y[:, 0], y[:, 1], marker='+', color='red')
    plt.colorbar()
    plt.xlabel('PC1')
    plt.ylabel('PC2')

    plt.savefig(filename)

def plot_grad_history(logs, true_logs, comments, true_comments, stats, filename):
    true_score = true_comments[0].score
    x = [comment.score - true_score for comment in comments]
    plt.plot(x)
    plt.xlabel('iteration')
    plt.ylabel('score - true_score')
    plt.tight_layout()
    plt.savefig(filename)

def plot_basis(logs, stats, filename):
    plt.figure(figsize=(10, 20))
    N = 6

    x = [cycle['id'] for cycle in stats['cycles']]
    lengths = [cycle['len'] for cycle in stats['cycles']]
    plt.subplot(N, 1, 1)
    plt.bar(x, lengths)

    # up/down: how many times the modification of i-th base occured and how it changed the score?
    up = [[] for _ in range(stats['cycle_summary']['n_cycles'])]
    down = [[] for _ in range(stats['cycle_summary']['n_cycles'])]
    for log in logs:
        d = diff_index(log.now_state, log.next_state)
        ds = log.next_score - log.now_score
        i, direction = d[0]
        if direction == 1:
            up[i].append(ds)
        else:
            down[i].append(ds)

    plt.subplot(N, 1, 2)
    count_up = [len(u) for u in up]
    count_down = [len(d) for d in down]
    plt.bar(x, count_up, label='up')
    plt.bar(x, count_down, bottom=count_up, label='down')
    plt.legend()

    plt.subplot(N, 1, 3)
    plt.scatter(x, [min(u) if len(u) > 0 else 0 for u in up])
    plt.scatter(x, [max(u) if len(u) > 0 else 0 for u in up])

    plt.subplot(N, 1, 4)
    plt.scatter(
        [lengths[diff_index(log.now_state, log.next_state)[0][0]] for log in logs],
        [log.next_score - log.now_score for log in logs]
    )

    plt.subplot(N, 1, 5)
    plt.scatter(lengths, count_up, alpha=0.2)

    plt.subplot(N, 1, 6)
    plt.scatter(lengths, count_down, alpha=0.2)

    plt.savefig(filename)

def plot_freq_em_history(logs, comments, stats, filename):
    ids = [log.id for log in logs]
    scores = [log.now_score for log in logs]
    true_score = [comment.score for comment in comments if comment.type == 'ONCE'][0]

    plt.plot(ids, scores, label='optimize')
    plt.plot(ids, np.ones(len(ids)) * true_score, label='true')
    plt.legend()
    plt.savefig(filename)

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('json_filename', type=str, help='json created by dbgphmm stats')
    parser.add_argument('tsv_filename', type=str, help='tsv created by dbgphmm optimize')
    parser.add_argument('--true_tsv_filename', type=str, help='true tsv created by dbgphmm optimize (from true)')
    parser.add_argument('--optimize_mode', type=str, choices=['annealer', 'grad', 'freq-em'], help='optimize method used by dbgphmm when producing tsv')
    args = parser.parse_args()

    # parse
    stats = parse_stats(args.json_filename)
    logs, comments = parse_logs(args.tsv_filename, with_comments=True)
    if args.true_tsv_filename:
        true_logs, true_comments = parse_logs(args.true_tsv_filename, with_comments=True)

    # plot
    if args.optimize_mode == 'annealer':
        plot_temp(logs, stats, args.tsv_filename + '.tempscore.png')
        plot_projected(logs, stats, args.tsv_filename + '.projected.png')
        plot_basis(logs, stats, args.tsv_filename + '.basis.png')
    elif args.optimize_mode == 'grad':
        if args.true_tsv_filename:
            # plot_grad_projected_with_true(logs, true_logs, stats, args.tsv_filename + '.projgrad.png', method='tsne')
            plot_grad_history(logs, true_logs, comments, true_comments, stats, args.tsv_filename + '.history.png')
        else:
            plot_grad_projected(logs, comments, stats, args.tsv_filename + '.projgrad.png')
    elif args.optimize_mode == 'freq-em':
        plot_freq_em_history(logs, comments, stats, args.tsv_filename + '.projfreqem.png')

if __name__ == '__main__':
    main()
