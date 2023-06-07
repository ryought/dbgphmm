from typing import List, Tuple, Dict
from dbgphmm import parse_bmap_file, MapInfoList
import matplotlib.pyplot as plt
import numpy as np


def forward_diff_ranking(bmap: List[List[MapInfoList]], bmap_true: List[List[MapInfoList]]):
    return [
        i
        for _, i in
        sorted(
            ((abs(x[-1].p_forward - y[-1].p_forward), i)
             for i, (x, y) in enumerate(zip(bmap, bmap_true))),
            reverse=True
        )
    ]


def backward_diff_ranking(bmap: List[List[MapInfoList]], bmap_true: List[List[MapInfoList]]):
    return [
        i
        for _, i in
        sorted(
            ((abs(x[0].p_backward - y[0].p_backward), i)
             for i, (x, y) in enumerate(zip(bmap, bmap_true))),
            reverse=True
        )
    ]


def draw_overviews(bmap: List[List[MapInfoList]], bmap_true: List[List[MapInfoList]], filename=None):
    """
    """
    plt.figure(figsize=(10, 10))

    dfs = []
    dbs = []
    dfbs = []
    dfbts = []
    ns = []
    nts = []
    for i, (x, y) in enumerate(zip(bmap, bmap_true)):
        df = abs(x[-1].p_forward - y[-1].p_forward)
        db = abs(x[0].p_backward - y[0].p_backward)
        dfb = abs(x[-1].p_forward - x[0].p_backward)
        dfbt = abs(y[-1].p_forward - y[0].p_backward)
        n = sum(w.n_active_nodes for w in x)
        nt = sum(w.n_active_nodes for w in y)
        dfs.append(df)
        dbs.append(db)
        dfbs.append(dfb)
        dfbts.append(dfbt)
        ns.append(n)
        nts.append(nt)

    n = len(bmap)
    plt.subplot(3, 2, 1)
    plt.bar(np.arange(n), dfs)
    plt.title('diff forward')

    plt.subplot(3, 2, 2)
    plt.bar(np.arange(n), dbs)
    plt.title('diff backward')

    plt.subplot(3, 2, 3)
    plt.bar(np.arange(n), dfbs)
    plt.title('forward vs backward in approx')

    plt.subplot(3, 2, 4)
    plt.bar(np.arange(n), dfbts)
    plt.title('forward vs backward in true')

    plt.subplot(3, 2, 5)
    plt.bar(np.arange(n) - 0.2, ns, label='approx')
    plt.bar(np.arange(n) + 0.2, nts, label='true')
    plt.legend()
    plt.title('# filled cells')

    if filename:
        plt.savefig(filename)
        plt.clf()
        plt.close()
    else:
        plt.show()


def draw_forward_backward_diff(i: int, bmap: List[List[MapInfoList]], bmap_true: List[List[MapInfoList]], filename=None, top=3):
    plt.figure(figsize=(10, 10))

    plt.subplot(5, 1, 1)
    pf = [mi.p_forward for mi in bmap[i]]
    pft = [mi.p_forward for mi in bmap_true[i]]
    plt.plot(pf, label='approx')
    plt.plot(pft, label='true')
    plt.legend()
    plt.xlabel('base')
    plt.ylabel('p_forward')

    plt.subplot(5, 1, 2)
    pf = [mi.p_backward for mi in bmap[i]]
    pft = [mi.p_backward for mi in bmap_true[i]]
    plt.plot(pf, label='approx')
    plt.plot(pft, label='true')
    plt.legend()
    plt.xlabel('base')
    plt.ylabel('p_backward')

    plt.subplot(5, 1, 3)
    plt.ylabel('active node id in true')
    for k in range(top):
        ss = [mi.infos_merged[k].edge_compact for mi in bmap_true[i]]
        plt.plot(ss, label=f'#{k}')
    plt.legend()

    plt.subplot(5, 1, 4)
    ns = [mi.n_active_nodes for mi in bmap[i]]
    nts = [mi.n_active_nodes for mi in bmap_true[i]]
    plt.plot(ns, label='approx')
    plt.plot(nts, label='true')
    plt.legend()
    plt.xlabel('base')
    plt.ylabel('n_active_nodes')

    plt.subplot(5, 1, 5)
    plt.ylabel('active node id in approx')
    for k in range(top):
        ss = [mi.infos_merged[k].edge_compact for mi in bmap[i]]
        plt.plot(ss, label=f'#{k}')
    plt.legend()

    plt.suptitle(f'read #{i}')

    if filename:
        plt.savefig(filename)
        plt.clf()
        plt.close()
    else:
        plt.show()


def draw_active_nodes_diff(i: int, bmap: List[List[MapInfoList]], bmap_true: List[List[MapInfoList]], filename=None, top=3, jaccard_top=10, ranks=[1, 5, 10, 20]):
    n_bases = len(bmap[i])
    plt.figure(figsize=(10, 10))

    plt.subplot(6, 1, 1)
    plt.xlim(0, n_bases)
    plt.title(f'read #{i}')
    pf = [mi.p_forward for mi in bmap[i]]
    pft = [mi.p_forward for mi in bmap_true[i]]
    plt.plot(pf, label='approx')
    plt.plot(pft, label='true')
    plt.xlabel('base')
    plt.ylabel('p_forward')

    plt.subplot(6, 1, 2)
    plt.xlim(0, n_bases)
    for k in range(top):
        ss = [mi.infos_merged[k].edge_compact for mi in bmap_true[i]]
        plt.plot(ss, label=f'#{k}')
    plt.legend()
    plt.ylabel('node id with top-k p_forward')

    plt.subplot(6, 1, 3)
    plt.xlim(0, n_bases)
    ns = [mi.n_active_nodes for mi in bmap[i]]
    nts = [mi.n_active_nodes for mi in bmap_true[i]]
    plt.plot(ns, label='approx')
    plt.plot(nts, label='true')
    plt.legend()
    plt.xlabel('base')
    plt.ylabel('n_active_nodes')

    plt.subplot(6, 1, 4)
    plt.xlim(0, n_bases)
    plt.ylim(0 - 0.1, 1 + 0.1)
    xs = []
    for j, (mi, mit) in enumerate(zip(bmap[i], bmap_true[i])):
        v = set(m.edge_compact for m in mi.infos_forward)
        w = set(m.edge_compact for m in mit.infos_forward)
        x = len(v & w) / len(v | w)
        xs.append(x)
    plt.ylabel('jaccard')
    plt.plot(xs)

    plt.subplot(6, 1, 5)
    plt.xlim(0, n_bases)
    plt.ylim(0 - 0.1, 1 + 0.1)
    xs = []
    for j, (mi, mit) in enumerate(zip(bmap[i], bmap_true[i])):
        vs = set(m.edge_compact for m in mi.infos_forward)
        ws = set(m.edge_compact for m in mit.infos_merged[:jaccard_top])
        x = len(vs & ws) / len(ws)
        xs.append(x)
    plt.ylabel('% active covered by forward')
    plt.plot(xs)

    plt.subplot(6, 1, 6)
    plt.xlim(0, n_bases)
    for rank in ranks:
        xs = [mi.infos_forward[0].prob -
              mi.infos_forward[rank].prob for mi in bmap_true[i]]
        plt.plot(xs, label=f'pf[0] - pf[{rank}]')
    plt.legend()

    if filename:
        plt.savefig(filename)
        plt.clf()
        plt.close()
    else:
        plt.show()
