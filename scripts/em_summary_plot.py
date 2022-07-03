#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import argparse
from dataclasses import dataclass
from enum import Enum
# https://github.com/r1chardj0n3s/parse
from parse import *
import numpy as np
import matplotlib.pyplot as plt

@dataclass(frozen=True)
class ExperimentParams:
    u: int
    n: int
    di: float
    dh: float
    c: int
    l: float
    e: float
    seed: int
    def from_str(s):
        parsed = parse('u{u}n{n}di{di}dh{dh}c{c}l{l}e{e}seed{seed}', s)
        return ExperimentParams(
            u=int(parsed['u']),
            n=int(parsed['n']),
            di=float(parsed['di']),
            dh=float(parsed['dh']),
            c=int(parsed['c']),
            l=float(parsed['l']),
            e=float(parsed['e']),
            seed=int(parsed['seed']),
        )

@dataclass(frozen=True)
class ClassificationResult:
    tp: int
    fp: int
    tn: int
    fn: int
    xx: int
    tp0: int
    fp0: int
    tn0: int
    fn0: int
    xx0: int
    def from_str(s):
        parsed = parse('TP={tp}/{tp0};FP={fp}/{fp0};TN={tn}/{tn0};FN={fn}/{fn0};XX={xx}/{xx0};', s)
        return ClassificationResult(
            tp=int(parsed['tp']),
            fp=int(parsed['fp']),
            tn=int(parsed['tn']),
            fn=int(parsed['fn']),
            tp0=int(parsed['tp0']),
            fp0=int(parsed['fp0']),
            tn0=int(parsed['tn0']),
            fn0=int(parsed['fn0']),
            xx=int(parsed['xx']),
            xx0=int(parsed['xx0']),
        )

class Algorithm(Enum):
    TRUE = 1
    RAW = 2
    V3 = 3
    def from_str(s):
        if s == 'TRUE':
            return Algorithm.TRUE
        elif s == 'RAW':
            return Algorithm.RAW
        elif s == 'V3':
            return Algorithm.V3
        else:
            raise Exception

@dataclass(frozen=True)
class ExpResult:
    params: ExperimentParams
    algorithm: str
    classification: ClassificationResult
    def from_str(s):
        rows = s.split('\t')
        params = ExperimentParams.from_str(rows[0])
        classification = ClassificationResult.from_str(rows[5])
        algorithm = Algorithm.from_str(rows[1])
        return ExpResult(
            params=params,
            algorithm=algorithm,
            classification=classification,
        )

def param_to_color_by_error(param):
    if param.e == 0.01:
        # 1%~3% error
        return 'red'
    elif param.e == 0.001:
        # 0.1%~0.3% error
        return 'blue'
    elif param.e == 0.02:
        # 2%~6% error
        return 'green'
    else:
        raise Exception

def param_to_color_by_lambda(param):
    if param.l == 0.1:
        return 'red'
    elif param.l == 0.01:
        return 'blue'
    elif param.l == 0.001:
        return 'green'
    else:
        raise Exception

def plot_fn_tn_by_unit_size(rs):
    for (u, marker) in [(20, '+'), (10, '.'), (100, 'x')]:
        rs_v3 = [r for r in rs if r.algorithm == Algorithm.V3 and r.params.u == u]
        tn = [r.classification.tn for r in rs_v3]
        fn = [r.classification.fn for r in rs_v3]
        c = [param_to_color_by_error(r.params) for r in rs_v3]
        plt.scatter(
            x=np.array(tn),
            y=np.array(fn),
            c=c,
            marker=marker,
        )
    plt.xlabel('TN')
    plt.ylabel('FN')
    plt.show()

def plot_fn_tn_by_lambda(rs):
    # for (e, marker) in [(0.01, '+'), (0.001, '.'), (0.02, 'x')]:
    # for (e, marker) in [(0.01, '+')]:
    for (e, marker) in [(0.02, 'x')]:
        rs_v3 = [r for r in rs if r.algorithm == Algorithm.V3 and r.params.e == e]
        tn = [r.classification.tn for r in rs_v3]
        fn = [r.classification.fn for r in rs_v3]
        c = [param_to_color_by_lambda(r.params) for r in rs_v3]
        plt.scatter(
            x=np.array(tn),
            y=np.array(fn),
            c=c,
            marker=marker,
            label='e={}'.format(e),
        )
    plt.xlabel('TN')
    plt.ylabel('FN')
    plt.legend()
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('summary_filename', type=str)
    args = parser.parse_args()

    rs = []
    with open(args.summary_filename) as f:
        for line in f:
            er = ExpResult.from_str(line)
            rs.append(er)

    plot_fn_tn_by_lambda(rs)

def test():
    ep = ExperimentParams.from_str('u20n20di0.01dh0.01c10l0.001e0.01seed4')
    print(ep)
    cr = ClassificationResult.from_str('TP=74/14;FP=0/0;TN=745/21;FN=0/0;XX=0/0;')
    print(cr)

if __name__ == '__main__':
    main()
