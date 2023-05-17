#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import argparse
from pathlib import Path
from dbgphmm import parse_inspect_file, CopyNums


def str_is_up(is_up):
    if is_up:
        return '+'
    else:
        return '-'


def copy_nums_to_history(copy_nums):
    ret = []
    x = copy_nums.copy_nums.copy()
    for info in reversed(copy_nums.infos):
        for update in info:
            if update.is_up:
                x[update.edge] -= 1
            else:
                x[update.edge] += 1
        update_x = ['{}{}'.format(x[u.edge], str_is_up(u.is_up)) for u in info]
        ret.append(update_x)
    ret.reverse()
    return ret


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('inspect', type=Path, help='INSPECT file')
    parser.add_argument('k', type=int, help='')
    parser.add_argument('depth', type=int, help='')
    args = parser.parse_args()
    inspects = parse_inspect_file(args.inspect, args.k)

    for copy_nums in inspects.copy_nums:
        if len(copy_nums.infos) == 0:
            print('init')
        elif len(copy_nums.infos) == args.depth:
            print('{}\t{}\t{}'.format(copy_nums.id, copy_nums.log_likelihood,
                  ', '.join([''.join(h) for h in copy_nums_to_history(copy_nums)])))


if __name__ == '__main__':
    main()
