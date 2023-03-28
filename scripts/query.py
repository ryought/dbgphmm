#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from dbgphmm import *
import argparse


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('prefix', type=Path, help='prefix of output data')
    parser.add_argument('k', type=int)
    parser.add_argument('x', type=str, choices=['x0', 'xT'])
    parser.add_argument('--update', type=str)
    args = parser.parse_args()

    k = args.k
    inspect = parse_inspect_file(
        str(args.prefix) + '.k{}.inspect'.format(k), k)

    print(inspect.k)

    if args.x == 'x0':
        for c in inspect.copy_nums:
            if len(c.infos) == 0:
                c0 = c
    elif args.x == 'xT':
        c0 = inspect.copy_nums[0]
    else:
        raise Exception()

    x0 = c0.copy_nums
    print('c0=', c0)
    print('x0=', x0)

    x = x0.copy()

    for update in Update.parse_list(args.update):
        if update.is_up:
            x[update.edge] += 1
        else:
            x[update.edge] -= 1

    print('x =', x)

    for c in inspect.copy_nums:
        if c.copy_nums == x:
            print('c =', c)
            break
    else:
        print('c not found')


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    main()
