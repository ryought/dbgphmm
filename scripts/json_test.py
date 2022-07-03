#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import json

def main():
    j = json.loads('[{"hoge": 10},{"hoge":11}]')
    print(j)

if __name__ == '__main__':
    main()
