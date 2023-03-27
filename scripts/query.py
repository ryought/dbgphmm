from dbgphmm import *

inspects = parse_inspect_files_by_prefix('n/p01_u500_n4_pi_v2')

for inspect in inspects:
    if inspect.k == 50:
        print(inspect.k)
        x0 = inspect.copy_nums[0].copy_nums
        for c in inspect.copy_nums:
            print(c.infos)
            if len(c.infos) == 0:
                xi = c.copy_nums

        x = inspect.copy_nums[0].copy_nums.copy()
        ups = [23]
        downs = [1, 55, 58, 21, 37]

        for i in (ups + downs):
            print('x[{}]={}'.format(i, x[i]))

        for i in ups:
            x[i] += 1
        for i in downs:
            x[i] -= 1

        for i in (ups + downs):
            print('x[{}]={}'.format(i, x[i]))

        xa = xi.copy()
        for i in ups:
            xa[i] += 1
        for i in downs:
            xa[i] -= 1

        print('x =', x)
        print('x0=', x0)
        print('xi=', xi)
        print('xa=', xa)

        for c in inspect.copy_nums:
            if c.copy_nums == x0:
                print('x0 found', c.id)
            if c.copy_nums == x:
                print('x found', c.id)
            if c.copy_nums == xi:
                print('xi found', c.id)
            if c.copy_nums == xa:
                print('xa found', c.id)
