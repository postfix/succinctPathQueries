from math import log
import sys

if __name__ == '__main__':
    n = len(sys.stdin.readline())/2
    q = dict()
    i = 0
    while i < n:
        numbers = map(int, raw_input().split())
        for x in numbers:
            if q.has_key(x):
                v = q.get(x)
                assert v != None
                q[x] = v+1
            else:
                q[x] = 1
        i += len(numbers)
    e = 0
    for key,val in q.items():
        p = val/(n+0.00)
        e = e+p*log(p)
    print -e/log(2.00)
