#! /usr/bin/env python3

N = 2
while N <= 256:

    with open("Move_{}.dot".format(N), "w") as f:

        z = []
        for i in range(N // 2):
            z.append(i)
            z.append(i + N // 2)

        print("digraph G {", file = f)

        for i in range(N):
            print("    {} -> {}".format(i, z[i]), file = f)

        print("}", file = f)

    N *= 2
