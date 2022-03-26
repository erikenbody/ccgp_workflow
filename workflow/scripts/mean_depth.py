#!/usr/bin/env python3
import sys

infile = sys.argv[1]

nums = []
with open(infile, "r") as f:
    for line in f:
        line = line.strip().split()[2]
        nums.append(int(line))

mean = sum(nums) // len(nums)

print(mean)
