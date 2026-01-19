#!/bin/bash

bs=(0.5 1 2 4 8 16)
fall=(0.001 0.01 0.1 1 10 100 1000)

for x in "${bs[@]}"; do
for f in "${fall[@]}"; do
  python absorbing_boundaries.py $f $x
done
done