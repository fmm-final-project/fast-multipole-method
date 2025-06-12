# Naive Gravity

This folder contains code that reads the generated particles from a binary file and computes the gravitational forces acting on each of them using the naive $O(N^2)$ algorithm.

The computed forces are output in both CSV and binary formats, with each particle's force represented as $[F_x, F_y, F_z]$.

To run the code, execute the following commands:
```
make clean
make
./naive_gravity
```