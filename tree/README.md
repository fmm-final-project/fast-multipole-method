# Tree Algorithm

This folder contains code that reads the generated particles from a binary file and computes the gravitational forces acting on each of them using the $O(N \log N)$ Tree Algorithm.

The computed forces are output in CSV format, with each particle's force represented as $[F_x, F_y, F_z]$.

To run the code, execute the following commands:
```
make clean
make
./tree
```