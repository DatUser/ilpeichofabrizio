# Physically Based Rendering on Sphere

Implementation of physically based rendering lightning that depends on
roughness and metalness of material.


## How to build and run

```
$ conan install .. -s build_type=Release --build missing
$ cmake -B build
$ cd build && make -j
$ ./main
```
