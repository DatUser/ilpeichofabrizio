# Physically Based Rendering on Sphere

Implementation of physically based rendering lightning that depends on
roughness and metalness of material.


## How to build and run

```
$ mkdir build; cd build
$ conan install .. -s build_type=Release --build missing
$ cmake .. -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DCMAKE_BUILD_TYPE=Release
$ make -j
$ ./main
```
