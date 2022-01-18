# Pathtracer

Mederic Carriat
Gomez Xavier

# Results

https://youtu.be/QbSvGH9mtfk


## How to build and run

```
$ mkdir build; cd build
$ conan install .. -s build_type=Release --build missing
$ cmake .. -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DCMAKE_BUILD_TYPE=Release
$ make -j
$ ./main
```
