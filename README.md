# Adyar

Approximate Common Substring

## Dependencies

* A modern, C++11 ready compiler such as `g++` version 4.7 or higher or `clang` version 3.2 or higher.
* The [cmake](www.cmake.org) build system (Version >= 2.8.11).
* A 64-bit operating system. Either Mac OS X or Linux are currently supported.
* [sdsl-lite](https://github.com/simongog/sdsl-lite) is included as a submodule.
* For increased performance of SDSL, the processor of the system should support fast bit operations available in `SSE4.2`

## Compiling the code

The program can be compiled in one of the following two ways

### Compiling without cmake
First, Install sdsl library at the default location (i.e., home directory). Then, compile the src/main.cpp as follows:

     g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib src/main.cpp -o adyar.x -lsdsl -ldivsufsort -ldivsufsort64

### Compiling with cmake

[sdsl-lite](https://github.com/simongog/sdsl-lite)  library is included as a submodule. Initialize the submodule as below, if they are not already initialized.

     git submodule init
     git submodule update

Next, Create a build directory, for e.g.,

     cd ..
     mkdir build
     cd build

Finally, build the executable as follows.

     cmake ..
     make

## Running the program

The program is run as follows:
`adayar.x <text_file_name> <value_of_k>`

Example:
`adyar.x text.txt 1`

Output:
Runtime, value of symmetric acs distance
