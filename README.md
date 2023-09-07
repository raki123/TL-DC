# TL;DC: "Too Long; Didn't Count"
A length limited path counter.

## Building
A simple `make` in the main directory should suffice to build TL;DC and (most) of its dependencies.
Running `make` will provide a binary `a.out` in the root directory of this repository and a binary (the same) `tldc` in the `build` directory.

### Dependencies and custom number types
The only dependency that this does not provide is GMP. You need to install it by hand.
On Ubuntu this works via `sudo apt install libgmp-dev`.

Alternatively, if you think you do not need integers that big, you can change the `Edge_weight` type to something you consider suitable [here](src/graph.h#L39).
In this case, you will need to additionally provide a function 
```
size_t get_offset(std::vector<Edge_weight> const& result);
```
that converts `result[0]` to `size_t`. 
See [here](src/graph.cpp#L46-L50) for examples. 

Note that this is also a good idea if you know that the result of path counting will be small, since GMP is fast but still much slower than a simple `uint64_t`.

## Usage:

TL;DC accepts path counting instances in DIMACS format from [ICGCA](https://afsa.jp/icgca/):

Each line begins with a letter that defines the rest of the line. The legal lines are:

    c Comment:the remainder of the line is ignored.
    p Problem: must be of form p edge n m where n is the number of nodes (to be numbered 1..n) and m the number of edges.
    e Edge: must be of form e v1 v2 where v1 and v2 are the endpoints of the edge.
    l Length: must be of form l len where len is the maximum path length (the number of edges on the path).
    t Terminals: must be of form t v1 v2 where v1 and v2 are the path terminals. When this line is not given, count the total number of paths between every vertex pair.

Instances must be given to TL;DC via stdin.

### Command line options
```
tldc [-s .] [-t .] [-m .] [-h]
       -s         --strategy STRATEGY     use STRATEGY.
                                          Must be one of:
                                          * auto: choose automatically (default)
                                          * pathfbs: frontier-based search along a path decomposition
                                          * nautyfbs: frontier-based search along a path decomposition with modulo automorphisms
                                          * treefbs: frontier-based search along a tree decomposition
                                          * nautydfs: depth first search with modulo automorphisms
       -t         --threads NUM           use NUM threads.
       -m         --memory_limit LIMIT    set a hardlimit of LIMIT MB memory.
       -h         --help                  print this help.
```

### Example invocation
Typical calls to TL;DC might look something like this:
```
./tldc -m 10000 < instances/one_pair/000.col
```
calls tldc with a memory limit of 10000 MB, all available threads, and an automatically chosen strategy on the instance in file `instances/one_pair/000.col`.


```
./tldc -t 1 -s nautydfs < instances/one_pair/002.col
```
calls tldc without a memory limit, one thread, and strategy nautydfs on the instance in file `instances/one_pair/000.col`.

Especially the strategy option is worth varying.

## Description
For a description of the solver and the techniques that were used see the [short abstract](solver_description.pdf) or the [slides](TL-DC_slides.pdf) of the presentation at FIT 2023.
