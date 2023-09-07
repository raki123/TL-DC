# TL;DC: "Too Long; Didn't Count"
A length limited path counter.

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
