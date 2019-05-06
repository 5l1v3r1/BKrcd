# BKrcd

This implementation is based on [D. Strash's work](https://github.com/darrenstrash/quick-cliques).

The algorithm is discussed in our paper *"Fast Maximal Clique Enumeration for Real-world Graphs"*.

## how to use

You could run `make` in algo\_imph directory first to build the program,

then the program is in `bin` directory and you can call it with parameters

For example:

```shell
$ cd algo_imph
$ make
$ cd bin
$ nqc --input-file=[your data file] --algorithm=newhybrid
```

There are three MCE approaches in our implementation. Everyone could try them by using `--algorithm=[newhybrid|naive|degeneracy]`

## Input Graph Format

The format of our codes is similar to edgelist format, but it has the total number of vertices in the first line.

For example:

```shell
nodenum
v0 v1
v1 v2
...
vi vj

```

For any questions, please contact liyinuo@hust.edu.cn or zyshao@hust.edu.cn .

