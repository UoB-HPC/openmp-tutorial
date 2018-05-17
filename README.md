# Programming Your GPU with OpenMP

This is a hands-on tutorial that introduces the basics of targetting GPUs with OpenMP 4.5 through a series of worked examples.

Starting with serial code, the tutorial takes you thorugh parallellising, exploring the performance characteristics, and optimising the following small programs:

* `vadd` – A simple vector addition program, often considered the "hello world" of GPU programming.
* `pi` – A numerical integration program that calculates and approximate value of π.
* `jac_solv` – A Jacobi solver.

## Usage

To build all the examples:

```bash
make
```

To run, submit jobs using your training account:

```bash
qsub submit_vadd     # For vector add
qsub submit_pi       # For pi
qsub submit_jac_solv # For Jacobi
```

## Publication history

This tutorial was presented [at SC'17](https://sc17.supercomputing.org/presentation/?id=tut127&sess=sess217).

