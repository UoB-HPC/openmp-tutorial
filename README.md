# Programming Your GPU with OpenMP

This is a hands-on tutorial that introduces the basics of targetting GPUs with OpenMP 4.5 through a series of worked examples.

Starting with serial code, the tutorial takes you thorugh parallellising, exploring the performance characteristics, and optimising the following small programs:

* `vadd` – A simple vector addition program, often considered the "hello world" of GPU programming.
* `pi` – A numerical integration program that calculates and approximate value of π.
* `jac_solv` – A Jacobi solver.
* `heat` - An explicit finite difference 5-point stencil code.

## Usage

To build all the examples:

```bash
make
```

To run, submit jobs using your training account:

```bash
qsub submit_vadd     # For vector add
qsub submit_pi       # For pi
qsub submit_jac_solv # For Jacobi
qsub submit_heat     # For heat
```

## Publication history

This tutorial was presented [at SC'17](https://sc17.supercomputing.org/presentation/?id=tut127&sess=sess217).
A version of this tutorial was presented [at UK OpenMP Users' Conference 2018](https://www.eventbrite.co.uk/e/uk-openmp-users-conference-2018-tickets-37685633745).

