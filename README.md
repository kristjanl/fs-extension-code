# CFlow

The tool CFlow is an extension of the contiuuous reachability part of the tool [Flow*](https://flowstar.org/), with the focus on solving systems in a compositional way.

## Installation
CFlow is implemented based on the same open-source software as Flow*.

* [GMP Library](http://gmplib.org/)
* [GNU MPFR Library](http://www.mpfr.org/)
* [GNU Scientific Library](http://www.gnu.org/software/gsl/)
* [GNU Linear Programming Kit](http://www.gnu.org/software/glpk/)
* [Bison - GNU parser generator](http://www.gnu.org/software/bison/)
* [flex: The Fast Lexical Analyzer](http://flex.sourceforge.net/)
* [Gnuplot](http://www.gnuplot.info/)

CFlow source can be compiled with the `make` tool.

When using user installation for dependencies the following lines need to be commented in `Makefile`
```
CFLAGS += -I $(DEP_HOME)/include
LINK_FLAGS += -L $(DEP_LIB) -Wl,-rpath,$(DEP_LIB)
```
and the `DEP_HOME` variable in `makefile.local` needs point to the appropriate directory. For example, if the dependencies are installed in home directory under `local/cflow`, then the variable needs to be
```
DEP_HOME = $(HOME)/local/cflowdep
```

## Usage

### Solving systems
CFlow can be used to integrate a system given by model file `reach.model` by using the command
```
./flowstar < reach.model
```

CFlow uses largely the same syntax for models as does the Flow*.

CFlow supports continuous reachability models with `poly ode 1` solving strategy of the tool Flow*. The differences between CFlow and Flow* model files are limited to options in setting scope

Whether to use native Flow* solving strategy or the compositional one is determined by the option
```
use cflow
```

When using compositional strategy model file also needs to specify the components. This can be done automatically (using as many components as small components as possible) with the option `auto components`, not using any component with `no components` or by specifying the components manally (for example `[[x1,x2],[y1,y2]]` specifies two components, one with variables `x1` and `x2` and one with variable `y1` and `y2`).

Cflow supports the following techniques to process the flowpipe between integration steps:
Name | example of the setting option
---|---
identity precondition | `identity precondition`
QR preconditioning | `QR precondition`
parallelepiped preconditioning | `parallelepiped precondition`
shrink wrapping after every `n`steps | `shrink wrapping 10`
shrink wrapping based on the size of the remainder | `shrink wrapping rem`
no processing | `no preocessing`

With CFlow, identity preconditioning and QR preconditioning can be either left model compositional or fully compositional (preceding the option with `left model compositional` or `fully compositional`, respectively)

An example of the settings part of the CFlow model file is the following:
```
setting
{
  fixed steps 1
  time 30
  remainder estimation 1e-2
  identity precondition
  gnuplot interval x1,x1
  fixed orders 5
  cutoff 1e-15
  precision 53
  output mult_comp_my_id
  use cflow
  auto components
  print on
}
```

### Results

Running CFlow on a model file produces a data file in csv directory containing the interval boxes of the flowpipes for each time step.

### Comparing methods

### As a library