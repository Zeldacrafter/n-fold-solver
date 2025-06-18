# n-fold ILP solver

An implementation of a solver for bounded integer linear programs with $n$-fold block structure.  
The algorithm is based on the paper [Faster Algorithms for Integer Programs with Block Structure](https://arxiv.org/abs/1802.06289) by Eisenbrand , Hunkenschröder and Klein  

For $n, r, s, r \in \mathbb{Z}_{>0}$ an $n$-fold is a matrix $A \in \mathbb{Z}^{(r+ns)\times nt}$ with the block structure

$$
A = \begin{pmatrix}
    A_1 & A_2 & \dots & A_n \\
    B_1 & 0 & \dots & 0 \\
    0 & B_2 & \dots & 0\\
    \vdots & \vdots & \ddots & \vdots \\
    0 & 0 & \dots & B_n
\end{pmatrix},
$$

where $A_1, \dots, A_n \in \mathbb{Z}^{r\times t}$ and $B_1, \dots, B_n \in \mathbb{Z}^{s \times t}$.

For $c,l,u \in \mathbb{Z}^{nt}$ and $b \in \mathbb{Z}^{r+ns}$ we want to find a vector $x \in \mathbb{Z}^{nt}$ with $l \leq x \leq u$ 
such that $Ax = b$ and the objective function $c^T x$ is maximized. This corresponds to the integer linear program

$$
    \max \{ c^T x | Ax = b, l \leq x \leq u, x \in \mathbb{Z}^{nt}\}.
$$

## Usage
The size of the n-fold needs  to be specified at compile time. This can be done manually by setting the 
`PARAM_N`, `PARAM_R`, `PARAM_S` and `PARAM_T` variables in the makefile or by specifying an input file in the variable
`INPUT_FILE`. By default, the input file `input/input_1.in` is specified.

The project can then be build using CMake with
```shell
cmake -S . -B build-dir
cmake --build build-dir
```
Then execute the produced executable and provide the input via stdin:
```shell
./build-dir/n-fold-solver.out < input/input_1.in
```

## Input File Format
The input file format is
```
n r s t
l_1 ... l_nt
u_1 ... u_nt
b_1 ... b_(r+st)
c_1 ... c_nt
A1_(1,1) ... A1_(1,t) A1_(2,1) ... A1_(r,t)
...
AN_(1,1) ... AN_(1,t) AN_(2,1) ... AN_(r,t)
B1_(1,1) ... B1_(1,t) B1_(2,1) ... B1_(s,t)
...
BN_(1,1) ... BN_(1,t) BN_(2,1) ... BN_(s,t)
```
## Used libraries

* [eigen](https://eigen.tuxfamily.org/index.php)
* [boost](https://www.boost.org/)
* [hopscotch-map](https://github.com/Tessil/hopscotch-map)
