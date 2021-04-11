# n-fold ILP solver

## Usage
Install `eigen` and `boost`. Under Arch Linux it would be
```shell
pacman -S boost eigen
```

Edit the first line in `compile.sh` to specify any of the libraries that are not in `$PATH`.  
For example if eigen is located in `/usr/include/eigen3` but not in `$PATH`:
```shell
local_opts="-I /usr/include/eigen3"
```

Compilation via the `compile.sh` script:
```shell
./compile.sh input.in                  # Compilation for input 'input.in'
./compile.sh input.in exec_name.out    # Produce 'exec_name.out' as executable
./compile.sh 1 2 3 4                   # Compilation for n = 1, r = 2, s = 3, t = 4
./compile.sh 1 2 3 4 exec_name.out     # Produce 'exec_name.out' as executable
```
Afterwards execute the produced executable and provide the input via stdin:
```shell
./exec_name.out < input.in
```

## Used libraries

* [eigen 3.3.9](https://eigen.tuxfamily.org/index.php)
* [boost 1.75.0](https://www.boost.org/)
* [hopscotch-map](https://github.com/Tessil/hopscotch-map)