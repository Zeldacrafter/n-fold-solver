#!/usr/bin/bash
NAME=$(echo "$1" | sed -e "s/.*\///g" -e "s/\..*//g")
./a.out < $1 > $NAME.lp
glpsol --max --cpxlp $NAME.lp -o $NAME.out
