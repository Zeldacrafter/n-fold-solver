#!/usr/bin/bash
NAME=$(echo "$1" | sed -e "s/.*\///g" -e "s/\..*//g")
if [ ! -f "$NAME.lp" ]; then
    ./a.out < $1 > $NAME.lp
fi
glpsol --max --cpxlp $NAME.lp -o $NAME.out
