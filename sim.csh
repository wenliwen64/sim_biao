#!/bin/csh

root -b <<EOF >& log/$1.log
.L simulation.C
main($1)
.q
EOF

