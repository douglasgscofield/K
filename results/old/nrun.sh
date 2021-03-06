#!/bin/sh

PREFIX="KN";

s0=1;
s1=0;

h0=0;
h1=0;

# SR="0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1";

SR="0.00 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16 0.18 0.20 0.22 0.24 0.26 0.28 0.30 0.32 0.34 0.36 0.38 0.40 0.42 0.44 0.46 0.48 0.50 0.52 0.54 0.56 0.58 0.60 0.62 0.64 0.66 0.68 0.70 0.72 0.74 0.76 0.78 0.80 0.82 0.84 0.86 0.88 0.90 0.92 0.94 0.96 0.98 1.00";



UU0="0.02 0.2 0.5 1 1.5 2";
U1=0;

for U0 in $UU0; do
echo "for U0=$U0 ...";

h0=0;
N="${PREFIX}_U=(${U0},${U1})_s=(${s0},${s1})_h=(${h0},${h1}).txt ...";
echo "creating file $N";
./K.exe -nested -tableheadingonly > $N;
for i in $SR ; do
    echo "U=(${U0},${U1}) s=(${s0},${s1}) h=(${h0},${h1}) S=$i ...";
    ./K.exe -nested -table -U0 $U0 -U1 $U1 -s0 $s0 -s1 $s1-h0 $h0 -h1 $h1 -S $i >> $N;
done
echo "done";
echo

h0=0.02;
N="${PREFIX}_U=(${U0},${U1})_s=(${s0},${s1})_h=(${h0},${h1}).txt ...";
echo "creating file $N";
./K.exe -nested -tableheadingonly > $N;
for i in $SR ; do
    echo "U=(${U0},${U1}) s=(${s0},${s1}) h=(${h0},${h1}) S=$i ...";
    ./K.exe -nested -table -U0 $U0 -U1 $U1 -s0 $s0 -s1 $s1-h0 $h0 -h1 $h1 -S $i >> $N;
done
echo "done";
echo

done
