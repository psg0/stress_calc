#!/bin/sh
l[1]="L111"
l[2]="L121"
l[3]="L221"
l[4]="L311"
l[5]="L321"
l[6]="L331"
l[7]="L113"
l[8]="L123"
l[9]="L313"
l[10]="L333"
for i in `seq 1 10`
do
fname6="./fig6/${l[i]}.dat";
gnuplot << EOF
set title "${l[i]} semi-infinite duration";
set xlabel "t";
plot "${fname6}" u 1:2 w lp t "ABC";
replot "${fname6}" u 1:3 w lp t "CDA";
replot "${fname6}" u 1:4 w lp t "ABC+CDA";
set terminal postscript eps color
set output "./fig6/${l[i]}.eps"
replot
quit
EOF
done
for i in `seq 1 10`
do
fname7="./fig7/${l[i]}.dat";
gnuplot << EOF
set title "${l[i]} finite duration";
set xlabel "t";
plot "${fname7}" u 1:2 w lp t "ABC";
replot "${fname7}" u 1:3 w lp t "CDA";
replot "${fname7}" u 1:4 w lp t "ABC+CDA";
set terminal postscript eps color
set output "./fig7/${l[i]}.eps"
replot
quit
EOF
done
