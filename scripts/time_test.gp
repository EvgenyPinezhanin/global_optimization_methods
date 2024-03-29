#! /usr/bin/gnuplot

datafile="output_data/mggsa_time_test.txt"

set xlabel "number of trial points"
set ylabel "time, c"
set grid

set title "Time test of mggsa" font "Helvetica Bold, 20"

if (ARG1 == 0) {
     plot datafile using 1:2 ls 5 lc rgb "red" title "calc mu", \
          datafile using 1:3 ls 5 lc rgb "green" title "calc z*", \
          datafile using 1:4 ls 5 lc rgb "blue" title "calc R", \
          datafile using 1:5 ls 5 lc rgb "yellow" title "trial in point", \
          datafile using 1:6 ls 5 lc rgb "black" title "add point in trials", \
          datafile using 1:7 ls 5 lc rgb "violet" title "add point in I"
}
if (ARG1 == 1) {
     plot datafile using 1:2 ls 5 lc rgb "red" title "calc mu", \
          datafile using 1:3 ls 5 lc rgb "green" title "calc z*", \
          datafile using 1:4 ls 5 lc rgb "blue" title "calc R", \
          datafile using 1:5 ls 5 lc rgb "yellow" title "trial in point", \
          datafile using 1:6 ls 5 lc rgb "black" title "add point in trials", \
          datafile using 1:7 ls 5 lc rgb "violet" title "add point in I", \
          datafile using 1:8 ls 5 lc rgb "black" title "calc P"
}

bind all "alt-End" "exit gnuplot"
pause mouse close
