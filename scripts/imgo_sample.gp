#!/usr/bin/gnuplot

taskName = "imgo_sample"
trialfile = "output_data/".taskName.".txt"
load "output_data/".taskName."_opt.txt"

f1(i, x) = (i == 0 ? sin(x) : 1/0)

f2(i, x) = (i == 0 ? -2.0 * x + 3.0 : \
            i == 1 ? sin(x) : 1/0)

f3(i, x) = (i == 0 ? 5.0 * x ** 2 + 3.0 * x - 1.0 : \
            i == 1 ? x ** 2 - 0.05 : \
            i == 2 ? -x + 0.1 : 1/0)

f(i, j, x) = (i == 1 ? f1(j, x) : \
              i == 2 ? f2(j, x) : \
              i == 3 ? f3(j, x) : 1/0)

set grid

set xlabel "x"
set ylabel "f(x)"

title(n) = sprintf("Graph of the sample function №%d, method imgo", n)
titlePng(n) = (ARG1 == 1) ? title(n) : sprintf("")
functionName(i) = (i == 0) ? sprintf("f(x)") : sprintf("g%d(x)", i)

if (ARG1 == 0) {
    ind = 3 * (ARG2 - 1)

    set title title(int(ARG2)) font "Helvetica Bold, 20"

    plot for [i = 0 : numberConstraints[int(ARG2)]] f(ARG2, i, x) title functionName(i), \
         trialfile index ind + 2 ls 4 lc rgb "green" lw 2 title "trials", \
         trialfile index ind     ls 7 lc rgb "red"   ps 2 title "X_{min}", \
         trialfile index ind + 1 ls 7 lc rgb "blue"  ps 1 title "X"

    bind all "alt-End" "exit gnuplot"
    pause mouse close
} else {
    set terminal pngcairo size 1440, 800 font "Helvetica Bold, 15"
    system "mkdir -p output_graph/".taskName

    do for [i = 1 : 3] {
        set output "output_graph/".taskName."/".taskName."_".i.".png"

        ind = 3 * (i - 1)

        set title titlePng(i) font "Helvetica Bold, 15"

        plot for [j = 0 : numberConstraints[i]] f(i, j, x) title functionName(j), \
             trialfile index ind + 2 ls 4 lc rgb "green" lw 2 title "trials", \
             trialfile index ind     ls 7 lc rgb "red"   lw 6 title "X_{min}", \
             trialfile index ind + 1 ls 7 lc rgb "blue"  lw 1 title "X"
    }
}
