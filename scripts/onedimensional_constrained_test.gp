#!/usr/bin/gnuplot

load "scripts/functions/onedimensional_constrained.txt"

fontName = "Helvetica, 16"

set samples 600

set grid

set xlabel "x" font "Helvetica, 16"
set ylabel "f(x)" font "Helvetica, 16" offset -1

set tics font fontName

set key box width -3 outside right top
set key font fontName spacing 1.3

title(block, n, method) = sprintf("Graph of the %s function №%d, %s method", block, n, method)
titlePng(block, n, method) = ARG2 == 1 ? title(block, n, method) : sprintf("")
functionName(i) = (i == 0 ? sprintf("f(x)") : sprintf("g%d(x)", i))

if (ARG2 == 0) {
    trialsFile = "output_data/".ARG1."_test/".blockNames[ARG3 + 1]."_".(ARG4 + 1)
    exist = system("[ ! -f ".trialsFile." ]; echo $?")
    if (exist == 1) {
        set lmargin 11

        set xrange [A(ARG3 + 1, ARG4 + 1) : B(ARG3 + 1, ARG4 + 1)]

        set title title(blockNames[ARG3 + 1], ARG4 + 1, ARG1) font fontName

        plot for [ind = 0 : numberConstraints(ARG3 + 1, ARG4 + 1)] f(ARG3 + 1, ARG4 + 1, x, ind) title functionName(ind), \
             trialsFile index 2 using 1:2 ls 4 lc rgb "green" lw 2 title "trial points", \
             trialsFile index 0 using 1:2 ls 7 lc rgb "red"   ps 2 title "X_{min}", \
             trialsFile index 1 using 1:2 ls 7 lc rgb "blue"  ps 1 title "X"

        bind all "alt-End" "exit gnuplot"
        pause mouse close
    } else {
        print("TrialsFile ".trialsFile." not exist")
    }
} else {
    set terminal pngcairo size 1640, 950 font fontName
    system "mkdir -p output_graph/".ARG1."_test"

    set lmargin 12
    set rmargin 18
    set tmargin 3
    set bmargin 3

    do for [i = 1 : numberBlocks] {
        do for [j = 1 : numberFunctions[i]] {
            trialsFile = "output_data/".ARG1."_test/".blockNames[i]."_".j
            exist = system("[ ! -f ".trialsFile." ]; echo $?")
            if (exist == 1) {
                set output "output_graph/".ARG1."_test/".blockNames[i]."_".j.".png"

                set xrange [A(i, j) : B(i, j)]

                set title titlePng(blockNames[i], j, ARG1) font fontName

                plot for [ind = 0 : numberConstraints(i, j)] f(i, j, x, ind) title functionName(ind), \
                     trialsFile index 2 ls 4 lc rgb "green" lw 2 title "trial points", \
                     trialsFile index 0 ls 7 lc rgb "red"   ps 2 title "X_{min}", \
                     trialsFile index 1 ls 7 lc rgb "blue"  ps 1 title "X"
            } else {
                print("TrialsFile ".trialsFile." not exist")
            }
        }
    }
}