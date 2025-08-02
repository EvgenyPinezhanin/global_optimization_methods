#! /usr/bin/gnuplot

sampleName = "fitting_family_problems_operational_characteristics"

load "output_data/new_model_class/".sampleName."/vars.txt"

set linetype  1  lc rgb "red"         lw 2 dt 1
set linetype  2  lc rgb "green"       lw 2 dt 1
set linetype  3  lc rgb "blue"        lw 2 dt 1
set linetype  4  lc rgb "orange"      lw 2 dt 1
set linetype  5  lc rgb "brown"       lw 2 dt 1
set linetype  6  lc rgb "dark-yellow" lw 2 dt 1
set linetype  7  lc rgb "dark-violet" lw 2 dt 1
set linetype  8  lc rgb "cyan"        lw 2 dt 1

set linetype  9  lc rgb "blue"         lw 2 dt 2
set linetype 10  lc rgb "green"       lw 2 dt 2
set linetype 11  lc rgb "blue"        lw 2 dt 2
set linetype 12  lc rgb "orange"      lw 2 dt 2
set linetype 13  lc rgb "brown"       lw 2 dt 2
set linetype 14  lc rgb "dark-yellow" lw 2 dt 2
set linetype 15  lc rgb "dark-violet" lw 2 dt 2
set linetype 16  lc rgb "cyan"        lw 2 dt 2

set linetype 17  lc rgb "dark-violet"         lw 2 dt 4
set linetype 18  lc rgb "green"       lw 2 dt 4
set linetype 19  lc rgb "blue"        lw 2 dt 4
set linetype 20  lc rgb "orange"      lw 2 dt 4
set linetype 21  lc rgb "brown"       lw 2 dt 4
set linetype 22  lc rgb "dark-yellow" lw 2 dt 4
set linetype 23  lc rgb "dark-violet" lw 2 dt 4
set linetype 24  lc rgb "cyan"        lw 2 dt 4

set linetype cycle 24

fontName = "Helvetica, 20"

set grid

set xlabel "K" font fontName
set ylabel "P_s(K)" font fontName offset -1

set tics font fontName

set key inside right bottom
set key font fontName spacing 1.3


title(familyName) = sprintf("Operational characteristics on the %s", familyName)
titlePng(familyName) = ARG1 == 1 ? title(familyName) : sprintf("")

dataFile(methodNumber, variantNumber) = methodNumber == 1 ? \
    sprintf("output_data/new_model_class/%s/%s/%s_%s_%s_%s", sampleName, methodNames[methodNumber], familyName, density, \
                                                             key[variantNumber], r[variantNumber]) : \
                                        methodNumber == 2 ? \
    sprintf("output_data/new_model_class/%s/%s/%s", sampleName, methodNames[methodNumber], familyName) : \
    sprintf("output_data/new_model_class/%s/%s/%s_%s", sampleName, methodNames[methodNumber], familyName, population[variantNumber])

titleGraph(methodNumber, variantNumber) = methodNumber == 1 ? \
    sprintf("%s, den = %s, r = %s, key = %s", methodNames[methodNumber], density, r[variantNumber], key[variantNumber]) : \
                                          methodNumber == 2 ? sprintf("%s", methodNames[methodNumber]) : \
                                                              sprintf("%s, pop = %s", methodNames[methodNumber], population[variantNumber])

if (ARG1 == 0) {
    set lmargin 10

    set title title(familyName) font fontName
    plot for [i = 1 : numberMethods] for [j = 1 : numberVariants[i]] dataFile(i, j) using 1:2 with lines title titleGraph(i, j)

    bind all "alt-End" "exit gnuplot"
    pause mouse close
} else {
    set terminal pngcairo size 1440, 700 font fontName
    system "mkdir -p output_graph/new_model_class/".sampleName

    set lmargin 10
    set rmargin 12
    set tmargin 3
    set bmargin 3

    set output "output_graph/new_model_class/".sampleName."/".familyName.".png"
    
    set title titlePng(familyName) font fontName
    plot for [i = 1 : numberMethods] for [j = 1 : numberVariants[i]] dataFile(i, j) using 1:2 with lines lt (i - 1) * 8 + j title titleGraph(i, j)
}
