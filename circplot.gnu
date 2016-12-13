reset
set terminal epslatex
set output "CircTraj.tex"
set multiplot layout 2,2
set format xy "$%g$"
set xlabel "$x$"
set ylabel "$y$"
set size square 
set parametric
unset key
set xrange [-1.5:1.5]
set yrange [-1.5:1.5]
#
plot [0:2*pi] cos(t),sin(t) lc rgb "black",\
"gnuCircle1.txt" with lines lc rgb "red" lt 1
#
plot [0:2*pi] cos(t),sin(t) lc rgb "black",\
"gnuCircle2.txt" with lines lc rgb "red lt 1
#
plot [0:2*pi] cos(t),sin(t) lc rgb "black",\
"gnuCircle3.txt" with lines lc rgb "red" lt 1
#
plot [0:2*pi] cos(t),sin(t) lc rgb "black",\
"gnuCircle4.txt" with lines lc rgb "red" lt 1
#
unset multiplot
unset output

set terminal x11
set multiplot layout 2,2
set format xy "%g"
set xlabel "x"
set ylabel "y"
set size square	
set parametric
unset key
set xrange [-1.5:1.5]
set yrange [-1.5:1.5]
#
plot [0:2*pi] cos(t),sin(t) lc rgb "black",\
"gnuCircle1.txt" with lines lc rgb "red" lt 1
#
plot [0:2*pi] cos(t),sin(t) lc rgb "black",\
"gnuCircle2.txt" with lines lc rgb "red" lt 1
#
plot [0:2*pi] cos(t),sin(t) lc rgb "black",\
"gnuCircle3.txt" with lines lc rgb "red" lt 1
#
plot [0:2*pi] cos(t),sin(t) lc rgb "black",\
"gnuCircle4.txt" with lines lc rgb "red" lt 1
#
unset multiplot
