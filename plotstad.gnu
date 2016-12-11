set terminal latex
set output "StadTrag.tex"
set format xy "$%g"
set xlabel "$x$"
set ylabel "$y$"
set parametric
a = 1
b = 2
set xrange [-0.5*(a+1):b+0.5*(a+1)]
set yrange [-0.5:a+0.5]
set arrow 1 from 0,0 to b,0 nohead lw 0.1
set arrow 2 from 0,a to b,a nohead lw 0.1
set size square
plot [-pi/2:pi/2] -0.5*a*cos(t),(0.5*a)+0.5*a*sin(t) lc rgb "black",\
b+0.5*a*cos(t),0.5*a+0.5*a*sin(t) lc rgb "black",\
"gnuStadium.txt" with lines
