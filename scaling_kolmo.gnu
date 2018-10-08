set term post eps enha color 24
set out 'scaling.eps'

set log y

set pointsize 1.5

tau=2./3.
cs2=1./3.
k2=(2.*pi/50.)**2
ni=cs2*(tau-0.5)
u0=0.1

set xlabel 't'
set ylabel 'u( y = 12, t )'
set format y '{10^{%T}}'

p './scaling_kolmo' u ($0*500):1 t 'Data' w p pt 4,\
   u0*exp(-ni*k2*x) t 'Theory' w l


