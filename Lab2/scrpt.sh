set terminal post enhanced colour solid font 20  # wybor formatu, w jakim chcemy utworzyc wynikowy rysunek

set output "v_x.eps" # nazwa pliku wynikowego
set title "V(x)" # tytul wykresu
set xlabel "x" # etykieta osi OX
set ylabel "V(x)" # etykieta osi OY
set grid # wlaczenie widoczności siatki pomocniczej

plot "out.dat" pt 7 ps 0.5 t "V(x)" \
, x/16 + 0.125 t "x/16 + 1/8" w l lt rgb "#ff0033" \
, x/16 - 0.125 t "x/16 - 1/8" w l lt rgb "#ff0033" \
, -x*x/2 - 7*x/16 t "-x*x/2 - 7*x/16" w l lt rgb "#254FA1" \
, x*x/2 - 7*x/16 t "x*x/2 - 7*x/16" w l lt rgb "#254FA1"
# plot - polecenie rysowania pliku o podanej nazwie "out.dat"
# w p == with points
# t "dt = 0.1" == title "dt = 0.1"
