set terminal post enhanced colour solid font 20  # wybor formatu, w jakim chcemy utworzyc wynikowy rysunek

set output "x_t.eps" # nazwa pliku wynikowego
set title "Wychylenie x(t)" # tytul wykresu
set xlabel "t" # etykieta osi OX
set ylabel "x(t)" # etykieta osi OY
set grid # wlaczenie widoczności siatki pomocniczej
set samples 200
plot "out.dat" pt 7 ps 0.5 t "x(t), dt = h = 0.1", cos(x) t "cos(t)" w l lt rgb "#254FA1"
# plot - polecenie rysowania pliku o podanej nazwie "out.dat"
# w p == with points
# t "dt = 0.1" == title "dt = 0.1"
