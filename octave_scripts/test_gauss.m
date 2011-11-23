#Function for comparing gaussian packet evoltion
# test_gauss(p, t)
# t >0 - time
# p - momentum
#
function [X,Ys, Ya, time]=test_gauss(p,t)

tic
# skrypt octave odpalający program
half_size= 512;
xstep = 0.5;
# stałe fizyczne

h_bar=1.0;
g_c=0;
g_r=0;
lam_c=0;
lam_r=0;
R=0;
m=1.0;

# ---ROZWIĄZYWANE RÓWNANIE
# dziedzina

X=linspace(-half_size*xstep, (half_size-1)*xstep, half_size * 2);
# przepisanie stałych

physical_constants=[h_bar, g_c, g_r, lam_c, lam_r, m, R, xstep, 0.01];
#disp("Wczytanio stałe, fizyczne i dane symulacji - physical_constants")

# potencjał
V=complex(zeros(1,half_size*2.0));

# początkowe natężenie w rezerwuarze  tu zero!
n_r= zeros(1,half_size*2);
# moc pompy
P_l= zeros(1,half_size*2) ;
# kork czasowy
time_steps=linspace(0.0, t , 2);
# startowa funkcja
Ys=  complex(gauss(X, p, 0.0, 1.0)) ;

#disp("Liczenie")
#wypluwanie wyników
[PSI_wynik, NR_wynik]=NLS_solver(Ys,V, n_r, P_l,time_steps, physical_constants );

Ys=abs( PSI_wynik(:,2));
#analityczne rozwiązanie
Ya= abs( gauss(X, p, t, 1.0) );

#axis([-half_size*xstep, half_size*xstep ])

#title('analitical vs simulated solution');

time=toc
endfunction
