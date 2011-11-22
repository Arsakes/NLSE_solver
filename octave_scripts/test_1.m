function symulacja(s,t )

#%skrypt octave odpalający program
half_size= 1024;
xstep = 0.1;
#%stałe fizyczne

h_bar=0.658211915;
g_c=0.006;
g_r=2.0*g_c;
lam_c=0.33;
lam_r=1.5*lam_c;
R=0.00000160022;
m=0.040685;
#%szerokość próbki
d=5.0;
a=1.0/(d*2.0*pi);

#%przeskalowanie ze wzgledu na jednowymiarowość
R=a*R;
g_c=a*g_c;
g_r=a*g_r;
#%stałe minimalne
P0=lam_c*lam_r/R;
n0=lam_c/R;

#%dziedzina
X=linspace(-half_size*xstep, (half_size-1)*xstep, half_size * 2);
window_size = (half_size * 0.8)
#%przepisanie stałych
physical_constants=[h_bar, g_c, g_r, lam_c, lam_r, m, R, xstep, 0.01];
disp("Wczytanio stałe, fizyczne i dane symulacji - physical_constants")

#%potencjał
V=complex(zeros(1,half_size*2.0));

#%początkowe natężenie w rezerwuarze 
n_r=s*n0 * exp( -0.5 * (X/5.0).^2 ) ;
#%moc pompy
P_l=s*P0 * exp( -0.5 * (X/5.0).^2 ) ;
#%kork czasowy
time_steps=linspace(0.0, t , 128);
#%startowa funkcja
rnd=stdnormal_rnd(half_size*2,1) * 0.05;
PSI=exp(-(X/window_size).^8 ) .* (sqrt(complex(rnd, rnd)))';

%wypluwanie wyników
[PSI_wynik, NR_wynik]=NLS_solver(PSI,V, n_r, P_l,time_steps, physical_constants );
[xx,yy]=meshgrid(time_steps, X);

#%pierwszy wykres
#subplot(2,1,1)
#pcolor(xx,yy,abs(PSI_wynik) ), shading("interp");
#axis([0, max(time_steps) ,-(half_size)*xstep*0.5, (half_size -1)*xstep*0.5]);
#colormap(hot)
#title('wave function evolution')
#print 'wave_function_evolution' -dpng;
#%plot(X, abs(PSI_wynik));
#%axis([-half_size*xstep, (half_size -1)*xstep]);
#subplot(2,1,2)
#pcolor(xx,yy,NR_wynik), shading("interp");
#colormap(hot)
#axis([0, max(time_steps) ,-(half_size)*xstep*0.5, (half_size -1)*xstep*0.5]);
#title('reservoir polariton density');
#print 'reservoir' -dpng;

#%size(PSI_wynik)

plot(X, abs(PSI_wynik(:,128)) );
title('wave function t=10 ps');
#print 'wave_function_1d' -dpng
endfunction
