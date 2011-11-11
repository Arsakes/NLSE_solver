#%skrypt octave odpalający program
half_size= 1024;
xstep = 0.05;
#%stałe fizyczne

h_bar=0.658211915;
g_c=0.006;
g_r=2.0*g_c;
lam_c=0.33;
lam_r=1.5*lam_c;
R=0.00000160022;
m=0.040685;
#%szerokość próbki
d=10.0;
a=1.0/(d*2.0*pi);

#%przepisanie stałych
physical_constants=[h_bar, a*g_c, a*g_r, lam_c, lam_r, m, a*R, xstep, 0.005];
#%potencjał

V=complex(zeros(1,half_size*2.0));
X=linspace(- half_size*xstep, (half_size-1)*xstep, half_size * 2);
#% tutaj PSI jest guassem PSI=complex(exp(-X.^2));
#%PSI = soliton(X,0.0,0.0);
n_r=real(zeros(1,half_size*2));
P_l=real(exp(-X.^2));
time_steps=linspace(0.0,0.1 , 16);
rnd=stdnormal_rnd(half_size*2,1) * 0.1;
PSI=complex(rnd, rnd);

%wypluwanie wyników
[PSI_wynik, NR_wynik]=NLS_solver(PSI,V, n_r, P_l,time_steps, physical_constants );
[xx,yy]=meshgrid(time_steps, X);
#%pierwszy wykres
subplot(2,1,1)
pcolor(xx,yy,abs(PSI_wynik) ), shading("interp");
axis([0, max(time_steps) ,-(half_size)*xstep*0.5, (half_size -1)*xstep*0.5]);
colormap(hot)
title('wave function evolution')
#%plot(X, abs(PSI_wynik));
#%axis([-half_size*xstep, (half_size -1)*xstep]);
subplot(2,1,2)
pcolor(xx,yy,NR_wynik);
axis([0, max(time_steps) ,-(half_size)*xstep*0.5, (half_size -1)*xstep*0.5]);
title('reservoir polariton density')

size(PSI_wynik)

