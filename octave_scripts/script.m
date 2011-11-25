#skrypt octave odpalający program
half_size= 1024;
xstep = 0.1;
#stałe fizyczne
# jednostki energia [J]
# przestrzeń mikrometr
# czas [ps]

#masa elektronu
m_e=9.10938e-19       #[J*(ps/μm)^2] #ok
h_bar=1.05457168e-22; #[J*ps]a  #OK

g_c=9.6131e-25;        #[J*μm^2]  #OK
g_r=2.0*g_c;

lam_c=0.33           #[1/ps]
lam_r=1.5*lam_c      #[1/ps]

P_th=0.025;         #[W] #OK 2.5e-14 [J/ps]
R=lam_c*lam_r/P_th;  #OK
m=m_e * 5e-5;
#szerokość próbki
d=5.0;
a=1.0/sqrt( (d.^2)*2.0*pi);
#przeskalowanie stałych ze względu na jednowymiarowaść problemu
R=a*R
g_c=g_c*a
g_r=g_r*a

#przepisanie stałych
physical_constants=[h_bar, g_c, g_r, lam_c, lam_r, m, R, xstep, 0.1];
#potencjał

#progowa wartość mocy lasera, progowe n_r
P_0=(lam_c*lam_r)/R;
n_0=lam_c/R;

V=complex(zeros(1,half_size*2.0));
X=linspace(- half_size*xstep, (half_size-1)*xstep, half_size * 2);
# tutaj PSI jest guassem PSI=complex(exp(-X.^2));
#PSI = soliton(X,0.0,0.0);
n_r= 3.0*n_0*real(exp( -0.5*(X./10.0).^2));
P_l= 3.0*P_0*real(exp( -0.5*(X./10.0).^2) ) ;
time_steps=linspace(0.0,500.0 , 512);
rnd=stdnormal_rnd(half_size*2,1) * 0.1;
#PSI=complex(exp(-(X.^2) );
PSI=sqrt(complex(rnd, rnd))';

#wypluwanie wyników
[PSI_wynik, NR_wynik]=NLS_solver(PSI,V, n_r, P_l,time_steps, physical_constants );
[xx,yy]=meshgrid(time_steps, X);
#pierwszy wykres
subplot(2,1,1)
    #pcolor(xx,yy,abs(PSI_wynik) ), shading("interp");
    plot( X, abs(PSI_wynik(:,512).^2 ), X, NR_wynik(:,512), X, P_l ) ;
    #axis([0, max(time_steps) ,-(half_size)*xstep*0.5, (half_size -1)*xstep*0.5]);
    title('wave function evolution |PSI|')
    colormap(hot)
#plot(X, abs(PSI_wynik));
#axis([-half_size*xstep, (half_size -1)*xstep]);
subplot(2,1,2)
    pcolor(xx,yy, NR_wynik ), shading("interp");
    axis([0, max(time_steps) ,-(half_size)*xstep*0.5, (half_size -1)*xstep*0.5]);
    title('reservoir polariton density')
    colormap(hot)

#zapisywanie wyników
nr=NR_wynik(:, 512);
psi=(abs(PSI_wynik(:, 512))).^2;

#
printf("Zapisywanie do pliku \n ");
save("-ascii", "X", "X");
save("-ascii", "psi", "psi");
save("-ascii", "nr", "nr");
save("-ascii", "P_l", "P_l");

size(PSI_wynik)

