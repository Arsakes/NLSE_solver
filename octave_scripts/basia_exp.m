# Symulacja eksperymentu Basi Piętki z FUW
# Dane próbki wzięte z artykułu Probing the Dynamics of Spontaneous Quantum Vortices
# in Polariton Superfluids
#
#

#skrypt octave odpalający program
half_size= 256;
xstep = 0.06;
time_step=0.0005;
s=4.0;             #stosunek mocy minimalnej do dostarczanej
#stałe fizyczne
# jednostki energia [J]
# przestrzeń mikrometr
# czas [ps]

#do przeliczanie jednostek
eV=1.602e-19;    #[jeden milielektronovolt = tyle dżuli J] 

#masa elektronu
m_e=9.10938e-19         #[J*(ps/μm)^2] #ok
h_bar=1.05457168e-22    #[J*ps]a  #OK
g=2.4e-6 * eV  /h_bar;        #[J*μm^2]  #OK
g_r=0.2e-6 *eV /h_bar ;       #OK
lam_c=0.65e-3*eV/h_bar       #[1/ps]
lam_a=10e-3 *eV /h_bar 	#[1/ps]
lam_i=1.3e-6*eV /h_bar
tau=1.0/(0.1e-3*eV/h_bar) 
#P_th=100e-18 / eV;    	# [J/ps] #Dane z eksperymentu basi
R=0.2e-3 * eV /h_bar ;  #OK

m=m_e * 1e-4         #OK
#szerokość próbki
d=1.8;
a=1.0/sqrt( (d.^2)*2.0*pi);
#przeskalowanie stałych ze względu na jednowymiarowaść problemu
R=a*R
g=g*a
g_r=g_r*a
#przepisanie stałych
physical_constants=[h_bar, g, g_r, lam_c, lam_a, lam_i, tau, m, R, xstep, time_step];
#potencjał

#progowa wartość mocy lasera, progowe n_a
P_0=(lam_c*lam_a)/R
n_0=lam_c/R

#---------------PRZESTRZEN--------------------------------------------------
X=linspace(- half_size*xstep, (half_size-1)*xstep, half_size * 2);


#-------------------POTENCJAŁ-----------------------------------------------
# periodyczny z okresem wziętym z obrazka
#
#

lattice = 3                  #[μm] okres przestrzenny
V0 = 1500.0*g;
V=complex( V0 * cos(2*pi* X / lattice) );
# tutaj PSI jest guassem PSI=complex(exp(-X.^2));
#PSI = soliton(X,0.0,0.0);
n_0= s*n_0*real(exp( -0.5*(X./5.0).^2));
P_l= s*P_0*real(exp( -0.5*(X./5.0).^2) ) ;
time_steps=linspace(0.0,100.0 , 256);
rnd=stdnormal_rnd(half_size*2,1) * 0.1;
#PSI=complex(exp(-(X.^2) );
PSI=sqrt(complex(rnd, rnd))';

#wypluwanie wyników
[PSI_wynik, NA_wynik, NI_wynik]=NLS_solver(PSI,V, n_0, n_0, P_l,time_steps, physical_constants );
[xx,yy]=meshgrid(time_steps, X);

#dla ułatwienia
PSI=PSI_wynik(:,255);
#już nie używane można przeskalować
V=V*max(abs(PSI).^2) / V0;
#pierwszy wykres
subplot(2,2,1)
    pcolor(xx,yy,abs(PSI_wynik).^2 ), shading("interp");
    #plot( X, real(PSI_wynik(:,512).^2 ),X, real(V)  ) ;
    axis([0, max(time_steps) ,-(half_size)*xstep, (half_size -1)*xstep]);
    title('wave function evolution |PSI|')
    colormap(jet)
#axis([-half_size*xstep, (half_size -1)*xstep]);
subplot(2,2,2)
#    pcolor(xx,yy, NA_wynik ), shading("interp");
#    axis([0, max(time_steps) ,-(half_size)*xstep*0.5, (half_size -1)*xstep*0.5]);
     plot(X, abs(PSI_wynik(:,255)).^2, X, real(PSI_wynik(:,255)),  X, V  );
     axis([-1,1]*half_size*xstep)
     title('PSI(x) norm -blue, real-blue, V(x) - red')
#    colormap(hot)
#subplot(2,2,4)
#     plot(X, NA_wynik(:,255), X, NI_wynik(:,255) )
#     axis([-1,1]*half_size*xstep)
#     title('reservoir polariton density n_a-blue, n_i-green ')
subplot(2,2,4)
      plot(X, arg(PSI))
      axis([-1,1]*half_size*xstep)
#zapisywanie wyników
#nr=NR_wynik(:, 512);
#psi=(abs(PSI_wynik(:, 512))).^2;

#
#printf("Zapisywanie do pliku \n ");
#save("-ascii", "X", "X");
#save("-ascii", "psi", "ps:w
#save("-ascii", "nr", "nr");
#save("-ascii", "P_l", "P_l");

#-------------------------DEBUG-----------------------------

phase = complex( (g_r*( NA_wynik(:,255) + NI_wynik(:,255) ) ));
phase -= i*0.5*(lam_c - R*NA_wynik(:,255));
phase +=  g* abs(PSI).^2;
phase *= time_step;
subplot(2,2,3)
     plot(X, real(phase), X, imag(phase) );
     axis([-1,1]*xstep*half_size)
     title('phase change in next step PSI(x)')
     disp(1.0 - lam_c/R /max( NA_wynik(:,255))  )
     disp('Maximal phase change from kinetic part of equation:')
     kinetic_part=( max( real(PSI))-min(real(PSI))) /(xstep .^2);
     kinetic_part*=h_bar *time_step/m
#size(PSI_wynik)

