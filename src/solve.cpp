#include"solve.h"

const double PI = 3.1415926535;

solution::solution(const std::vector<complex>& temp, const std::vector<complex>& temp_V, 
     const std::vector<double>& temp_n, const std::vector<double>& temp_p,
     double xstep =0.01, double tstep=0.001, physical_constants C = physical_constants()   )
{

   step_spatial=xstep;
   step_temporal=tstep;
   //inicjacja stałych fizycznych
   CONST.h_bar = C.h_bar;
   CONST.g_r = C.g_r;
   CONST.g_c = C.g_c;
   CONST.lam_c = C.lam_c;
   CONST.lam_r = C.lam_r;
   CONST.m = C.m;
   CONST.R = C.R;

   data_size = temp.size();
   V = std::vector<complex>(temp_V);
   n_r = std::vector<double>(temp_n);
   P_l = std::vector<double>(temp_p);

   //inicjacja fftw
   FT=fourier_transform();
   psi = reinterpret_cast<complex*>( fftw_malloc(sizeof(complex) * data_size) );
   FT.init( psi, psi, data_size);
   
   for(int i = 0; i < data_size; ++i)
   {
       psi[i] = temp[i];
   }
   
}//spoko



//
// ROZWIĄZYWANE ZAGADNIENIE
//  równanie zacytowano dla ustelenia znaków w problemie
//  1)  dF_dt = i*dF_dxx*[h_bar/2m] - i*[g/h_bar]*|F|^2 *F - i[V/h]F
//
//  2)  dN_dt = Pl (x,t ) − γR n_R (x,t ) − R(n_R (x,t ))|ψ(x,t )|2 
//  
solution::~solution()
{
   fftw_free(psi);
   //fftw_free(A_n);
}

void inline solution::first_step()
{
   //jedna część równania z pochodną przestrzenną jest
   //ewoluowana do następnej chwili czasowej
   FT.compute( forward );
   for(int i =0; i<data_size ; ++i)
   {
      //normalizacja
      psi[i] = psi[i] / complex(data_size,0.0);
      //coś co wyskakuje przy liczeniu pochodnej
      //numer indeksu po którym częstotliwość zaczyna być dodatnia
      int i_max = (data_size-data_size%2)/2 - 1 + data_size%2;
      int k= i;
      if( i > i_max )
      { 
          k = i- data_size;
      }
      double k_square = std::pow( PI*k*2.0 / (data_size*step_spatial), 2.0);
       //std::cout<<std::sqrt(k_square)<<std::endl;
       psi[i] = psi[i] * std::exp(complex(0.0,-k_square ) * step_temporal * 0.5*CONST.h_bar/CONST.m ) ;   
    }


    FT.compute( backward );
}//spoko

void inline solution::second_step()
{ 
    complex phase(0.0,0.0);
    for(int i = 0 ; i < data_size;  ++i)
    { 
        //EWOLUCJA "psi"
        //ewolucja napędzane nieliniowością i potencjałem
         phase = complex(0.0, CONST.g_c * std::norm(psi[i]) ); 
         phase += complex(0.0,1.0) * V[i] ;
        
        //ewolucja napędzana przez rozpraszanie i zmiany w ilości polarytonów
         phase += complex( 0.5* (CONST.lam_c - R(n_r[i]))*CONST.h_bar, CONST.g_r *n_r[i] );  // plus działanie funkcji R na n_r
        //aplikacja do funkcji falowej
         psi[i] = psi[i] * std::exp( -phase * step_temporal /CONST.h_bar );

        //EWOLUCJA "n_r"
        n_r[i] += ( P_l[i] - CONST.lam_r * n_r[i] - R(n_r[i] ) * std::norm(psi[i]) )*step_temporal;

    }
}//spoko

void solution::make_time_step()
{
    //jak wskazywałaby nazwa krok drugi jest wykonywany pierwszy
    first_step();
    second_step();
}//spoko

void solution::evolution(int steps )
{
    for(int i = 0; i < steps; i++)
    {
            make_time_step();
            apply_boundary_condition(complex(0.0,0.0), complex(0.0,0.0), 0.1 );
    }
}//spoko

//
void solution::apply_boundary_condition(complex a=complex(0.0,0.0), complex b=complex(0.0,0.0), double c = 0.1)
{
   int falloff = std::floor( c * static_cast<double>(data_size) );
   int k;
   double x;
   //jak długim obszarze rociąga się "schodek" od 0 do 1
   
   //założenie 
   for( int i=0; i < falloff; ++i)
   { 
       //translacja po y, mnożenie, translacja powrotna
       //lewy warunek brzegowy 
       k=i;
       x = 1.8 * (static_cast<double>( i)/static_cast<double>(falloff-1) ) - 1.8;
       psi[k]=(psi[k] - a)*std::exp( -std::pow( x, 4)) + a;

       //prawy warunek brzegowy, exp zanika w lewo, idziemy w lewo
       k=data_size-1-i;
       x = -1.8* (static_cast<double>(-i)/static_cast<double>(falloff-1) ) + 1.8;
       psi[k]=(psi[k] - b)*std::exp( -std::pow( x, 4)) + b;
       
   }
}
//na razie funkcja rozumie tylko warunki brzegowe postaci:
// -stałe na początek i koniec

//test rozwiązania polegający na policzeniu jego  normy(nie kwadratu normy)
double solution::abs_val( )
{
    double wynik=0;
    for(int i =0; i< data_size ; ++i)
 	wynik+= std::norm(psi[i]) * step_spatial;
    return std::sqrt(wynik);
}//spoko loko
const std::vector<complex> solution::output_psi()
{
    std::vector<complex> out = std::vector<complex>(0);
    for(int i = 0; i < data_size; ++i )
        out.push_back(psi[i]);
    return out;    
}//spoko loko

const complex solution::output_psi(int i)
{
    return psi[i];
}
//spoko loko
const std::vector<double> solution::output_n_r()
{
    std::vector<double> out = std::vector<double>(0);
    for(int i = 0; i < data_size; ++i )
        out.push_back(n_r[i]);
    return out;    
}//spoko loko

const double solution::output_n_r(int i)
{
    return n_r[i];
}
//-------------------------PHYSICAL_CONSTANTS------------------------------
physical_constants::physical_constants()
{
    h_bar = 1.0;
    g_c = -1.0;
    g_r = 0.0;
    m = 0.5;
    lam_c = 0.0;
    lam_r = 0.0;
    R=0.0;
}

//-------------------FFTW3-------------------------------------------------

//fftw może używać tego samego planu jak długo rozmiar tablicy jest ten sam 
fourier_transform::~fourier_transform()
{
    fftw_destroy_plan( forward_plan );
    fftw_destroy_plan( backward_plan );
}

