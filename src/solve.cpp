#include"solve.h"

//c++11:: SPRAWDŹ CZY MOŻLIWA JEST KONWERSJA
//static_assert( sizeof(fftw_complex) == sizeof(std_complex), "fftw complex not compatibile with std::complex<double>" );


const double PI = 3.1415926535;

solution::solution(const std::vector<std_complex>& temp, const std::vector<std_complex>& temp_V, 
     const std::vector<double>& temp_na, const std::vector<double>& temp_ni,
     const std::vector<double>& temp_p,
     double xstep =0.01, double tstep=0.001, physical_constants C = physical_constants()   )
{

   step_spatial=xstep;
   step_temporal=tstep;
   //inicjacja stałych fizycznych
   CONST.h_bar = C.h_bar;
   CONST.g_r = C.g_r;
   CONST.g = C.g;
   CONST.lam_c = C.lam_c;
   CONST.lam_a = C.lam_a;
   CONST.lam_i = C.lam_i;
   CONST.tau = C.tau;
   CONST.m = C.m;
   CONST.R = C.R;

   data_size = temp.size();
   V = std::vector<std_complex>(temp_V);
   n_a = std::vector<double>(temp_na);
   n_i = std::vector<double>(temp_ni);
   P_l = std::vector<double>(temp_p);

   //inicjacja fftw
   FT=fourier_transform();
   
   //FIXME -- czy ta konwersja jest bezpieczna? znaleźć inny sposób
   psi = reinterpret_cast<std_complex*>( fftw_malloc(sizeof(std_complex) * data_size) );
   FT.init( psi, psi, data_size);
   
   for(int i = 0; i < data_size; ++i)
   {
       psi[i] = temp[i];
   }
   
   step = 0;
   
}//spoko


solution::~solution()
{
   fftw_free(psi);
}
//spoko loko


void solution::first_step()
{
   //jedna część równania z pochodną przestrzenną jest
   //ewoluowana do następnej chwili czasowej
   FT.compute( forward );
   for(int i =0; i<data_size ; ++i)
   {
      //normalizacja
      psi[i] = psi[i] / std_complex(data_size,0.0);
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
       psi[i] = psi[i] * std::exp(std_complex(0.0,-k_square ) * step_temporal * 0.5*CONST.h_bar/CONST.m ) ;   
    }


    FT.compute( backward );
}//spoko

void solution::second_step()
{ 
    std_complex phase(0.0,0.0);
    for(int i = 0 ; i < data_size;  ++i)
    { 
        //EWOLUCJA "psi"
        //ewolucja napędzane nieliniowością i potencjałem
        phase = std_complex(0.0, CONST.g * std::norm(psi[i]) ); 
        phase += std_complex(0.0,1.0) * V[i] ;
        
        //ewolucja napędzana przez rozpraszanie i zmiany w ilości polarytonów
        phase += std_complex( 0.5* (CONST.lam_c - R(n_a[i])) , CONST.g_r*(n_a[i]+n_i[i]));

        //aplikacja do funkcji falowej
        psi[i] = psi[i] * std::exp( -phase * step_temporal);
        //spoko loko

        //EWOLUCJA "n_i" - czyli nieaktywnej populacji ekscytonów
        n_i[i] += ( -n_i[i]*(CONST.lam_i + 1.0/CONST.tau) + P_l[i] ) * step_temporal;
        //spoko loko
       
        //EWOLUCJA "n_a" - czyli aktywowanej populacji ekscytonów
        n_a[i] += ( -n_a[i] *(CONST.lam_a + CONST.R*std::norm(psi[i])) + n_i[i]/CONST.tau ) * step_temporal;
        //spoko loko
    }
}//spoko


//LEAP-FROG
void solution::make_time_step()
{
    //jak wskazywałaby nazwa krok drugi jest wykonywany pierwszy
    second_step();
    first_step();
    ++step;
}//spoko

void solution::evolution(int steps )
{
    for(int i = 0; i < steps; i++)
    {
            make_time_step();
            apply_boundary_condition(std_complex(0.0,0.0), std_complex(0.0,0.0), 0.1 );
    }
}//spoko

void solution::apply_boundary_condition(std_complex a=std_complex(0.0,0.0), std_complex b=std_complex(0.0,0.0), double c = 0.1)
{
   int falloff = std::floor( c * static_cast<double>(data_size) );
   int k;
   double x;
   //jak długim obszarze rociąga się "schodek" od 0 do 1
   
   //założenie 
   for( int i=0; i < falloff; ++i)
   { 
       //FIXME możliwe odwołanie się do nieistniejącego elementu, patrz niżej
       
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
 	wynik+= std::norm(psi[i]);
    return std::sqrt(wynik*step_spatial);
}//spoko loko

const std::vector<std_complex> solution::output_psi()
{
    std::vector<std_complex> out = std::vector<std_complex>(0);
    for(int i = 0; i < data_size; ++i )
        out.push_back(psi[i]);
    return out;    
}//spoko loko

const std_complex solution::output_psi(int i)
{
    return psi[i];
}//spoko loko

const std::vector<double> solution::output_n_a()
{
    return n_a;    
}//spoko loko

const std::vector<double> solution::output_n_i()
{
    return n_i;    
}//spoko loko

const double solution::output_n_a(int i)
{
    return n_a[i];
}//spoko loko

const double solution::output_n_i(int i)
{
    return n_i[i];
}//spoko loko

const num_exception::exception solution::report_exception()
{
   //FIXME -- czy takie operacje dają wynik niezależny od implementacji "int"
   int error_state = 0;
   std::ostringstream report;

   std::string comment("NAN or infinite value found in: ");
   //TODO -- dodać też sprawdzanie nan w n_i
   for(int i = 0;i<data_size; ++i)
   {
       if( num_exception::is_finite(n_a[i])==0 && (error_state & num_exception::nan_n_a) == 0)
       {
          error_state =  error_state | num_exception::nan_n_a;
          report << comment;
          report << std::string("n_a");
	  report << std::endl;
       }
       if( num_exception::is_finite(psi[i])==0 && (error_state & num_exception::nan_psi) == 0 )
       {
          error_state =  error_state | num_exception::nan_psi;
          report << comment;
          report << std::string("psi");
	  report << std::endl;
       }
       if( num_exception::is_finite(V[i])==0 && ( error_state & num_exception::nan_V) == 0 )
       {
          error_state =  error_state | num_exception::nan_V;
          report << comment;
          report << std::string("V");
	  report << std::endl;
       }
       if( num_exception::is_finite(P_l[i])==0 && (error_state & num_exception::nan_P_l) == 0 )
       {
          error_state =  error_state | num_exception::nan_P_l;
	  report << comment;
	  report << std::string("P_l");
	  report << std::endl;
       }
   }
   
   if(error_state == 0)
   {
       report << std::string("No exception found ");
   }
   
   report << std::string("in step: ");
   report << step;
   report << std::endl;

   num_exception::exception temp;
   temp.report = report.str();
   temp.error_state = error_state;
   
   return temp;   
}//spoko loko

//-------------------------PHYSICAL_CONSTANTS------------------------------
physical_constants::physical_constants()
{
    h_bar = 1.0;
    g = -1.0;
    g_r = 0.0;
    m = 0.5;
    lam_c = 0.0;
    lam_a = 0.0;
    lam_i = 0.0;
    tau = 1.0;
    R=0.0;  
}

//-------------------FFTW3---------------------------------------------------

//fftw może używać tego samego planu jak długo rozmiar tablicy jest ten sam 
fourier_transform::~fourier_transform()
{
    fftw_destroy_plan( forward_plan );
    fftw_destroy_plan( backward_plan );
}


//------------------num_exception namespace-----------------------------------
bool num_exception::is_finite( std_complex  x )
{
   if( std::isfinite( x.real() ) && std::isfinite( x.imag()  )  )
       return true;
   else 
       return false;
}

