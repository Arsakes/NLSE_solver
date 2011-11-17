//
//		PROGRAM MA SŁUŻYĆ JAKO MODUŁ DO OCTAVE
//
//   Ewolucja wedle nieliniowego równania Sh, z uwzględniniem zmiennej liczby cząstek
//
//
//

#include<vector>
#include<complex>

//funkcje matematyczne
#include<cmath>
//transformata fouriera
#include <fftw3.h>


enum direction
{
   forward =  -1,
   backward =  1
};


//PRZESTRZEŃ NAZW - chyba nie potrzebna jak na razie
typedef std::complex<double> complex;

//TYP DO MACIERZY LICZB RZECZYWISTYCH
//typedef std::vector< std::vector<double> > matrix_double;

//TRANSFORMACJA FOURIERA 
class fourier_transform
{
  private:
    //plan liczenia transformaty fouriera - liczy się przy pierwszej transformacie i jest 
    //używany wielokrotnie
    //wielkość dyskretnego zbioru argumentów
    //RZUTOWANIE
    fftw_plan forward_plan;
    fftw_plan backward_plan;
  public:
    void inline init( complex* in , complex* out, int size)
    {
           forward_plan = fftw_plan_dft_1d( size, reinterpret_cast<fftw_complex*>(in), 
	   		reinterpret_cast<fftw_complex*>(out), FFTW_FORWARD, FFTW_PATIENT);
           backward_plan = fftw_plan_dft_1d( size, reinterpret_cast<fftw_complex*>(in), 
	   		reinterpret_cast<fftw_complex*>(out), FFTW_BACKWARD, FFTW_PATIENT);
    }
    void inline compute( direction d)
    {
        if( d == forward )
	   fftw_execute( forward_plan );
	else if( d == backward )
	   fftw_execute( backward_plan);
    }
    ~fourier_transform();
};

//STAŁE FIZYCZNE OBECNE W PROBLEMIE, OPAKOWANE W STRUKTURE
//mają nadpisany konstruktor domyślny dla potrzeb testów
struct physical_constants
{
   double h_bar;
   double g_c;
   double g_r;
   double lam_c;
   double lam_r;
   double m;
   double R;
   physical_constants();
};

class solution
{
  private:
    //zmienne przechowujące wartość rozwiązania dla ustalonego czasu "t" 
    //warunki brzegowe
    //complex boundary_begin;
    //complex boundary_end;
    complex* psi;
    std::vector<double>  n_r;
    
    //DANE RÓWNANIA 
    //długość kroku
    double step_spatial;
    double step_temporal; //krok czasowy może być urojony odpowiada to operatorowi nie
    int data_size;
    std::vector<complex> V;
    std::vector<double> P_l;
 
    //STAŁE FIZYCZNE
    physical_constants CONST;
  
    fourier_transform FT;
    
    //FUNKCJA, która może być dostarczona z zewnątrz, ale na razie po prostu jest metodą
    //żeby być bardziej "explicit"
    double inline R( double s )
    { 
        return s*CONST.R;
    }
 
    //METODY
    void first_step();
    void second_step();
    void make_time_step();
    void apply_boundary_condition(complex start, complex end,  double falloff);
    
  public:
    solution(const std::vector<complex>& temp, const std::vector<complex>& temp_V, 
          const std::vector<double>& temp_n, const std::vector<double>& temp_p,
          double xstep, double tstep, physical_constants C );
    ~solution();
    void evolution(int steps  );
    const std::vector<complex> output_psi();
    const std::vector<double> output_n_r();
    const complex output_psi(int i);
    const double output_n_r(int i );
    double abs_val();
};
