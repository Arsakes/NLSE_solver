//
//		PROGRAM MA SŁUŻYĆ JAKO MODUŁ DO OCTAVE
//
//   Ewolucja wedle nieliniowego równania Sh, z uwzględniniem zmiennej liczby cząstek
//
//
//
#include<vector>
#include<complex>
#include<cmath>

//raportowanie błędów
#include<sstream>

//transformata fouriera
#include <fftw3.h>
//PRZESTRZEŃ NAZW - chyba nie potrzebna jak na razie
typedef std::complex<double> std_complex;



enum direction
{
   forward =  -1,
   backward =  1
};

//--------------NAMESPACE::nan_exception----------------------------------------------
namespace num_exception
{

enum Enum
{
   nan_psi = 1,
   nan_n_a = 2,
   nan_P_l = 4,
   nan_V   = 8,
};

struct exception
{
   std::string report;
   int error_state;
};

//sprawdza czy liczba jest skończona (tzn nie jest nieskończona albo NaN)
bool is_finite( std_complex  x );

}
//--------------ENDOF_NAMESPACE::nan_exception---------------------------------------

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
    //FIXME -- czy coś takiego rzutowanie jest bezpieczne, może da się lepiej?
    void inline init( std_complex* in , std_complex* out, int size)
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
   double g;
   double g_r;
   double lam_c;
   double lam_a;
   double lam_i;
   double m;
   double R;
   double tau; //czas życia nieaktywnych ekscytonów
   physical_constants();
};

class solution
{
  private:
    //zmienne przechowujące wartość rozwiązania dla ustalonego czasu "t" 
    //warunki brzegowe
    //std_complex boundary_begin;
    //std_complex boundary_end;
    
    //FIXME - może jest bezpieczniejszy sposób żeby pogodzić fftw i wskaźniki?
    std_complex* volatile psi;
    
    std::vector<double>  n_a;
    std::vector<double>  n_i; 
    
    //DANE RÓWNANIA 
    //długość kroku
    double step_spatial;
    double step_temporal; //krok czasowy może być urojony odpowiada to operatorowi nie
    int data_size;
    std::vector<std_complex> V;
    std::vector<double> P_l;
 
    //STAŁE FIZYCZNE
    physical_constants CONST;
  
    fourier_transform FT;
    
    
    //to zanczy dzielenie przez zero albo wyjście poza zakres
   
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
    void apply_boundary_condition(std_complex start, std_complex end,  double falloff);
    
    //NUMER KROKU
    int step;
    
  public:
    solution(const std::vector<std_complex>& temp, const std::vector<std_complex>& temp_V, 
          const std::vector<double>& temp_na, const std::vector<double>& temp_ni, const std::vector<double>& temp_p,
          double xstep, double tstep, physical_constants C );
    ~solution();
    void evolution(int steps  );
    const std::vector<std_complex> output_psi();
    const std::vector<double> output_n_a();
    const std::vector<double> output_n_i();
    const std_complex output_psi(int i);
    const double output_n_a(int i );
    const double output_n_i(int i );
    double abs_val();

    //FUNKCJA RAPORTUJĄCA O WSZELKICH WARTOŚCIACH NAN
    const num_exception::exception report_exception();
};
