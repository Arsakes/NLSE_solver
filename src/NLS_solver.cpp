#include<octave/oct.h>
#include"solve.h"
#include"input_test.h"

//NLS_solver( PSI, V, na0, ni0, P, T, CONSTS)
// 1) sprawdzanie wejścia
// 2) przepisywanie wejścia na typy c++
// 3) symulacja
// 4) przepisywanie wyjścia z typów c++ na typy octave'a

DEFUN_DLD( NLS_solver, args, nargout , "Non-linear Schrodinger equation solver." )
{
//------------------------CZĘŚĆ INTERAKCJI Z UŻYTKOWNIKIM------------------------
    ComplexNDArray psi0;
    ComplexNDArray V;
    NDArray P_l; 			
    NDArray n_a0;
    NDArray n_i0;
    NDArray time_steps;
    NDArray constants;
    //format stałych  [h_bar, g, g_r, lam_c, lam_a, lam_i, tau,  m, R, xstep,tstep] ;
    
    //args - lista argumentów 
    //    args(0) - wektor wierszowy warunek początkowy psi, zespolony wektor
    //    args(1) - wektor wierszowy potencjał,              zespolony wektor
    //    args(2) - wektor wierszowy warunek początkowy n_a  rzeczywisty wektor
    //    args(3) - wektor wierszowy warunek początkowy n_i  rzeczywisty wektor
    //    args(4) - wektor wierszowy oświetlanie przez laser rzeczywisty wektor
    //    args(5) - wektor wierszowy kroki czasowe           rzeczywisty wektor
    //    args(6) - wektor wierszowy stałe w problemie,      rzeczywisty wektor
    
//----------------------KONTROLA WEJŚCIA ZE WZGLĘDU NA TYP--------------------------
    input_test_val error_state_str;
    if(error_state)
    { 
        error("some unidentyfied error");
        return octave_value_list();
    }
 
    error_state_str = type_test( args(0), complex_array );
    if(  error_state_str.value ) 
    {
        error( ("psi0:" +error_state_str.error_code).c_str() );
        return octave_value_list();
    }
    else psi0 = args(0).complex_array_value();
  
    error_state_str = type_test( args(1), complex_array );
    if(  error_state_str.value )	 {
        error( ("V: "+error_state_str.error_code).c_str() );
        return octave_value_list();
    }
    else V = args(1).complex_array_value();
 
    error_state_str = type_test( args(2), real_array );
    if(  error_state_str.value )	{
        error( ("n_a0: "+error_state_str.error_code).c_str() );
        return octave_value_list();
    }
    else n_a0 = args(2).array_value();

    error_state_str = type_test( args(3), real_array );
    if(  error_state_str.value )	{
        error( ("n_i0: "+error_state_str.error_code).c_str() );
        return octave_value_list();
    }
    else n_i0 = args(3).array_value();

    error_state_str = type_test( args(4), real_array );
    if(  error_state_str.value ) 	 {
        error( ("P_l: " + error_state_str.error_code  ).c_str()  );
        return octave_value_list();
    }
    else P_l = args(4).array_value();

    error_state_str = type_test( args(5), real_array );
    if(  error_state_str.value ) 	{
        error( ("time_steps: "+error_state_str.error_code).c_str()  );
        return octave_value_list();
    }
    else  time_steps = args(5).array_value();
 
    error_state_str = type_test( args(6), real_array );
    if(  error_state_str.value )	{
        error( ("physical_constants: "+error_state_str.error_code).c_str()  );
        return octave_value_list();
    }
    else constants = args(6).array_value();

    
    //przepisywanie do wyjścia
    int spatial_size = psi0.nelem();
    if ( spatial_size != P_l.nelem() || spatial_size != n_a0.nelem() || spatial_size != V.nelem() || spatial_size != n_i0.nelem() ) 
         error("Inconsistent size of vector arguments.");
    if ( constants.nelem() != 11 ) 
         error("Bad number of simulation constants provided");

    //punkt czasowy z którego zaczynamy

    std_complex x; 
    std::vector<std_complex> in_psi0( 0);
    std::vector<std_complex> in_V(0);
    std::vector<double>  in_n_a0(0);
    std::vector<double>  in_n_i0(0);
    std::vector<double>  in_P_l(0);
    physical_constants in_consts;

    for(int i =0; i< spatial_size ;i++)  
    {
        in_psi0.push_back( std_complex( psi0(i).real(), psi0(i).imag() ));
        in_V.push_back( std_complex( V(i).real(), V(i).imag() ));
        in_n_a0.push_back( n_a0(i) );
        in_n_i0.push_back( n_i0(i) );
        in_P_l.push_back( P_l(i) );
    }
    in_consts.h_bar = constants(0);
    in_consts.g = constants(1);
    in_consts.g_r = constants(2);
    in_consts.lam_c = constants(3);
    in_consts.lam_a = constants(4);
    in_consts.lam_i = constants(5);
    in_consts.tau = constants(6);
    in_consts.m = constants(7);
    in_consts.R = constants(8);
    double xstep = constants(9);
    double tstep = constants(10);
    
//------------------------KONTROLA WEJŚCIA ZE WZGLĘDU NA ZAKRES WARTOŚCI-------------------
if(xstep == 0) 
{
    error("Spatial step must be positve value.");
    return octave_value_list();
}
if(tstep == 0)
{
    error("Time step must be positve value.");
}
if(constants(0) == 0) 
{
    error(" 'h_bar' must be non-zero value.");
    return octave_value_list();
}
if(constants(5) == 0)
{
    error(" 'm' must be non-zero value.");
    return octave_value_list();
}


//------------------CZĘŚĆ SYMULACYJNA FUNKCJI--------------------
    
     solution NLS_solver=solution( in_psi0, in_V, in_n_a0, in_n_i0, in_P_l,  xstep, tstep, in_consts  );   

     //FIXME -- kompilator mówi że taka deklaracja jest przestarzała
     dim_vector dv (2);
     dv(0) = spatial_size; 
     dv(1) = time_steps.nelem();
     octave_value_list retval;
     ComplexNDArray output_psi(dv);
     ComplexNDArray output_n_a(dv);
     ComplexNDArray output_n_i(dv);
     double time = 0.0;
     
     //ZMIENNA PRZECHOWUJĄCA INFORMACJE O BŁĘDACH WEWNĄTRZ SYMULACJI 
     num_exception::exception sim_error;

     for( int index=0; index < time_steps.nelem(); index++)
     { 
         //raportuj o wszystkich napotkanych nieskończonościach i NANach
         sim_error = NLS_solver.report_exception();
         if( sim_error.error_state )	
         {
             //moze disp?
             printf( sim_error.report.c_str() );
         }
         
         double delta_time = time_steps(index) - time ; 
         // dopuszczone sa puste iteracje (ciało symulacji nie jest wykonywane a krok nie jest
         // liczony, tak się zdarza dla chujowo małych czasów
         //iteruje aż do momentu wyplucia danych w momencie zadanym przez tablice time_steps
         NLS_solver.evolution( std::floor( delta_time/ std::abs(tstep)  ) );
            
         for(int i=0; i< spatial_size; ++i)
         { 
            output_psi(i + spatial_size * index) = NLS_solver.output_psi(i);
            output_n_a(i + spatial_size * index) = NLS_solver.output_n_a(i);
            output_n_i(i + spatial_size * index) = NLS_solver.output_n_i(i);
         }
        
        //żeby poprawnie liczyć czas
        time += std::floor(delta_time/std::abs(tstep)) * std::abs(tstep);
    }
  
//--------------------------CZĘŚĆ INTERAKCJI WYJŚCIE--------------------------------------
    retval.append( output_psi );
    retval.append( output_n_a );
    retval.append( output_n_i );
    
    return retval;
}

