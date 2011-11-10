#include<octave/oct.h>
#include<string>
#include"input_test.h"

input_test_val type_test(const octave_value&  to_be_test, data_type reference_type )
{   
    input_test_val A = input_test_val();
    A.value = true;
    switch( reference_type )
    {
        case real_scalar:
             A.value = !to_be_test.is_real_scalar(); 
             if(A.value)      
                 A.error_code="Argument isn't real scalar";
        break;
        case complex_scalar:
             A.value = !(to_be_test.is_complex_scalar() );
             if(A.value) 
                 A.error_code="Argument isn't complex scalar";
        break;
        case complex_array:
	     A.value = !(to_be_test.is_complex_type() && (to_be_test.rows() > 1 || to_be_test.columns() > 1) );
             if(A.value)
                 A.error_code="Argument isn't complex array";
        break;
        case real_array:  
             A.value = !(to_be_test.is_real_type() && (to_be_test.rows() > 1 || to_be_test.columns() > 1) );
             if(A.value)
                 A.error_code="Argument isn't real array";
        break;
    }

    return A;
}

