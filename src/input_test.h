//PLIK Z FUNKCJĄ SPRAWDZAJĄCĄ POPRAWNOŚĆ WEJŚCIA
enum data_type
{
    real_scalar,
    complex_scalar,
    real_array,
    complex_array,

};

struct input_test_val
{
    bool value;
    std::string error_code;
};

//sprawdza czy wartość to_be_test typu reference_type 
//zwraca: kod błędu i wartość logiczną
input_test_val type_test(const octave_value&  to_be_test, data_type reference_type );
