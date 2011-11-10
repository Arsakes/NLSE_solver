#FUNNY BUNNY
function return_val = soliton(x, v, t)
    return_val = sqrt(2)*exp(0.5*i*v.*x + i*(1 - 0.25* x.^2 )*t ).* sech(x.- v*t);
endfunction
