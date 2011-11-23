# Y=gauss(X,p,t,sig0)
# t -time 
# X - arguments space
# p - momentum at t = 0
# m = 1.0
# h = 1.0
#
# d - half width(?)
function Y=gauss(X, p, t, sig0 )
Y=[];
#just for simpliyfying equations
sig = sig0* sqrt(1.0 + (t.^2)/(4.0 * sig0.^4) ); #ok

Y=exp( -((X - 2.0*p*(sig0.^2) ).^2) / (4.0*(sig0.^2 + i*t*0.5 ) ) ); #ok
Y=Y.* exp( -(sig0*p).^2 )/ sqrt( sig0/sqrt(2.0 * pi) ); #ok
Y=Y./sqrt(sig0.^2 + 0.5*i*t ); #ok
 
endfunction 
