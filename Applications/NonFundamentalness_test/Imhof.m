function y = Imhof(u)

global w;
global x;


a = sin(0.5*sum(atan(w*u),1)-0.5*x*u);
b = prod((1+w.^2*u.^2).^0.25,1);
y = a./(u.*b);