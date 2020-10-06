function y = Par(x)
z = (pi/6)*x;
y = (1 - 6*z^2 + 6*z^3)*(z <= 0.5) + (2 * (1-z)^3)*(z>0.5)*(z<=1);
end