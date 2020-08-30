function y = COVX(x,j,N,z)
for t = j+1:N
    z =  z + x(t,:)'*x(t-j,:);
end
z = z';
vech = z(triu(true(size(z))));
y = sum(vech)/(N-j);
end
