function x_answer = hw07_conjugate(A, b, n, epsilon, kmax)

k = 0;
for i = 1:n
    x(i) = 0;
end
x_answer = x';
r = b;
rho(1) = norm(r)^2;

while(sqrt(rho(k+1)) > epsilon*norm(b) && k<kmax)
    k = k+1;
    if(k==1)
        p = r;
    else
        Beta(k+1) = rho(k)/rho(k-1);
        p = r + Beta(k+1)*p;
    end
    w = A*p;
    alpha(k+1) = rho(k)/((p')*w);
    x_answer = x_answer + alpha(k+1)*p;
    r = r - alpha(k+1)*w;
    rho(k+1) = norm(r)^2;
end
