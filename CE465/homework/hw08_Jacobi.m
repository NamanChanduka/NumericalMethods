function[eigenval, k, offList] = hw08_Jacobi(A, n, eps)
k = 0;
V = eye(n);

[max_values, Index]=max(abs(A-diag(diag(A))));
[~, q]=max(max_values);
p = Index(q);                                                               %Finding indices of the absolute max term in the matrix i.e p and q
    
offnorm = sqrt(sum(sum((A-diag(diag(A))).^2)));                             %Offdiagonal norm of the matrix
offList = offnorm;
frobenius = sqrt(sum(sum(A.^2)));                                           %Frobenius norm of the matrix

while(offnorm > eps*frobenius)
    k = k + 1;                                                              %Tracking number of iterations
    
    if ( A(p,q) ~= 0)
        tau = (A(q,q) - A(p,p))/(2*A(p,q));
        t = sign(tau)/(abs(tau) + sqrt(1 + tau^2));
        c = 1/(sqrt(1+t^2));
        s = t*c;
    else
        c = 1;
        s = 0;                                                              %Calculating the c and s for the J(p,q,theta) matrix
    end
    J = eye(n); 
    J(p,p) = c;     
    J(q,q) = c;
    J(p,q) = s;
    J(q,p) = -s;                                                            %Defining J matrix
    A = J'*A*J;
    V = V*J;                                                                %Finding next and more diagonal similarity transformed matrix A
    
    offnorm = sqrt(sum(sum((A-diag(diag(A))).^2)));
    offList = [offList; offnorm];
    
    [max_values, Index]=max(abs(A-diag(diag(A))));
    [~, q]=max(max_values);
    p = Index(q);
end

eigenval = diag(A);                                                         %When all offdiagonal terms reach 0, the diagonal elements are eigenvalues itself.