function b = hw06_solve(A, n, b, pivot)


for k = 1:n-1
    i = pivot(k);
    if pivot(k) ~= k                                                    %perform pivoting on b matrix according to Gauss Elim
        temp = b(i);
        b(i) = b(k);
        b(k) = temp;
    end
    for i = k+1:n                                                       %copying matrix calculations from A using multipliers stored in A since it is already factorised
        b(i) = b(i) - A(i,k)*b(k);
    end
end
b(n) = b(n)/A(n,n);                                                     %storing value of last variable x8 and then using back substitution for others
for i = (n-1):-1:1
    sum = 0;
    for j = i+1:n
        sum = sum + A(i ,j)*b(j);
    end
    b(i) = (1/A(i,i))*(b(i) - sum);                                     %finding x in Ax = b using back substitution.
end
end