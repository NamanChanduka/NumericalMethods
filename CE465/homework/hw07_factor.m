function [A, pivot, determinant, ier] = hw07_factor(A, n, pivot)

determinant = 1;
k = 1;                                                                      %looping variable initialized
ier = -1;
while k < n
    Ck = abs(A(k,k));                                                       %initialize Ck
    i0 = k;                                                                 %initialize index to k
    partialpivot = k;
    while partialpivot < n+1                                                %looping to find max value and index
        if abs(A(partialpivot,k)) > Ck
            Ck = A(partialpivot,k);
            i0 = partialpivot;
        end
        partialpivot = partialpivot+1;
    end
    pivot(k) = i0;                                                          %pivot = index of max val
    if abs(Ck) == 0                                                         %checking if max value of column is 0
        ier = 1;                                                            %changing error flag to 1
        determinant = 0;
        return;
    end
    if i0 ~= k                                                              %checking if index of max value is equal to current value of looping variable
        determinant = -determinant;
        j = k-1;
        while j < n+1                                                       %loop to interchange rows.
            t = A(k, j);
            A(k, j) = A(i0, j);
            A(i0, j) = t;
            j = j+1;
        end
    end
    for i = k+1:n
        m = A(i,k)/A(k,k);                                                  %finding row multiplier
        A(i,k) = m;
        for j = k+1:n
            A(i,j) = A(i,j) - m*A(k,j);
        end
    end
    determinant = A(k, k) * determinant;
    k = k+1;
end
determinant = determinant * A(n,n);
ier = 0;                                                                    %changing error flag to 0
end