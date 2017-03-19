function [U] = heath_correct_complex_2d(k)
    n = 50;
%     k = k*ones(1,n^2);
k = reshape(k,[1,n^2]);
    A = zeros(n^2,n^2);
    delta_x = 0.01/n;
    Q = 2.5*ones(n^2,1)/0.01^2;
    for i = n+1:n^2-n
        A(i,i-n) = k(i) + k(i-n);
        A(i,i-1) = k(i) + k(i-1);
        A(i,i) = -(4*k(i) + k(i-n) + k(i-1) + k(i+1) + k(i+n));
        A(i,i+1) = k(i) + k(i+1);
        A(i,i+n) = k(i) + k(i+n);
    end
    
    % RVW 1: linkerrand geisoleerd
    for i = 2:n-1
        A(i,:) = zeros(1,n^2);
        A(i,i-1) = k(i) + k(i-1);
        A(i,i) = -(3*k(i) + k(i-1) + k(i+1) + k(i+n));
        A(i,i+1) = k(i) + k(i+1);
        A(i,i+n) = k(i) + k(i+n);
    end
    
    % RVW 2: onderrand geisoleerd
    for i = n+1:n:n*(n-2)+1
        A(i,:) = zeros(1,n^2);
        A(i,i-n) = k(i) + k(i-n);
        A(i,i) = -(3*k(i) + k(i-n) + k(i+1) + k(i+n));
        A(i,i+1) = k(i) + k(i+1);
        A(i,i+n) = k(i) + k(i+n);
    end

    % RVW 3: rechterrand geisoleerd
    for i = n^2-n+2:n^2-1
        A(i,:) = zeros(1,n^2);
        A(i,i-n) = k(i) + k(i-n);
        A(i,i-1) = k(i) + k(i-1);
        A(i,i) = -(3*k(i) + k(i-n) + k(i-1) + k(i+1));
        A(i,i+1) = k(i) + k(i+1);
    end
   
    % RVW 4: bovenrand geisoleerd
    for i = 2*n:n:n*(n-1)
        A(i,:) = zeros(1,n^2);
        A(i,i-n) = k(i) + k(i-n);
        A(i,i-1) = k(i) + k(i-1);
        A(i,i) = -(3*k(i) + k(i-n) + k(i-1) + k(i+n));
        A(i,i+n) = k(i) + k(i+n);
    end
    
    % RO
    i = n*(n-1)+1;
    A(i,:) = zeros(1,n^2);
    A(i,i-n) = k(i) + k(i-n);
    A(i,i) = -(2*k(i) + k(i-n) + k(i+1));
    A(i,i+1) = k(i) + k(i+1);
    % RB
    i = n^2;
    A(i,:) = zeros(1,n^2);
    A(i,i-n) = k(i) + k(i-n);
    A(i,i-1) = k(i) + k(i-1);
    A(i,i) = -(2*k(i) + k(i-n) + k(i-1));
    % LO
    i = 1;
    A(i,:) = zeros(1,n^2);
    A(i,i) = -(2*k(i) + k(i+1) + k(i+n));
    A(i,i+1) = k(i) + k(i+1);
    A(i,i+n) = k(i) + k(i+n);
    % LB
    i = n;
    A(i,:) = zeros(1,n^2);
    A(i,i-1) = k(i) + k(i-1);
    A(i,i) = -(2*k(i) + k(i+n) + k(i-1));
    A(i,i+n) = k(i) + k(i+n);
    
    A = -A./(2*delta_x^2);
    
    % Randvoorwaarde 1: linkerrand = T
    start = 0.003/delta_x+1;
    stop = 0.007/delta_x-1;
    for i = start:stop
        A(i,:) = zeros(1,n^2);
        A(i,i) = 1;
        Q(i) = 300;
    end
    A = sparse(A);
    U = A\Q;
    U = reshape(U,[n,n]);
    surf(U)
end