
function [U] = heath_verkeed(k)
    n=100;
    dx = 0.01/n;
    A = zeros(n^2);
    k = k*ones(n,n);
    mu = 1/dx^2;

    for x = 1:(n)
        for y = 1:(n)
            if x==1
                if y==1
                    %linker rand, bovenhoek
                    A((y-1)*n+x,(y-1)*n+x) = -mu*(2*k(y,x)+k(y,x+1)+k(y+1,x))/2;
                    A((y-1)*n+x,(y)*n+x) = (k(y+1,x)+k(y,x))*mu/2;
                    A((y-1)*n+x,(y-1)*n+(x+1)) = (k(y,x+1)+k(y,x))*mu/2;
                else
                    if y==n
                        %linker rand, onder hoek
                        A((y-1)*n+x,(y-2)*n+x) = (k(y-1,x)+k(y,x))*mu/2;
                        A((y-1)*n+x,(y-1)*n+x) = -mu*(k(y-1,x)+2*k(y,x)+k(y,x+1))/2;
                        A((y-1)*n+x,(y-1)*n+(x+1)) = (k(y,x+1)+k(y,x))*mu/2;
                    else
                        if y*dx>0.003 && (y-1)*dx<0.007
                            A((y-1)*n+1,(y-1)*n+1) = -1;
                        else
                            %linker rand, geen hoeken
                            A((y-1)*n+x,(y-2)*n+x) = (k(y-1,x)+k(y,x))*mu/2;
                            A((y-1)*n+x,(y-1)*n+x) = -mu*(k(y-1,x)+3*k(y,x)+k(y,x+1)+k(y+1,x))/2;
                            A((y-1)*n+x,(y)*n+x) = (k(y+1,x)+k(y,x))*mu/2;
                            A((y-1)*n+x,(y-1)*n+(x+1)) = (k(y,x+1)+k(y,x))*mu/2;
                        end
                    end
                end
            else
                if x==n
                    if y==n
                            %rechter rand, onderhoek
                            A((y-1)*n+x,(y-2)*n+x) = (k(y-1,x)+k(y,x))*mu/2;
                            A((y-1)*n+x,(y-1)*n+(x-1)) = (k(y,x-1)+k(y,x))*mu/2;
                            A((y-1)*n+x,(y-1)*n+x) = -mu*(k(y-1,x)+2*k(y,x)+k(y,x-1))/2;
                    else
                        if y==1
                            %rechter rand,boven hoek
                            A((y-1)*n+x,(y-1)*n+(x-1)) = (k(y,x-1)+k(y,x))*mu/2;
                            A((y-1)*n+x,(y-1)*n+x) = -mu*(2*k(y,x)+k(y,x-1)+k(y+1,x))/2;
                            A((y-1)*n+x,(y)*n+x) = (k(y+1,x)+k(y,x))*mu/2;
                        else
                            %rechter rand, geen hoeken
                            A((y-1)*n+x,(y-2)*n+x) = (k(y-1,x)+k(y,x))*mu/2;
                            A((y-1)*n+x,(y-1)*n+(x-1)) = (k(y,x-1)+k(y,x))*mu/2;
                            A((y-1)*n+x,(y-1)*n+x) = -mu*(k(y-1,x)+3*k(y,x)+k(y,x-1)+k(y+1,x))/2;
                            A((y-1)*n+x,(y)*n+x) = (k(y+1,x)+k(y,x))*mu/2;
                        end
                    end
                else
                    if y==n
                        %onder rand, geen hoeken
                        A((y-1)*n+x,(y-2)*n+x) = (k(y-1,x)+k(y,x))*mu/2;
                        A((y-1)*n+x,(y-1)*n+(x-1)) = (k(y,x-1)+k(y,x))*mu/2;
                        A((y-1)*n+x,(y-1)*n+x) = -mu*(k(y-1,x)+3*k(y,x)+k(y,x-1)+k(y,x+1))/2;
                        A((y-1)*n+x,(y-1)*n+(x+1)) = (k(y,x+1)+k(y,x))*mu/2;
                    else
                        if y ==1
                            %boven rand, geen hoeken
                            A((y-1)*n+x,(y-1)*n+(x-1)) = (k(y,x-1)+k(y,x))*mu/2;
                            A((y-1)*n+x,(y-1)*n+x) = -mu*(3*k(y,x)+k(y,x-1)+k(y,x+1)+k(y+1,x))/2;
                            A((y-1)*n+x,(y)*n+x) = (k(y+1,x)+k(y,x))*mu/2;
                            A((y-1)*n+x,(y-1)*n+(x+1)) = (k(y,x+1)+k(y,x))*mu/2;
                        else
                            %normaal geval
                            A((y-1)*n+x,(y-2)*n+x) = (k(y-1,x)+k(y,x))*mu/2;
                            A((y-1)*n+x,(y-1)*n+(x-1)) = (k(y,x-1)+k(y,x))*mu/2;
                            A((y-1)*n+x,(y-1)*n+x) = -mu*(k(y-1,x)+4*k(y,x)+k(y,x-1)+k(y,x+1)+k(y+1,x))/2;
                            A((y-1)*n+x,(y)*n+x) = (k(y+1,x)+k(y,x))*mu/2;
                            A((y-1)*n+x,(y-1)*n+(x+1)) = (k(y,x+1)+k(y,x))*mu/2;
                        end
                    end
                end
            end
        end
    end
    A=-A;
    Q=ones(n^2,1);
    Q = 2.5*Q/(0.01^2);
    for y=1:n
        if y*dx>0.003 && (y-1)*dx<0.007
            Q((y-1)*n+1) = 300;
        end
    end
    A = sparse(A);
    U = A\Q;
    U = reshape(U,[n,n]);
end