x = [1;1];
del2f=compute_hessian(x)

function delF=compute_hessian(x)
    delF=zeros(2,2);
    h = 0.001;
    f_x = func(x);
    delF(1,1) = (func(x+[h;0]) - 2*f_x + func(x-[h;0]))/(h^2);
    delF(2,2) = (func(x+[0;h]) - 2*f_x + func(x-[0;h]))/(h^2);
    A = func(x+[h;h]);
    B = func(x-[h;h]);
    C = func(x+[-h;h]);
    D = func(x+[h;-h]);
    delF(1,2) = (A+B-C-D)/(4*h^2);
    delF(2,1) = delF(1,2);

    function f=func(X)
    f=2*X(2,1)*exp(X(1,1))+3*X(1,1)*(X(2,1)^2);
    end
end