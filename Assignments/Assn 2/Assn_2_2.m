x = [1;1];
f_grad = compute_gradient(x)

function delF=compute_gradient(x)
    delF = [0;0];
    h = 0.001;
    delF(1,1) = (func(x+[h;0])-func(x-[h;0]))/(2*h);
    delF(2,1) = (func(x+[0;h])-func(x-[0;h]))/(2*h);

    function f=func(X)
    f=2*X(2,1)*exp(X(1,1))+3*X(1,1)*(X(2,1)^2);
    end
end