x = [1;1];
fvalue=func(x)

function y=func(x)
    y=2*x(2,1)*exp(x(1,1))+3*x(1,1)*(x(2,1)^2);
end