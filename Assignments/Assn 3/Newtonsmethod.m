%fun=@(x)x(1)^2+x(2)^2+(0.5*x(1)+x(2))^2+(0.5*x(1)+x(2))^4;
%fun=@(x)(x(1)+2*x(2)-7)^2+(2*x(1)+x(2)-5)^2;
%fun=@(x)(-0.0001* (abs( sin(x(1))*sin(x(2))*exp( abs( 100 - (sqrt( x(1)^2+x(2)^2 ))/pi ) ) ))^(0.1));
fun=@(x)(x(1)-1)^2+(x(2)-1)^2-x(1)*x(2);
[X,F_val,delta_fx]= Newtons_method(fun)

function [X,F_val,delta_fx]= Newtons_method(fun)
    N=15000;
    x_intguess=[1.5 1.5]';
    epsilon=10^-6;

    % Newton Method
    X = x_intguess;
    delta_fx = grad_f(X);
    x_vals = zeros(N,2);
    x_vals(1,:) = X;
    i = 1;
    while i<=N+1 && norm(delta_fx)>epsilon
        delta_fx=grad_f(X);
        del2F=Hessian_f(X);
        X=X-del2F\delta_fx;
        x_vals(i+1,:) = X';
        i = i+1;
    end
    if i == N+1
        fprintf('Maximum iterations reached but convergence did not happen.\n');
    end
    F_val=fun(X);
    figure
    plot(0:i-1,x_vals(1:i,1),'r-o',0:i-1,x_vals(1:i,2),'b-o')
    xlabel('Iteration Number')
    ylabel('x values')
    legend('x(1)','x(2)')

    function delF=grad_f(xstar) % write your logic for gradient computing function below
        delF = [0;0];
        h = 0.001;
        delF(1,1) = (fun(xstar+[h;0])-fun(xstar-[h;0]))/(2*h);
        delF(2,1) = (fun(xstar+[0;h])-fun(xstar-[0;h]))/(2*h);
    end
    
    function del2F=Hessian_f(xstar) % wrie your logic for Hessian computing function below
        del2F=zeros(2,2);
        h = 0.001;
        f_x = fun(xstar);
        del2F(1,1) = (fun(xstar+[h;0]) - 2*f_x + fun(xstar-[h;0]))/(h^2);
        del2F(2,2) = (fun(xstar+[0;h]) - 2*f_x + fun(xstar-[0;h]))/(h^2);
        A = fun(xstar+[h;h]);
        B = fun(xstar-[h;h]);
        C = fun(xstar+[-h;h]);
        D = fun(xstar+[h;-h]);
        del2F(1,2) = (A+B-C-D)/(4*h^2);
        del2F(2,1) = del2F(1,2);
    end
end