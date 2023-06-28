% Write the test function (given in assignment 5 problem statement pdf) below to check correctness of your code.
% Define your function as follows:
%func=@(x)(x(1)-1)^2+(x(2)-1)^2-x(1)*x(2);
%func=@(x)x(1)^2+x(2)^2+(0.5*x(1)+x(2))^2+(0.5*x(1)+x(2))^4;
%func=@(x)(-0.0001* (abs( sin(x(1))*sin(x(2))*exp( abs( 100 - (sqrt( x(1)^2+x(2)^2 ))/pi ) ) ))^(0.1));
func=@(x)(x(1)+2*x(2)-7)^2+(2*x(1)+x(2)-5)^2;
x_initial=[1.5 1.5]';
% Also take intial guess
[ x_output, f_output, grad_output]=FR_CG(func, x_initial)


function [ x_output, f_output, grad_output]=FR_CG(func, x_initial)

%     -----------------------------------------------------------------------------------------------------------------------------
%     Write your logic for FR-CG using in-exact line search. 

%     Input parameters:  
%         func : input function to be evaluated
%         x_initial: initial value of x, a column vector
%     Returns:
%         x_output : converged x value, a column vector
%         f_output : value of f at x_output
%         grad_output : value of gradient at x_output, a column vector
%     -----------------------------------------------------------------------------------------------------------------------------       
%    Start your code here
    N=15000;
    epsilon=10^-6;

    X = x_initial;
    delta_fx0 = grad_f(X, func);
    delta_fx1 = delta_fx0;
    x_vals = zeros(N,2);
    f_vals = zeros(N,1);
    x_vals(1,:) = X;
    f_vals(1,:) = func(X);
    p = -delta_fx0;
    i = 1;
    while i<=N+1 && norm(delta_fx0)^2>epsilon
        delta_fx0 = delta_fx1;
        alpha = backtrack(delta_fx0, X, p, func);
        X = X + alpha*p;
        delta_fx1 = grad_f(X, func);
        beta = ((delta_fx1)'*delta_fx1)/((delta_fx0)'*delta_fx0);
        p = -delta_fx1 + beta*p;
        x_vals(i+1,:) = X;
        f_vals(i+1,:) = func(X);
        i = i+1;
    end
    if i == N+1
        fprintf('Maximum iterations reached but convergence did not happen.\n');
    end
    x_output = x_vals(i-1,:)';
    f_output=f_vals(i-1,:);
    grad_output = delta_fx0;
    figure
    plot(0:i-2,x_vals(1:i-1,1),'r-o',0:i-2,x_vals(1:i-1,2),'b-o')
    xlabel('Iteration Number')
    ylabel('x values')
    legend('x(1)','x(2)')

    figure
    plot(0:i-1,f_vals(1:i,1),'r-o')
    xlabel('Iteration Number')
    ylabel('f values')
    legend('f')

    function delF=grad_f(xstar, fun) % write your logic for gradient computing function below
        delF = [0;0];
        h = 0.001;
        delF(1,1) = (fun(xstar+[h;0])-fun(xstar-[h;0]))/(2*h);
        delF(2,1) = (fun(xstar+[0;h])-fun(xstar-[0;h]))/(2*h);
    end
    
    function alpha=backtrack(delF, xk, pk,func)
        alpha0 = 5;
        rho = 0.8;
        c = 0.1;
        alpha = alpha0;
        while(func(xk + alpha*pk) > func(xk) + c*alpha*(delF'*pk))
            alpha = alpha*rho;
        end
    end
end