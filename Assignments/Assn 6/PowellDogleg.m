% Write the test function (given in assignment 5 problem statement pdf) below to check correctness of your code.
% Define your function as follows:
%func=@(x)(x(1)-1)^2+(x(2)-1)^2-x(1)*x(2);
func=@(x)x(1)^2+x(2)^2+(0.5*x(1)+x(2))^2+(0.5*x(1)+x(2))^4;
%func=@(x)(-0.0001* (abs( sin(x(1))*sin(x(2))*exp( abs( 100 - (sqrt( x(1)^2+x(2)^2 ))/pi ) ) ))^(0.1));
%func=@(x)abs(x(1))^2+abs(x(2))^3;
%func=@(x)(x(1)+2*x(2)-7)^2+(2*x(1)+x(2)-5)^2;
x_initial=[1.5 1.5]';
% Also take intial guess
[x_output, f_output, grad_output]=TRPD(func, x_initial)


function [ x_output, f_output, grad_output]=TRPD(func, x_initial)

%     -----------------------------------------------------------------------------------------------------------------------------
%     Write your logic for Powell Dogleg Trust Region method. 

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
    delta_fx = grad_f(X, func);
    x_vals = zeros(N,2);
    f_vals = zeros(N,1);
    x_vals(1,:) = X;
    f_vals(1,:) = func(X);
    delta_bar = 1;
    delta = 0.5;
    eta = 0.2;
    i = 1;
    while i<=N+1 && norm(delta_fx)^2>epsilon
        delta_fx = grad_f(X,func);
        hess_fx = Hessian_f(X);
        
        ps = -delta * delta_fx / norm(delta_fx);
        
        b = delta_fx'*hess_fx*delta_fx;
        if(b <= 0)
            tau = 1;
        else
            tau = (delta_fx'*delta_fx)*norm(delta_fx) / (delta * b);
            if(norm(tau*ps) >= delta)
                tau = 1;
            end
        end

        pc = ps * tau;

        pn = -inv(hess_fx) * delta_fx;

        if(norm(pn)<=delta)
            p = pn;
        elseif(norm(pc)>=delta)
            p = delta * pc / norm(pc);
        else
            a = (-(pn-pc)'*pc + sqrt(((pn-pc)'*pc)^2+(delta^2 - (norm(pc))^2)*(norm(pn-pc))^2)) / (norm(pn-pc))^2;
            p = a*pn + (1-a)*pc;
        end

        rho = (f_vals(i)-func(X+p))/(-delta_fx'*p - p'*hess_fx*p);
        if(rho<0.25)
            delta = 0.25*norm(p);
        else
            if rho>0.75 && norm(p)==delta
                delta = min(2*delta, delta_bar);
            end
        end
        if rho>eta
            X = X + p;
        end
        x_vals(i+1,:) = X;
        f_vals(i+1) = func(X);
        i = i+1;
    end
    if i == N+1
        fprintf('Maximum iterations reached but convergence did not happen.\n');
    end
    x_output = x_vals(i-1,:);
    f_output= f_vals(i-1,:);
    grad_output = delta_fx;
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

    function del2F=Hessian_f(xstar) % write your logic for Hessian computing function below
        del2F=zeros(2,2);
        h = 0.001;
        f_x = func(xstar);
        del2F(1,1) = (func(xstar+[h;0]) - 2*f_x + func(xstar-[h;0]))/(h^2);
        del2F(2,2) = (func(xstar+[0;h]) - 2*f_x + func(xstar-[0;h]))/(h^2);
        A = func(xstar+[h;h]);
        B = func(xstar-[h;h]);
        C = func(xstar+[-h;h]);
        D = func(xstar+[h;-h]);
        del2F(1,2) = (A+B-C-D)/(4*h^2);
        del2F(2,1) = del2F(1,2);
    end
    
%     End your code here
end