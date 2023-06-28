x_initial = [1.5 1.5]';
func=@(x)(x(1)-1)^2+(x(2)-1)^2-x(1)*x(2);
%func=@(x)x(1)^2+x(2)^2+(0.5*x(1)+x(2))^2+(0.5*x(1)+x(2))^4;
%func=@(x)(-0.0001* (abs( sin(x(1))*sin(x(2))*exp( abs( 100 - (sqrt( x(1)^2+x(2)^2 ))/pi ) ) ))^(0.1));
%func=@(x)(x(1)+2*x(2)-7)^2+(2*x(1)+x(2)-5)^2;

[x_SD, f_SD, grad_SD, x_NM, f_NM, grad_NM, x_QN, f_QN, grad_QN]= iterative_methods(func, x_initial);


% ""=================================================== Assignment 4 ===================================================

% Some instructions:
%     * You can write seperate function for gradient and hessian computations.
%     * You can also write any extra function as per need.
%     * Use in-build functions for the computation of inverse, norm, etc. 




function [x_SD, f_SD, grad_SD, x_NM, f_NM, grad_NM, x_QN, f_QN, grad_QN]= iterative_methods(func, x_initial)


%%  A function to call your steepest descent, newton method and quasi-newton method.
[ x_SD, f_SD, grad_SD]=steepest_descent(func, x_initial)
[ x_NM, f_NM, grad_NM]=newton_method(func, x_initial)
[ x_QN, f_QN, grad_QN]=quasi_newton_method(func, x_initial)




end



function [ x_output, f_output, grad_output]=steepest_descent(func, x_initial)
%-----------------------------------------------------------------------------------------------------------------------------
%Write your logic for steepest descent using in-exact line search. 

%Input parameters: 
            %func : input function to be evaluated
            % x_initial: initial value of x, a column vector 

%Returns:
%         x_output : converged x value, a column vector 
%         f_output : value of f at x_output
%         grad_output : value of gradient at x_output, a column vector
 %  -----------------------------------------------------------------------------------------------------------------------------

 %   # Start your code here
epsilon = 10e-6;
 N = 15000;
 x_output = x_initial;
 grad_output = grad_f(x_output,func);
 num_iter = 0;
 x_vals = zeros(N+1,2);
 x_vals(1,:) = x_initial;
 f_vals = zeros(N+1);
 f_vals(1) = func(x_initial);
while(norm(grad_output) > epsilon && num_iter<N)
    num_iter=num_iter+1;
    grad_output = grad_f(x_output,func);
    p = -grad_output;
    alpha = backtrack(grad_output, x_output, p,func);
    x_output = x_output + alpha*p;
    x_vals(num_iter+1,:) = x_output;
    f_vals(num_iter+1) = func(x_output);
end
f_output = func(x_output);
if(num_iter==N && (norm(grad_output))^2 > epsilon)
    fprintf('Maximum iterations reached but convergence did not happen.\n');
end

figure(1)
plot(0:num_iter,x_vals(1:num_iter+1,1),'r-o',0:num_iter,x_vals(1:num_iter+1,2),'b-o')
title("Steepest Descent Method")
xlabel('Iteration Number')
ylabel('x values')
legend('x(1)','x(2)')

figure(2)
plot(0:num_iter,f_vals(1:num_iter+1,1),'r-o')
title("Steepest Descent Method")
xlabel('Iteration Number')
ylabel('f values')


    function delF=grad_f(xstar,func) % write your logic for gradient computing function below
        delF = [0;0];
        h = 0.001;
        delF(1,1) = (func(xstar+[h;0])-func(xstar-[h;0]))/(2*h);
        delF(2,1) = (func(xstar+[0;h])-func(xstar-[0;h]))/(2*h);
    end
        
    function del2F=Hessian_f(xstar,func) % wrie your logic for Hessian computing function below
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
    
    function alpha=backtrack(delF, xk, pk,func)
        alpha0 = 5;
        rho = 0.8;
        c = 0.1;
        alpha = alpha0;
        fxk = func(xk);
        while(func(xk + alpha*pk) > fxk + c*alpha*(delF'*pk))
            alpha = alpha*rho;
        end
    end
 %   # End your code here

end

function [ x_output, f_output, grad_output]=newton_method(func, x_initial)
%-----------------------------------------------------------------------------------------------------------------------------
% Write your logic for Newton method using in-exact line search. 

%Input parameters: 
            %func : input function to be evaluated
            % x_initial: initial value of x, a column vector

%Returns:
%         x_output : converged x value, a column vector 
%         f_output : value of f at x_output
%         grad_output : value of gradient at x_output, a column vector
 %  -----------------------------------------------------------------------------------------------------------------------------

 %   # Start your code here
epsilon = 10e-6;
 N = 15000;
 x_output = x_initial;
 grad_output = grad_f(x_output,func);
 num_iter = 0;
 x_vals = zeros(N,2);
 x_vals(1,:) = x_initial;
 f_vals = zeros(N+1);
 f_vals(1) = func(x_initial);
while(norm(grad_output) > epsilon && num_iter<N)
    num_iter=num_iter+1;
    grad_output = grad_f(x_output,func);
    p = - Hessian_f(x_output,func) \ grad_output;
    alpha = backtrack(grad_output, x_output, p,func);
    x_output = x_output + alpha*p;
    x_vals(num_iter+1,:) = x_output;
    f_vals(num_iter+1) = func(x_output);
end
f_output = func(x_output);
if(num_iter==N && (norm(grad_output))^2 > epsilon)
    fprintf('Maximum iterations reached but convergence did not happen.\n');
end

figure(3)
plot(0:num_iter,x_vals(1:num_iter+1,1),'r-o',0:num_iter,x_vals(1:num_iter+1,2),'b-o')
title("Newton's Method")
xlabel('Iteration Number')
ylabel('x values')
legend('x(1)','x(2)')

figure(4)
plot(0:num_iter,f_vals(1:num_iter+1,1),'r-o')
title("Newton's Method")
xlabel('Iteration Number')
ylabel('f values')


    function delF=grad_f(xstar,func) % write your logic for gradient computing function below
        delF = [0;0];
        h = 0.001;
        delF(1,1) = (func(xstar+[h;0])-func(xstar-[h;0]))/(2*h);
        delF(2,1) = (func(xstar+[0;h])-func(xstar-[0;h]))/(2*h);
    end
        
    function del2F=Hessian_f(xstar,func) % wrie your logic for Hessian computing function below
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
    
    function alpha=backtrack(delF, xk, pk,func)
        alpha0 = 5;
        rho = 0.8;
        c = 0.1;
        alpha = alpha0;
        fxk = func(xk);
        while(func(xk + alpha*pk) > fxk + c*alpha*(delF'*pk))
            alpha = alpha*rho;
        end
    end

 %   # End your code here

end

function [ x_output, f_output, grad_output]=quasi_newton_method(func, x_initial)
%-----------------------------------------------------------------------------------------------------------------------------
% Write your logic for Quasi newton method using in-exact line search. 

%Input parameters: 
            %func : input function to be evaluated
            % x_initial: initial value of x, a column vector 

%Returns:
%         x_output : converged x value, a column vector 
%         f_output : value of f at x_output
%         grad_output : value of gradient at x_output, a column vector
 %  -----------------------------------------------------------------------------------------------------------------------------

 %   # Start your code here
C = eye(2);
epsilon = 10e-6;
N = 15000;
x_output = x_initial;
grad_output = grad_f(x_output,func);
grad_next = grad_f(x_output,func);
num_iter = 0;
x_vals = zeros(N,2);
x_vals(1,:) = x_initial;
f_vals = zeros(N+1);
f_vals(1) = func(x_initial);
while(norm(grad_output) > epsilon && num_iter<N)
    num_iter=num_iter+1;
    grad_output = grad_next;
    p = - C * grad_output;
    alpha = backtrack(grad_output, x_output, p,func);
    s = alpha * p;
    x_output = x_output + s;
    grad_next = grad_f(x_output,func);
    y = grad_next - grad_output;
    C = (eye(2) - (s*y')/(y'*s))*C*(eye(2) - (y*s')/(y'*s)) + (s*s')/(y'*s);
    x_vals(num_iter+1,:) = x_output;
    f_vals(num_iter+1) = func(x_output);
end
if(num_iter==N && (norm(grad_output))^2 > epsilon)
    fprintf('Maximum iterations reached but convergence did not happen.\n');
end
f_output = func(x_output);

figure(5)
plot(0:num_iter,x_vals(1:num_iter+1,1),'r-o',0:num_iter,x_vals(1:num_iter+1,2),'b-o')
title("Quasi Newton's Method")
xlabel('Iteration Number')
ylabel('x values')
legend('x(1)','x(2)')

figure(6)
plot(0:num_iter,f_vals(1:num_iter+1,1),'r-o')
title("Quasi Newton's Method")
xlabel('Iteration Number')
ylabel('f values')


    function delF=grad_f(xstar,func) % write your logic for gradient computing function below
        delF = [0;0];
        h = 0.001;
        delF(1,1) = (func(xstar+[h;0])-func(xstar-[h;0]))/(2*h);
        delF(2,1) = (func(xstar+[0;h])-func(xstar-[0;h]))/(2*h);
    end
        
    function del2F=Hessian_f(xstar,func) % wrie your logic for Hessian computing function below
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
    
    function alpha=backtrack(delF, xk, pk,func)
        alpha0 = 5;
        rho = 0.8;
        c = 0.1;
        alpha = alpha0;
        fxk = func(xk);
        while(func(xk + alpha*pk) > fxk + c*alpha*(delF'*pk))
            alpha = alpha*rho;
        end
    end

 %   # End your code here

end
