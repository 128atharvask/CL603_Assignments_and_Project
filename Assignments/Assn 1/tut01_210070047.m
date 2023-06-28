% This is the script file. Any function that you create should be defined at the end of the file.
% IMPORTANT NOTE: Don't write,clc,clear all,close all commands.



n=60; % number of values of x1 that will be created.
% Define relevant variable as following variable name:
% X1 = column vector of length n (i.e. n x 1) corresponding to equally spaced points in the interval [-1,1]. Note that first element of X1 will be -1 and the last element will be +1. Store X1 as a column vector.
X1 = transpose(linspace(-1,1,n));

a=X1'*X1; % This is being computed for evaluation purposes.
% Evaluate the test function at each of the n values in X1 with x2 fixed at -1. Store the function values in column vector fX_1 (i.e. it is of size n x 1). To evaluate the test
% function, create a matlab function named "givenfun" at the end of the script file and call it here. See instructions at the end of this file.
fX_1 = givenfun(X1,-1);



b=fX_1'*fX_1; % This is being computed for evaluation purposes. 

% Now Evaluate the function at each of the n values in X1 but with x2 fixed at +1. Store the function values in column vector fX_2 (i.e. it is of size n x 1). Once again call
% the created "givenfun" for this purpose. 
fX_2 = givenfun(X1,1);



c=fX_2'*fX_2; % This is being computed for evaluation purposes. 
%% -------------------Generate the two Line Plots as asked in the question ------------------
 
figure(1)
plot(X1,fX_1,'g-o',X1,fX_2,'r-x')
title('Test function versus X_1 for fixed X_2')
legend('X_2=-1','X_2=1')
xlabel('X_1')
ylabel('Test function')




%% -------------------Generate the Surface plot----------------
fX_3 = zeros(n,n);
%size(X1)
for i = 1:1:n
    fX_3(:,i) = givenfun(X1(i,1),X1);
end

figure(2)
surf(X1,X1,fX_3)
title('Surface Plot of test function')
xlabel('X_1')
ylabel('X_2')
zlabel('Test function')




%% ----------------- Generate the Contour plot ----------------
figure(3)
contour(X1,X1,fX_3,50,'ShowText','on')
title('Contour Plot of test function')
xlabel('X_1')
ylabel('X_2')



% ---------------- create function -----------------
% NOTE : User defined function in script file must be written at the end of code.
% Create function below with function name "givenfun",output argument "fx",input argument "x" where x=[x1 x2]^T. The test function given in the problem should be coded in this "givenfun".
% Thus "givenfun" should take a vector x as input (with components x1, x2) and give the value of test function as the output. 

function fx = givenfun(x1,x2)
    fx = (x2-x1).^4 + 8*x1.*x2 - x1 + x2 + 3;
end