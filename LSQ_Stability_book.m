% This script illustrates the stability of various
% least-squares algorithms
clear all
clc
m = 100; n = 15;
disp(sprintf('We use a degree %d polynomial to fit the function exp(sin(4*x))', n-1))
disp(sprintf('at %d equally spaced points in [0,1].\n', m))
%scrnsz = get(0, 'Screensize');
xi = linspace(0,1,2*m+1);
eta = exp(sin(4*xi));
%figure('Position', [5 10 scrnsz(3)-5 scrnsz(4)-60])
%figure
%PolynomialLSQ(m,n-1);
clc
t = linspace(0,1,m)';
A = [];
for i=1:n
    A = [A t.^(i-1)];
end
b = exp(sin(4*t));
% normalized so that x(15) below should be 1
b = b/2006.787453080206;
% Solve LSQ problem and compute condition number
x = A \ b;
y = A*x;
kappa = cond(A);
theta = asin(norm(b-y)/norm(b));
eta = norm(A)*norm(x)/norm(y);
disp(sprintf('cond(A) = %e', kappa))
disp(sprintf('A large condition number shows that the problem (matrix)'))
disp(sprintf('is ill-conditioned.\n'))
disp(sprintf('theta = %e', theta))
disp(sprintf('A small value of theta shows that it is possible'))
disp(sprintf('to fit the data well with a degree %d polynomial.\n', n-1))
disp(sprintf('eta = %e', eta))
disp(sprintf('eta can be anywhere between 1 and cond(A)\n'))
disp(sprintf('Condition of least squares problem'))
illcond = kappa + kappa^2*tan(theta)/eta;
disp(sprintf('kappa + kappa^2*tan(theta)/eta = %e', illcond))
loss = floor(log10(illcond));
disp(sprintf('We can expect to lose %d digits due to ill-conditioning', loss))
pause
clc
disp(sprintf('LSQ solution with degree %d polynomial.',n-1))
disp(sprintf('The problem is normalized so that the coefficient'))
disp(sprintf('computed below should be 1.\n'))
% Check classical Gram-Schmidt
[Q, R] = cgs_qr(A);
x = R \ (Q'*b);
disp(sprintf('via classical Gram-Schmidt: c(%d) = %16.14f', n, x(15)))
disp(sprintf('                    error: %e', abs(1-x(15))))
loss_err = -floor(log10(abs(1-x(15))));
if loss_err < 7    % lost more than 10 digits
    disp(sprintf('Unstable algorithm\n'))
else
    disp(sprintf('Stable algorithm\n'))
end
pause
% Check modified Gram-Schmidt
[Q, R] = mgs_qr(A);
x = R \ (Q'*b);
disp(sprintf('via modified Gram-Schmidt: c(%d) = %16.14f', n, x(15)))
disp(sprintf('                    error: %e', abs(1-x(15))))
loss_err = -floor(log10(abs(1-x(15))));
if loss_err < 7    % lost more than 10 digits
    disp(sprintf('Unstable algorithm\n'))
else
    disp(sprintf('Stable algorithm\n'))
end
pause
% Check Householder QR
[Q, R] = qr(A, 0);
x = R \ (Q'*b);
disp(sprintf('via Householder QR: c(%d) = %16.14f', n, x(15)))
disp(sprintf('                    error: %e', abs(1-x(15))))
loss_err = -floor(log10(abs(1-x(15))));
if loss_err < 7    % lost more than 10 digits
    disp(sprintf('Unstable algorithm\n'))
else
    disp(sprintf('Stable algorithm\n'))
end
pause
% Check SVD
[U, S, V] = svd(A,0);
x = V*(S \ (U'*b));
disp(sprintf('via SVD: c(%d) = %16.14f', n, x(15)))
disp(sprintf('                    error: %e', abs(1-x(15))))
loss_err = -floor(log10(abs(1-x(15))));
if loss_err < 7    % lost more than 10 digits
    disp(sprintf('Unstable algorithm.\n'))
else
    disp(sprintf('Stable algorithm.\n'))
end