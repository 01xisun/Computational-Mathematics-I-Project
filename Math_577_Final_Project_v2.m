clear all
clc

%% Checking the accuracy
A = rand(50);
x = [];
for i=1:50
    x = [x 1];
end
x = x';
b = A*x;
% Classical Gram-Smidth
[Q1 R1] = cgs_qr(A);
x1 = R1\(Q1'*b);
% Modified Gram-Smidth
[Q2 R2] = mgs_qr(A);
x2 = R2\(Q2'*b);
% Householder
[Q3 R3] = householder_qr(A);
x3 = R3\(Q3'*b);
% SVD
[U S V] = svd(A,0);
x4 = V*(S\(U'*b));
% Put everything together
format long
X = [x1 x2 x3 x4];
Xerr = [norm(x1-x) norm(x2-x) norm(x3-x) norm(x4-x)];
disp(sprintf('The norm between the actual x and CGS is %d.\n',Xerr(1)))
disp(sprintf('The norm between the actual x and MGS is %d.\n',Xerr(2)))
disp(sprintf('The norm between the actual x and HH is %d.\n',Xerr(3)))
disp(sprintf('The norm between the actual x and SVD is %d.\n',Xerr(4)))




