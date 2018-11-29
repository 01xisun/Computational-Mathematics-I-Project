clear all
clc

%% Stability for QR Decompositions
% Q is any orthogonal matrix
[Q F] = qr(rand(50));

% R is any upper triangular matrix
R = triu(rand(50));

% Set matrix A by Q and R
A = Q*R;

% Condition number
kappa = cond(A);
%disp(sprintf('The condition number for A is %d',kappa))

% Construct the Q, R, in variety of ways%
[Q1 R1] = cgs_qr(A);
[Q2 R2] = mgs_qr(A);
[Q3 R3] = householder_qr(A);

% Verifying the forward error
% CGS
normq1 = norm(Q1-Q);
normr1 = norm(R1-R);
%disp(sprintf('The forward error for Q in CGS is %d',normq1))
%disp(sprintf('The forward error for R in CGS is %d.\n',normr1))
% MGS
normq2 = norm(Q2-Q);
normr2 = norm(R2-R);
%disp(sprintf('The forward error for Q in MGS is %d',normq2))
%disp(sprintf('The forward error for R in MGS is %d.\n',normr2))
% HH
normq3 = norm(Q3-Q);
normr3 = norm(R3-R);
%disp(sprintf('The forward error for Q in HH is %d',normq3))
%disp(sprintf('The forward error for R in HH is %d.\n',normr3))


% Verifying the backward error
% CGS
berrcgs = norm(A-Q1*R1)/norm(A);
disp(sprintf('The backward error for CGS is %d.\n',berrcgs))
% MGS
berrmgs = norm(A-Q2*R2)/norm(A);
disp(sprintf('The backward error for MGS is %d.\n',berrmgs))
% HH
berrhh = norm(A-Q3*R3)/norm(A);
disp(sprintf('The backward error for HH is %d.\n',berrhh))

% Constructiong the matrix with small perturbation
QQ1 = Q1+1e-15*rand(50);
QQ2 = Q2+1e-15*rand(50);
QQ3 = Q3+1e-15*rand(50);
RR1 = R1+1e-15*rand(50);
RR2 = R1+1e-15*rand(50);
RR3 = R1+1e-15*rand(50);

% Verifying the backward error again for accuracy
% CGS
berrcgs2 = norm(A-QQ1*RR1)/norm(A);
%disp(sprintf('The backward error for CGS with small perturbation is %d',berrcgs2))
% MGS
berrmgs2 = norm(A-QQ2*RR2)/norm(A);
%disp(sprintf('The backward error for MGS with small perturbation is %d',berrmgs2))
% HH
berrhh2 = norm(A-QQ3*RR3)/norm(A);
%disp(sprintf('The backward error for HH with small perturbation is %d',berrhh2))

%% Stability for SVD
% U is any orthogonal matrix
[U F] = qr(rand(50));

% V is another orthogonal matrix
[V F] = qr(rand(50));

% Sigma is any diagonal matrix
S = diag(sort(diag(rand(50)),'descend'));

% Constructing B
B = U*S*V';

% Construct U, Sigma, V
[U1 S1 V1] = svd(B);

% Verifying the forward error
normu1 = norm(U1-U);
norms1 = norm(S1-S);
normv1 = norm(V1-V);
%disp(sprintf('The forward error for U in SVD is %d',normu1))
%disp(sprintf('The forward error for Sigma in SVD is %d',norms1))
%disp(sprintf('The forward error for V in SVD is %d.\n',normv1))

% Verifying the backward error
berrsvd = norm(B-U1*S1*V1')/norm(B);
disp(sprintf('The backward error for SVD is %d.\n',berrsvd))

% Constructing matrix with small purterbations
UU1 = U1+1e-15*rand(50);
SS1 = S1+1e-15*rand(50);
VV1 = V1+1e-15*rand(50);

% Verifying the stability again
berrsvd2 = norm(B-UU1*SS1*VV1')/norm(B);
%disp(sprintf('The backward error for SVD with small perturbation is %d',berrsvd2))
