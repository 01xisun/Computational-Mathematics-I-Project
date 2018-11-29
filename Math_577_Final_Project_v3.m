clear all
clc
% Fix row=500, variate number of column
% CGS
normq1 = [];
normr1 = [];
backerr1 = [];
for n = 1:200
    [Q F] = qr(rand(200,n),0);
    R = triu(rand(n));
    A = Q*R;
    [Q1 R1] = cgs_qr(A);
    normq1 = [normq1 norm(Q1-Q)];
    normr1 = [normr1 norm(R1-R)];
    backerr1 = [backerr1 norm(A-Q1*R1)/norm(A)];
end
figure(1)
plot(normq1,'o')
title('Forward Error for CGS, Q')
xlabel('Number of Coulmns')
ylabel('Error')
pause
figure(2)
plot(normr1,'o')
title('Forward Error for CGS, R')
xlabel('Number of Coulmns')
ylabel('Error')
pause
figure(3)
plot(backerr1,'o')
title('Backward Error for CGS')
xlabel('Number of Coulmns')
ylabel('Error')
pause

% MGS
normq2 = [];
normr2 = [];
backerr2 = [];
for n = 1:200
    [Q F] = qr(rand(200,n),0);
    R = triu(rand(n));
    A = Q*R;
    [Q2 R2] = mgs_qr(A);
    normq2 = [normq2 norm(Q2-Q)];
    normr2 = [normr2 norm(R2-R)];
    backerr2 = [backerr2 norm(A-Q2*R2)/norm(A)];
end
figure(4)
plot(normq2,'o')
title('Forward Error for MGS, Q')
xlabel('Number of Coulmns')
ylabel('Error')
pause
figure(5)
plot(normr2,'o')
title('Forward Error for MGS, R')
xlabel('Number of Coulmns')
ylabel('Error')
pause
figure(6)
plot(backerr2,'o')
title('Backward Error for MGS')
xlabel('Number of Coulmns')
ylabel('Error')
pause

% Householder
normq3 = [];
normr3 = [];
backerr3 = [];
for n = 1:200
    [Q F] = qr(rand(200,n));
    R = triu(rand(200,n));
    A = Q*R;
    [Q3 R3] = householder_qr(A);
    normq3 = [normq3 norm(Q3-Q)];
    normr3 = [normr3 norm(R3-R)];
    backerr3 = [backerr3 norm(A-Q3*R3)/norm(A)];
end
figure(7)
plot(normq3,'o')
title('Forward Error for Householder, Q')
xlabel('Number of Coulmns')
ylabel('Error')
pause
figure(8)
plot(normr3,'o')
title('Forward Error for Householder, R')
xlabel('Number of Coulmns')
ylabel('Error')
pause
figure(9)
plot(backerr3,'o')
title('Backward Error for Householder')
xlabel('Number of Coulmns')
ylabel('Error')
pause

% SVD
normu = [];
norms = [];
normv = [];
backerr4 = [];
for n = 1:200
    [U F] = qr(rand(200,n),0);
    S = diag(sort(diag(rand(n)),'descend'));
    [V F] = qr(rand(n));
    B = U*S*V';
    [U1 S1 V1] = svd(B,0);            
    normu = [normu norm(U1-U)];
    norms = [norms norm(S1-S)];
    normv = [normv norm(V1-V)];
    backerr4 = [backerr4 norm(B-U1*S1*V1')/norm(B)];
end
figure(10)
plot(normu,'o')
title('Forward Error for SVD, U')
xlabel('Number of Coulmns')
ylabel('Error')
pause
figure(11)
plot(norms,'o')
title('Forward Error for SVD, Sigma')
xlabel('Number of Coulmns')
ylabel('Error')
pause
figure(12)
plot(normv,'o')
title('Forward Error for SVD, V')
xlabel('Number of Coulmns')
ylabel('Error')
pause
figure(13)
plot(backerr4,'o')
title('Backward Error for SVD')
xlabel('Number of Coulmns')
ylabel('Error')

