clear all
clc
% Fix column=50, variate number of row
% CGS
normq1 = [];
normr1 = [];
backerr1 = [];
Number.of.Rows = [];
for m = 50:200
    [Q F] = qr(rand(m,50),0);
    R = triu(rand(50));
    A = Q*R;
    [Q1 R1] = cgs_qr(A);
    normq1 = [normq1 norm(Q1-Q)];
    normr1 = [normr1 norm(R1-R)];
    backerr1 = [backerr1 norm(A-Q1*R1)/norm(A)];
    Number.of.Rows = [Number.of.Rows m];
end
figure(14)
plot(Number.of.Rows,normq1,'o')
title('Forward Error for CGS, Q')
xlabel('Number of Rows')
ylabel('Error')
pause
figure(15)
plot(Number.of.Rows,normr1,'o')
title('Forward Error for CGS, R')
xlabel('Number of Rows')
ylabel('Error')
pause
figure(16)
plot(Number.of.Rows,backerr1,'o')
title('Backward Error for CGS')
xlabel('Number of Rows')
ylabel('Error')
pause

% MGS
normq2 = [];
normr2 = [];
backerr2 = [];
for m = 50:200
    [Q F] = qr(rand(m,50),0);
    R = triu(rand(50));
    A = Q*R;
    [Q2 R2] = mgs_qr(A);
    normq2 = [normq2 norm(Q2-Q)];
    normr2 = [normr2 norm(R2-R)];
    backerr2 = [backerr2 norm(A-Q2*R2)/norm(A)];
end
figure(17)
plot(Number.of.Rows,normq2,'o')
title('Forward Error for MGS, Q')
xlabel('Number of Rows')
ylabel('Error')
pause
figure(18)
plot(Number.of.Rows,normr2,'o')
title('Forward Error for MGS, R')
xlabel('Number of Rows')
ylabel('Error')
pause
figure(19)
plot(Number.of.Rows,backerr2,'o')
title('Backward Error for MGS')
xlabel('Number of Rows')
ylabel('Error')
pause

% Householder
normq3 = [];
normr3 = [];
backerr3 = [];
for m = 50:200
    [Q F] = qr(rand(m,50));
    R = triu(rand(m,50));
    A = Q*R;
    [Q3 R3] = householder_qr(A);
    normq3 = [normq3 norm(Q3-Q)];
    normr3 = [normr3 norm(R3-R)];
    backerr3 = [backerr3 norm(A-Q3*R3)/norm(A)];
end
figure(20)
plot(Number.of.Rows,normq3,'o')
title('Forward Error for Householder, Q')
xlabel('Number of Rows')
ylabel('Error')
pause
figure(21)
plot(Number.of.Rows,normr3,'o')
title('Forward Error for Householder, R')
xlabel('Number of Rows')
ylabel('Error')
pause
figure(22)
plot(Number.of.Rows,backerr3,'o')
title('Backward Error for Householder')
xlabel('Number of Rows')
ylabel('Error')
pause

% SVD
normu = [];
norms = [];
normv = [];
backerr4 = [];
for m = 50:200
    [U F] = qr(rand(m,50),0);
    S = diag(sort(diag(rand(50)),'descend'));
    [V F] = qr(rand(15));
    B = U*S*V';
    [U1 S1 V1] = svd(B,0);            
    normu = [normu norm(U1-U)];
    norms = [norms norm(S1-S)];
    normv = [normv norm(V1-V)];
    backerr4 = [backerr4 norm(B-U1*S1*V1')/norm(B)];
end
figure(23)
plot(Number.of.Rows,normu,'o')
title('Forward Error for SVD, U')
xlabel('Number of Rows')
ylabel('Error')
pause
figure(24)
plot(Number.of.Rows,norms,'o')
title('Forward Error for SVD, Sigma')
xlabel('Number of Rows')
ylabel('Error')
pause
figure(25)
plot(Number.of.Rows,normv,'o')
title('Forward Error for SVD, V')
xlabel('Number of Rows')
ylabel('Error')
pause
figure(26)
plot(Number.of.Rows,backerr4,'o')
title('Backward Error for SVD')
xlabel('Number of Rows')
ylabel('Error')






