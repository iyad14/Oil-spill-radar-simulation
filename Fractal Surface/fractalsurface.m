function [row,col,R] = fractalsurface(n,H)
% inputs:
% n - number of iterations
% H - Hurst exponent, fractal dimension D = 3-H
% outputs:
% row - row location of surface
% col - column location of surface
% R - surface matrix
% Authors: Hui Yang
% Affiliations:
% The Pennsylvania State University
% 310 Leohard Building, University Park, PA
% Email: yanghui@gmail.com
% If you find this toolbox useful, please cite the following papers:
% [1]	Y. Chen and H. Yang, Numerical simulation and pattern characterization of
% nonlinear spatiotemporal dynamics on fractal surfaces for the whole-heart modeling
% applications, European Physical Journal, DOI: 10.1140/epjb/e2016-60960-6
% [2]	B. Yao, F. Imani, A. Sakpal, E. W. Reutzel, and H. Yang*, �Multifractal 
% analysis of image profiles for the characterization and detection of defects 
% in additive manufacturing,� ASME Journal of Manufacturing Science and Engineering, 
% Vol. 140, No. 3, p031014-13, 2017, DOI: 10.1115/1.4037891
% [3]	F. Imani, B. Yao, R. Chen, P. Rao, and H. Yang*, �Joint multifractal 
% and lacunarity analysis of image profiles for manufacturing quality control�. 
% ASME Journal Manufacturing Science and Engineering, Vol. 141, No. 4, 
% p 044501-7, 2019. DOI: 10.1115/1.4042579
N=2^n; % size of surface is (N+1)*(N+1)
L=zeros(N+1,N+1); % initial of surface
%% first iteration
R=randn(3,3);
Initial = [1,2^(n-1)+1,2^n+1];
for i=1:3
    for j=1:3
        L(Initial(i),Initial(j))=R(i,j);
    end
end
[row,col]=find(L); % location of first iteration
r=sqrt(length(row));
%% generation of the rest iterations
for l=1:n-1
    for i=1:r-1
        mid(i)=(row(i)+row(i+1))/2;
        for j=1:r
            L(mid(i),row(j))=mean([L(row(i),row(j)),L(row(i+1),row(j))])+randn*2^(-H*(l+1));
            L(row(j),mid(i))=mean([L(row(j),row(i)),L(row(j),row(i+1))])+randn*2^(-H*(l+1));
        end
    end
    for i=1:r-1
        for j=1:r-1
            L(mid(i),mid(j))=mean([L(row(i),row(j)),L(row(i),row(j+1)),L(row(i+1),row(j)),L(row(i+1),row(j+1))])+randn*2^(-H*(l+1));
        end
    end
    
    [row,col]=find(L);
    r=sqrt(length(row));
    
end
R=zeros(r,r);
for i=1:r
    for j=1:r
        R(i,j)=L(row(i),row(j));
    end
end