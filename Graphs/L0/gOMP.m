function [A]=gOMP(D,X,L,S)
%=============================================
% Sparse coding of a group of signals based on a given 
% dictionary and specified number of atoms to use. 
% input arguments: 
%       D - the dictionary (its columns MUST be normalized).
%       X - the signals to represent
%       L - the max. number of coefficients for each signal.
% output arguments: 
%       A - sparse coefficient matrix.
%=============================================

if nargin < 4
    S = round(max(L/4, 1));
end
[n,P]=size(X);
[n,K]=size(D);
for k=1:1:P
    a=[];
    x=X(:,k);
    residual=x; 
    indx0=[];
    for j=1:1:L
        proj=D'*residual;
        [~,pos]=sort(abs(proj),'descend');
        indx = union(indx0,pos(1:S));
        a=pinv(D(:,indx))*x;
        residual=x-D(:,indx)*a;
        indx0=indx;
        if sum(residual.^2)<1e-6,  break;   end
    end
    temp=zeros(K,1);
    temp(indx0)=a;
    A(:,k)=sparse(temp);
end
return;