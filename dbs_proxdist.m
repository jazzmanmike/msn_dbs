function D=dbs_proxdist(A,k)
% This function calculates the distance matrix
% from a proximity matrix A

%D=((1./A)-1)^k;

[n,m]=size(A);

for i=1:n
    for j=1:m
        D(i,j)=((1/A(i,j))-1)^k;
    end
end