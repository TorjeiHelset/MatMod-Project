N = 50;
a = 1:N-1;
b = 2:N;

c = transpose(ones(N,1));
A = c'*a;
A

d = 0:N-1;
e = N.*d;
A = A + e'*transpose(ones(N-1,1));

A_ = transpose(A)
index1 = reshape(A_,1,[]);
index2 = reshape(A_+1,1,[]);
length(index1)