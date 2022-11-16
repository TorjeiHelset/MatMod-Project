N = 4;
M = 2;

a = 1:N-1;
c = transpose(ones(M,1));
A = c'*a;
d = 0:M-1;
e = N.*d;
A = A + e'*transpose(ones(N-1,1));
A_ = transpose(A);
index1 = reshape(A_,1,[]);
index1