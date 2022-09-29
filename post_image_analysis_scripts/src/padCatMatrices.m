function C = padCatMatrices(A, B)

%Concatenates two matrices (M x N and M x O) into
% a matrix of size M x O, where O >=N. 
%values are padded with NaNs. Zeros in the 
%original data are temporarily stored as -Inf. 

%we'll save these for testing later
A_old = A;
B_old = B;

A(A == 0) = -Inf;
B(B == 0) = -Inf;

[sA1, sA2] = size(A);
[sB1, sB2] = size(B);
C(sA1+1:sA1+sB1, 1:sB2) = B;
C(1:sA1, 1:sA2)         = A;

C(C==0) = NaN;
C(C == -Inf) = 0;

%Check whether C has the desired dimensions
% ( M x O )
assert(size(C, 1) == size(A_old, 1) + size(B_old, 1));
assert(size(C, 2) == max(size(A_old, 2), size(B_old, 2)));

%make sure the chunks in this matrix correspond 
%to the original matrices
assert( isequaln( C(1:sA1, 1:sA2), A_old) );
assert( isequaln( C(sA1+1:sA1+sB1, 1:sB2), B_old) );


end