function K_R_product = Khatri_Rao(A_mat, B_mat)
% Generating a Khatri-Rao product, that is to say a column-wise Kronecker product

[row_A,col_A] = size(A_mat);
[row_B,col_B] = size(B_mat);

if col_A == col_B
    col_n = col_A;
else
    error('Dimensions of input matrices must agree, such as numbers of column for the input matrices');
end
K_R_product = zeros(row_A*row_B,col_n);
for ii = 1:col_n
    K_R_product(:,ii) = kron(A_mat(:,ii),B_mat(:,ii));
end