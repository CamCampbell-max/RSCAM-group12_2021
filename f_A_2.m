function F = f_A_2(A, beta, dt)
% Compute the matrix f_A return the internal nodes

% get dimension of the matrix
[n,n] = size(A)

% init empty matrix
M = zeros(n);

% i indices
main_ind_i = 1:n^2;
sub_ind_i = 2:n^2;
sup_ind_i = 1:(n^2 - 1);
far_sub_ind_i = n+1:n^2;
far_sup_ind_i = 1:(n^2 - n);

% j indices
main_ind_j = 1:n^2;
sup_ind_j = 2:n^2;
sub_ind_j = 1:(n^2 - 1);
far_sub_ind_j = 1:(n^2 - n);
far_sup_ind_j = n+1:n^2;

i = cat(2, main_ind_i, sub_ind_i, sup_ind_i, far_sub_ind_i, far_sup_ind_i);
j = cat(2, main_ind_j, sub_ind_j, sup_ind_j, far_sub_ind_j, far_sup_ind_j);

%%
main_diag = ones([1,n^2]) * -16 * beta
off_diag = ones([1,n^2 - 1]) * 4*beta;
far_off_diag = ones([1,n^2 - n]) * 4*beta;
values = cat(2, main_diag, off_diag, off_diag, far_off_diag, far_off_diag);

% reshape A into vector
A_long = reshape(A, [n^2,1]);

% first sparse matrix
mat_1 = sparse(i, j, values);

% get first matrix
A1 = mat_1 * A_long
A1 = (A1 .* (A_long.^(-1)))

%%

main_diag = zeros([1,n^2])
off_diag = ones([1,n^2 - 1]);
far_off_diag = zeros([1,n^2 - n]);
values = cat(2, main_diag, -off_diag, off_diag, far_off_diag, far_off_diag);

% second sparse matrix
mat_2 = sparse(i, j, values);

A2 = mat_2 * A_long
A2 = A2.^2

%%

main_diag = zeros([1,n^2])
off_diag = zeros([1,n^2 - 1]);
far_off_diag = ones([1,n^2 - n]);
values = cat(2, main_diag, off_diag, off_diag, -far_off_diag, far_off_diag);

% second sparse matrix
mat_3 = sparse(i, j, values);

A3 = mat_3 * A_long
A3 = A3.^2

A2 = - beta * (A_long .* (A2 + A3))

final = (A1 + A2 + dt*A_long)
final = reshape(final, [n,n])
F = final(2:end-1, 2:end-1)
end

%%

% % compute the interior nodes
% for i = 2:n-1
%     for j = 2:n-1
%         M(i,j) = 4 * beta / (A(i,j)) * (A(i+1,j) + A(i-1,j) -4*A(i,j) + A(i,j+1) + A(i,j-1))
%         M(i,j) = M(i,j) - beta / ((A(i,j)^2)) * ((A(i+1,j) + A(i-1,j))^2 + (A(i,j+1) + A(i,j-1))^2)
%         M(i,j) = M(i,j) + dt * A(i,j)
%     end
% end
% 
% % match to initial condition
% M(1:end, 1) = A(1:end, 1)
% M(1:end, end) = A(1:end, end)
% M(1, 1:end) = A(1, 1:end)
% M(end, 1:end) = A(end, 1:end)

% return value
%F = M;
%end