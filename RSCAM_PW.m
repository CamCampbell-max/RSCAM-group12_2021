%% Simulation parameters 

% number of days of simulation
no_days = 700;

% dimension of grid 128x128, n denote number of internal points
n = 62;

make_gif = 0;

%% Main Calculations

% Define variables - figure (b)

eta = 0.03;
theta = 0.56;
Gamma = 0.0019;
% number of adjacent nodes
z = 4;
% note, dx = l
dx = 1;
dt = 0.01;
w = 1/15;
epsilon = theta * dt;

% D = l^2 / dt^2
D = dx^2 / dt;
% redefine gamma
gamma = Gamma / dx^2;
mu = (eta * D * dt) / (2 * z * dx^2);
alpha = epsilon * D * dt;
beta = (D * dt) / (2 * z * dx^2);

B_hat = epsilon * D * gamma / w;
A_hat = B_hat + 1/30;
%n_hat = 41/(n^2);
n_hat = 1;

% n_hat = Gamma * dt / (1 - exp(- A_hat * dt));

% boundary condition constant 
bc_B = 1/40;
bc_A = 1/40;
bc_rho = n_hat;

% rho - initial rho matrix !!! needs to be done
rho = ones([n+2,n+2]) * n_hat;

% initial A
% A needs to be full i.e include BCs 
A0 = ones([n+2,n+2])*1/30;
B0 = ones([n+2,n+2])*B_hat;

% slight perturbation
perturb_B = zeros([n+2,n+2]);

% perturb_B(16,16) = +0.06;
% perturb_B(30,28) = +0.06;
% perturb_B(5,9) = +0.06;
% perturb_B(55,49) = +0.06;
% perturb_B(25,22) = +0.06;
% perturb_B(50,14) = +0.06;
% perturb_B(50,14) = +0.06;


% set these values so loop works
A = A0;
B = B0 + perturb_B;

% 1D A0 for later
A0_alt = reshape(A0(2:end-1, 2:end-1), [n^2,1]);

%%% C1 constant for all time

% constructing C1 - requires D1 and D2
% D1
main_diag = ones([1,n]) * (1 + 4*mu);
off_diag = ones([1,n-1]) * -mu;

% construct D1
D1 = diag(main_diag) + diag(off_diag,1) + diag(off_diag, -1);
D1 = sparse(D1);

% D2
main_diag = ones([1,n]) * -mu;
D2 = diag(main_diag);
D2 = sparse(D2);

% Construct block matrix C1
C1 = blktridiag(D1, D2, D2, n);

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

% G1 constant for all time

% H1
main_diag = ones([1,n]) * (1 + 4*beta);
off_diag = ones([1,n-1]) * -beta;

H1 = diag(main_diag) + diag(off_diag,1) + diag(off_diag, -1);
H1 = sparse(H1);

% H2
main_diag = ones([1,n]) * -beta;
H2 = diag(main_diag);

% construct G1
G1 = blktridiag(H1, H2, H2, n);

% adjusting for boundary conditions
% B_bc_adj = zeros([n,n]);
% B_bc_adj(1:end,1) = 2 * mu * bc_B;
% B_bc_adj(1:end,end) = 2*mu*bc_B;
% B_bc_adj(1,1:end) = 2*mu*bc_B;
% B_bc_adj(end,1:end) = 2*mu*bc_B;
% B_bc_adj(1,1) = 4*mu *bc_B;
% B_bc_adj(1,end) = 4*mu*bc_B;
% B_bc_adj(end,1) = 4*mu*bc_B;
% B_bc_adj(end,end) = 4*mu*bc_B;
% B_bc_adj = reshape(B_bc_adj, [n^2,1]);
% 
% % adjusting for boundary conditions
% rho_bc_adj = zeros([n,n]);
% rho_bc_adj(1:end,1) = 2*beta*bc_rho;
% rho_bc_adj(1:end,end) = 2*beta*bc_rho;
% rho_bc_adj(1,1:end) = 2*beta*bc_rho;
% rho_bc_adj(end,1:end) = 2*beta*bc_rho;
% rho_bc_adj(1,1) = 4*beta*bc_rho;
% rho_bc_adj(1,end) = 4*beta*bc_rho;
% rho_bc_adj(end,1) = 4*beta*bc_rho;
% rho_bc_adj(end,end) = 4*beta*bc_rho;
% rho_bc_adj = reshape(rho_bc_adj, [n^2,1]);

%intialise cells to store matrices at each step
B_list = cell([1, no_days]);
A_list = cell([1, no_days]);
rho_list = cell([1, no_days]);

B_list{1} = B;
A_list{1} = A;
rho_list{1} = rho;

for day = 2:no_days
    
    % reshape
    B = reshape(B(2:end-1,2:end-1), [n^2,1]);
    rho = reshape(rho(2:end-1,2:end-1), [n^2,1]);
    
    % compute f_matrix
    F = f_A_2(A, beta, dt);

    % D4
    % main_diag = ones([1,n]) * mu;
    % D4 = diag(main_diag);
    % D4 = sparse(D4);

    % construct C2

    % alternate approach

    % reshape rho 
    rho_1D = reshape(rho, [1, n^2]);
    
    main_diag = ones([1,n^2]) * (1 - 4*mu - w*dt) + rho_1D * alpha;
    %main_diag = ones([1,n^2]) * (1 - 4*mu - w*dt);
    
    off_diag = ones([1,n^2 - 1]) * mu;
    far_off_diag = ones([1,n^2 - n]) * mu;

    % concatenate each
    i = cat(2, main_ind_i, sub_ind_i, sup_ind_i, far_sub_ind_i, far_sup_ind_i);
    j = cat(2, main_ind_j, sub_ind_j, sup_ind_j, far_sub_ind_j, far_sup_ind_j);
    values = cat(2, main_diag, off_diag, off_diag, far_off_diag, far_off_diag);

    % construct the block matrix
    C2 = sparse(i, j, values);

    % now the rho discretisation

    % G2 0 requires H3 and H4

    % reshape f(A)
    F_1D = reshape(F, [1, n^2]);
    main_diag = ones([1,n^2]) * (1 - 4*beta) + F_1D;

    % NOTE: these are solely the iternal nodes of A1, A2, A3, A4, may have to
    % use periodicity in order to get boundary conditions?

    % create A matrices (previously A1)
    A1 = zeros([n,n]);
    A2 = zeros([n,n]);
    A3 = zeros([n,n]);
    A4 = zeros([n,n]);
    for k = 2:(n+1)
        for l = 2:(n+1)
            A1(k-1,l-1) = 1 - 1/A(k,l) * (A(k,l+1) - A(k,l-1));
            A2(k-1,l-1) = 1 + 1/A(k,l) * (A(k,l+1) - A(k,l-1));
            A3(k-1,l-1) = 1 - 1/A(k,l) * (A(k+1,l) - A(k-1,l));
            A4(k-1,l-1) = 1 + 1/A(k,l) * (A(k+1,l) - A(k-1,l));
        end
    end

%     %create A2
%     A2 = zeros([n,n]);
%     for k = 2:n+1
%         for l = 2:n+1
%             A2(k-1,l-1) = 1 + 1/A(k,l) * (A(k,l+1) - A(k,l-1));
%         end
%     end
% 
%     % create A3
%     A3 = zeros([n,n]);
%     for k = 2:n+1
%         for l =2:n+1
%             A3(k-1,l-1) = 1 - 1/A(k,l) * (A(k+1,l) - A(k-1,l));
%         end
%     end
% 
%     % create A4
%     A4 = zeros([n,n]);
%     for k = 2:n+1
%         for l =2:n+1
%             A4(k-1,l-1) = 1 + 1/A(k,l) * (A(k+1,l) - A(k-1,l));
%         end
%     end

    % sub diag
    A2_new = reshape(A2, [1, n^2]);
    sub_diag = A2_new(2:end) * beta;

    % sup diag
    A1_new = reshape(A1, [1, n^2]);
    sup_diag = A1_new(1:end-1) * beta;

    % far sub diag
    far_sub_diag = beta * reshape(A4(2:n,1:end), [1,n^2 - n]);

    % far sup diag
    far_sup_diag = beta * reshape(A3(1:n-1,1:end), [1,n^2 - n]);

    % assemble values as before, can use same indexing
    values = cat(2, main_diag, sub_diag, sup_diag, far_sub_diag, far_sup_diag);

    % construct G2
    G2 = sparse(i, j, values);

    % can use this to plot sparsity - spy(G2)

    % SOLVING

    rho_vec = reshape(rho, [n^2,1]);

    % solving b system
    % B_next = C1 \ (C2*B + alpha*(rho_vec.*A0_alt) + B_bc_adj); boundary
    % cond
    B_next = C1 \ (C2*B + alpha*(rho_vec.*A0_alt));

    % create vector
    gamma_vec = ones([n^2,1]) * gamma * dt;

    % solve rho system
    % rho_next = G1 \ (G2*rho_vec + gamma_vec + rho_bc_adj);
    rho_next = G1 \ (G2*rho_vec + gamma_vec);

    B_mesh = reshape(B_next, [n,n]);
    rho_mesh = reshape(rho_next, [n,n]);

    % apply boundary conditions
    B = zeros([n+2, n+2]);
    B(2:end-1, 2:end-1) = B_mesh;
    B(:,1) = bc_B;
    B(:,end) = bc_B;
    B(1,:) = bc_B;
    B(end,:) = bc_B;
    
    % apply boundary conditions
    rho = zeros([n+2, n+2]);
    rho(2:end-1, 2:end-1) = rho_mesh;
    rho(:,1) = bc_rho;
    rho(:,end) = bc_rho;
    rho(1,:) = bc_rho;
    rho(end,:) = bc_rho;
    
    % Calculate A
    A = A0 + B;
    
    % apply boundary conditions
    A(:,1) = bc_A;
    A(:,end) = bc_A;
    A(1,:) = bc_A;
    A(end,:) = bc_A;
    
    % Store results in cells
    A_list{day} = A;
    B_list{day} = B;
    rho_list{day} = rho;

end

%% Plotting, Animation, and Output

% suppress figures
set(0,'DefaultFigureVisible','off');

% animate B
loops = no_days;
M(loops) = struct('cdata',[],'colormap',[]);

if make_gif == 1

    for frame = 1:no_days

        hfig=figure;

        contourf(cell2mat(A_list(frame)));
        colorbar;
        
        %caxis([0.02, 0.3]);

        M(frame) = getframe(hfig);
        im{frame} = frame2im(M(frame));
    end

    % movie(M)



    %% Saving gif

    filename = 'PW_method.gif'; % Specify the output file name
    for idx = 1:no_days
        [A,map] = rgb2ind(im{idx},256);
        if idx == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.02);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.02);
        end
    end

    % suppress figures
    set(0,'DefaultFigureVisible','on');
    
else
    hfig = figure;
    contourf(cell2mat(A_list(no_days)));
    colorbar;
    drawnow;
   
    exportgraphics(hfig,'PW_Output.png','Resolution',300);
end