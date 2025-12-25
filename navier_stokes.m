
clc; clear; close all;

% ParamÃ¨tres
rho = 2052; mu0 = 0.0056; Ri = 0; Ro = 0.0125; L = 0.02; U0 = 0.1; 
mu_inf = 0.0035; lambda = 1.9; n = 0.3526;  n_iter = 10; a = 2.0; 

%% 2) Maillage P2-P1 avec PDE Toolbox
n_r = 8;   % Subdivisions radiales
n_z = 4;  % Subdivisions axiales
model = createpde(2);
R1 = [3,4, [0, Ro, Ro, 0], [0, 0, L, L]]';
g = decsg(R1);
geometryFromEdges(model, g);
mesh = generateMesh(model, 'Hgrad', 1.2, 'Hmax', Ro/n_r, 'GeometricOrder', 'quadratic');

nodes_P2 = mesh.Nodes';          
tri_P2 = mesh.Elements(1:6,:)';  
numT = size(tri_P2, 1);          
tri_P1 = mesh.Elements(1:3,:)';  
all_P1_nodes = unique(tri_P1(:));
nodes_P1 = nodes_P2(all_P1_nodes, :);
numP1 = size(nodes_P1, 1);
numP2 = size(nodes_P2, 1);

fprintf('Nombre d''Ã©lÃ©ments : %d\n', numT);
fprintf('Nombre de nÅ“uds P1 : %d\n', numP1);
fprintf('Nombre de nÅ“uds P2 : %d\n', numP2);

%% 3) Fonctions de forme de rÃ©fÃ©rence sur le triangle iso (P1 & P2)
psi1 = @(x,y) 1 - x - y;     
psi2 = @(x,y) x;             
psi3 = @(x,y) y;             
piG = @(x,y) [psi1(x,y).*(2*psi1(x,y)-1); psi2(x,y).*(2*psi2(x,y)-1); ...
    psi3(x,y).*(2*psi3(x,y)-1); 4*psi1(x,y).*psi2(x,y); ...
    4*psi2(x,y).*psi3(x,y); 4*psi3(x,y).*psi1(x,y)];

gradP2 = @(x,y) [ ...
    [4*x + 4*y - 3, 4*x + 4*y - 3]; 
    [4*x - 1, 0];                   
    [0, 4*y - 1];                   
    [4 - 8*x - 4*y, -4*x];          
    [4*y, 4*x];                     
    [-4*y, 4 - 4*x - 8*y]];         

%% 4) Initialisation des solutions
u_r = zeros(numP2,1);   
u_z = zeros(numP2,1);   
%u_z = U0 * (1 - (nodes_P2(:,1)/Ro).^2);  % Profil parabolique
p = zeros(numP1,1);     
U = [u_r; u_z; p]; % Vecteur initial

%% 5) Quadrature de Gauss Ã  7 points
quad_pts = [...
    1/3, 1/3; 
    (6-sqrt(15))/21, (6-sqrt(15))/21;
    (6-sqrt(15))/21, (9+2*sqrt(15))/21;
    (9+2*sqrt(15))/21, (6-sqrt(15))/21;
    (6+sqrt(15))/21, (6+sqrt(15))/21;
    (6+sqrt(15))/21, (9-2*sqrt(15))/21;
    (9-2*sqrt(15))/21, (6+sqrt(15))/21];
quad_wts = [...
    9/80; 
    (155 - sqrt(15))/2400; 
    (155 - sqrt(15))/2400; 
    (155 - sqrt(15))/2400; 
    (155 + sqrt(15))/2400; 
    (155 + sqrt(15))/2400; 
    (155 + sqrt(15))/2400];

%% 6) Initialisation des matrices globales
K = sparse(2*numP2, 2*numP2);   
M = sparse(2*numP2, 2*numP2); 
B = sparse(numP1, 2*numP2);          
f = zeros(2*numP2, 1);    
rr = zeros(numT, 7);

%% 7) Assemblage des matrices Ã©lÃ©mentaires

for e = 1:numT
    nodes_e_P2 = tri_P2(e, :);
    nodes_e_P1 = tri_P1(e, :);
    coords_P2 = nodes_P2(nodes_e_P2, :);
    r = coords_P2(:, 1); z = coords_P2(:, 2);
    
    Bmat = [r(2)-r(1), r(3)-r(1); z(2)-z(1), z(3)-z(1)];
    detB = abs(det(Bmat));
    invB = inv(Bmat);

    near_axis = any(r < 1e-5 * Ro);
    
    M_e = zeros(12, 12); 
    K_e = zeros(12, 12);  
    B_e = zeros(3, 12);   
    f_e = zeros(12, 1);

    for q = 1:7
        xi = quad_pts(q, 1); 
        eta = quad_pts(q, 2);

        % Quadrature standard
        rq = sum(r.*piG(xi, eta));
        rr(e, q) = rq;
        zq = z(1) + (z(2)-z(1))*xi + (z(3)-z(1))*eta;
        wq = quad_wts(q) * detB;
        phi = piG(xi, eta);
        dphi_ref = gradP2(xi, eta);
        dphi = (invB' * dphi_ref')';
        psi = [psi1(xi, eta); psi2(xi, eta); psi3(xi, eta)];
        
        % Assemblage 
        fact = 2 * pi * rq * wq;
        M_e(1:6, 1:6) = M_e(1:6, 1:6) + fact * (phi * phi');
        M_e(7:12, 7:12) = M_e(7:12, 7:12) + fact * (phi * phi');
        K_e(1:6, 1:6) = K_e(1:6, 1:6) + fact * mu0 * (2 * (dphi(:,1) * dphi(:,1)') + 2*(phi * phi') / (rq^2) ...
            + (dphi(:,2) * dphi(:,2)'));
        K_e(7:12, 7:12) = K_e(7:12, 7:12) + fact * mu0 * (2 * (dphi(:,2) * dphi(:,2)') + ...
            (dphi(:,1) * dphi(:,1)'));
        K_e(1:6, 7:12) = K_e(1:6, 7:12) + fact * mu0 * (dphi(:,2) * dphi(:,1)');
        K_e(7:12, 1:6) = K_e(7:12, 1:6) + fact * mu0 * (dphi(:,1) * dphi(:,2)');
        div_v = dphi(:,1) + phi/rq;
        for k = 1:3
            B_e(k, 1:6) = B_e(k, 1:6) - fact * psi(k) * div_v';
            B_e(k, 7:12) = B_e(k, 7:12) - fact * psi(k) * dphi(:,2)';
        end
        

    end

    idx = [nodes_e_P2, nodes_e_P2 + numP2];
    K(idx, idx) = K(idx, idx) + K_e;
    M(idx, idx) = M(idx, idx) + M_e;
    B(nodes_e_P1, idx) = B(nodes_e_P1, idx) + B_e;
    f(idx) = f(idx) + f_e;
end



%% 8) Conditions aux limites
tol_bc = 1e-5;
bc_r_wall = find(abs(nodes_P2(:,1) - Ro) < tol_bc);  % Paroi r = Ro
bc_r_axis = find(abs(nodes_P2(:,1)) < tol_bc);      % Axe r = 0
bc_z_in = find(abs(nodes_P2(:,2)) < tol_bc);        % Entrée z = 0
bc_z_out = find(abs(nodes_P2(:,2) - L) < tol_bc);   % Sortie z = L
bc_p_out = find(abs(nodes_P1(:,2) - L) < tol_bc, 1); % Un seul nœud pour p = 0

% Conditions de Dirichlet
dirichlet_nodes = [bc_r_wall; bc_r_wall + numP2; ...  % u_r = 0, u_z = 0 sur paroi
                  bc_r_axis; ...                     % u_r = 0 sur axe
                  bc_z_in; bc_z_in + numP2; ...      % u_r = 0, u_z parabolique à l'entrée
                  bc_z_out;  ...                    % u_r = 0
                  bc_p_out + 2*numP2];              % p = 0 sur un nœud

dirichlet_values = [zeros(size(bc_r_wall)); ...       % u_r = 0 sur paroi
                   zeros(size(bc_r_wall)); ...       % u_z = 0 sur paroi
                   zeros(size(bc_r_axis)); ...       % u_r = 0 sur axe
                   zeros(size(bc_z_in)); ...         % u_r = 0 à l'entrée
                   U0 * (1 - (nodes_P2(bc_z_in,1)/Ro).^2); ...  % u_z à l'entrée
                   zeros(size(bc_z_out)); ...        % u_r = 0 à la sortie
                   0];                               % p = 0 sur un nœud
%% 9) Assemblage et rÃ©solution

%% 9) Schéma semi-implicite avec Newton-Raphson

dt = 0.1;
tol = 1e-14;
U_old = U;  % Sauvegarde du pas de temps précédent
t_vec =  0:dt:1;
solutions = sparse(length(t_vec), length(U));

% Stockage des résultats
n_saves = length(t_vec);  % Sauver tous les 0.5s
vorticity_history = zeros(n_saves, numP2);
streamlines_history = zeros(n_saves, numP2);
velocity_history = zeros(n_saves, 2*numP2);
time_saves = t_vec;
save_idx = 1;
for t = t_vec
    fprintf('\n--- Pas de temps t = %.3f ---\n', t+dt);
    
    for iter = 1:10
        % Terme convectif et sa jacobienne
        [C, res_c] = term_c(u_r, u_z, tri_P2, numP2, nodes_P2, numT, piG, gradP2, quad_pts, quad_wts, rho);
        
        % Jacobienne : (M/dt + K + C)
        J = [(rho * M/dt + K + C), B'; B, sparse(numP1, numP1)];
        
        % Résidu : M*(u^{n+1} - u^n)/dt + K*u^{n+1} + N(u^{n+1}) + B'*p = 0
        u_current = [u_r; u_z];
        u_previous = U_old(1:2*numP2);
        residu = [rho * M*(u_current - u_previous)/dt + K*u_current + res_c + B'*p; 
                  B*u_current];
        
        % Conditions aux limites
        for i = 1:length(dirichlet_nodes)
            idx = dirichlet_nodes(i);
            J(idx, :) = 0;
            J(idx, idx) = 1;
            if idx <= 2*numP2
                residu(idx) = u_current(idx) - dirichlet_values(i);
            else
                residu(idx) = p(idx - 2*numP2) - dirichlet_values(i);
            end
        end
        
        % Newton : J * delta_U = -residu
        delta_U = -J \ residu;
        
        % Mise à jour
        du_r = delta_U(1:numP2);
        du_z = delta_U(numP2+1:2*numP2);
        dp = delta_U(2*numP2+1:end);
        
        u_r = u_r + du_r;
        u_z = u_z + du_z;  % CORRECTION
        p = p + dp;
        
        % Test convergence
        norme = norm([du_r; du_z]);
        fprintf('  iter=%d, norm(delta_u)=%e\n', iter, norme);
        
        if norme < tol
            break;
        end
    end
     
    % Mise à jour pour le pas suivant
    if ismember(t+dt, time_saves) && save_idx <= n_saves
        fprintf('Calcul vorticité et lignes de courant...\n');
        [omega, psi] = compute_vorticity_streamlines(u_r, u_z, nodes_P2, tri_P2, numT, piG, gradP2, quad_pts, quad_wts, Ro);
        
        vorticity_history(save_idx, :) = omega';
        streamlines_history(save_idx, :) = psi';
        velocity_history(save_idx, :) = [u_r; u_z]';
        
        save_idx = save_idx + 1;
    end

    U_old = [u_r; u_z; p];
    if mod(t, 10) == 1  % Sauver 1 solution sur 10
        solutions(t,:) = U_old';  % Matrice de stockage à préallouer
    end
end


 %% 10) VÃ©rifications
% Divergence
div_L2 = 0;
for e = 1:numT
    nodes_e_P2 = tri_P2(e, :);
    nodes_e_P1 = tri_P1(e, :);
    coords_P2 = nodes_P2(nodes_e_P2, :);
    r = coords_P2(:, 1); z = coords_P2(:, 2);
    
    Bmat = [r(2)-r(1), r(3)-r(1); z(2)-z(1), z(3)-z(1)];
    detB = abs(det(Bmat));
    invB = inv(Bmat);
    
    div_integral = 0;
    for q = 1:7
        xi = quad_pts(q, 1); 
        eta = quad_pts(q, 2);
        wq = quad_wts(q) * detB;
        phi = piG(xi, eta);
        dphi_ref = gradP2(xi, eta);
        dphi = (invB' * dphi_ref')';
        rq = sum(r .* phi);
        rq = max(rq, 1e-10 * Ro);
        
        dr_ur = dphi(:,1)' * u_r(nodes_e_P2);
        dz_uz = dphi(:,2)' * u_z(nodes_e_P2);
        ur = phi' * u_r(nodes_e_P2);
        div_v = dr_ur + ur / rq + dz_uz;
        
        div_integral = div_integral + wq * div_v^2;
    end
    div_L2 = div_L2 + div_integral;
end
div_L2 = sqrt(div_L2);
fprintf('Norme L2 de la divergence : %.2e\n', div_L2);

% VÃ©rification sur l'axe
axis_nodes = find(abs(nodes_P2(:,1)) < 1e-7);
fprintf('u_z sur l''axe (devrait Ãªtre ~%.2f):\n', U0);
%disp(u_z(axis_nodes));

% Erreur Ã  l'entrÃ©e
inlet_nodes = find(abs(nodes_P2(:,2)) < 1e-5);
u_z_inlet = u_z(inlet_nodes);
u_z_inlet_analytique = U0 * (1 - (nodes_P2(inlet_nodes,1)/Ro).^2);
err_inlet = norm(u_z_inlet - u_z_inlet_analytique) / norm(u_z_inlet_analytique);
fprintf('Erreur u_z Ã  l''entrÃ©e: %.2e\n', err_inlet);

% Erreurs globales
u_z_analytique = U0 * (1 - (nodes_P2(:,1)/Ro).^2);
u_r_analytique = zeros(numP2,1);
err_u_r = norm(u_r - u_r_analytique);
err_u_z = norm(u_z - u_z_analytique) / norm(u_z_analytique);
fprintf('Erreur absolue pour u_r : %.2e\n', err_u_r);
fprintf('Erreur relative pour u_z : %.2e\n', err_u_z);

%% 11) Visualisation
figure;
trisurf(tri_P2(:,1:3), nodes_P2(:,2), nodes_P2(:,1), u_z);
title('Vitesse axiale u_z');
xlabel('z [m]'); ylabel('r [m]'); zlabel('u_z [m/s]');
colorbar;

figure;
trisurf(tri_P2(:,1:3), nodes_P2(:,2), nodes_P2(:,1), u_r);
title('Vitesse radiale u_r');
xlabel('z [m]'); ylabel('r [m]'); zlabel('u_r [m/s]');
colorbar;
omega_final = vorticity_history(end-1, :)';
omega_theo = -2*U0 * nodes_P2(:,1) / Ro^2;
figure;
trisurf(tri_P1, nodes_P1(:,2), nodes_P1(:,1), p);
title('Pression p');
xlabel('z [m]'); ylabel('r [m]'); zlabel('p [Pa]');
colorbar;
figure;
plot(nodes_P2(:,1), omega_final, 'b.', 'DisplayName', 'Numérique');
hold on;
plot(nodes_P2(:,1), omega_theo, 'r-', 'DisplayName', 'Théorique');
xlabel('r [m]'); ylabel('ω_θ [1/s]');
legend; title('Comparaison vorticité');
psi_final = streamlines_history(end-1, :)';
r_nodes = nodes_P2(:,1);
psi_theo = U0 * r_nodes.^2 .* (1 - r_nodes.^2/(2*Ro^2)) / 2;

% Normaliser psi_theo pour avoir même référence que psi_final
psi_theo = psi_theo - mean(psi_theo);
psi_final = psi_final - mean(psi_final);
figure;
plot(nodes_P2(:,1), psi_final, 'b.', 'DisplayName', 'Numérique');
hold on;
plot(nodes_P2(:,1), psi_theo, 'r-', 'DisplayName', 'Théorique');
xlabel('r [m]'); ylabel('ψ [m³/s]');
legend; title('Comparaison fonction de courant');
%% vorticity and streamlines
figure('Position', [100, 100, 1200, 400]);

for i = 1:n_saves-1
    clf;
    
    % Sous-graphique 1 : Vorticité
    subplot(1,3,1);
    omega_i = vorticity_history(i, :)';
    
    % Interpolation sur grille régulière
    z_grid = linspace(min(nodes_P2(:,2)), max(nodes_P2(:,2)), 50);
    r_grid = linspace(min(nodes_P2(:,1)), max(nodes_P2(:,1)), 25);
    [Z_grid, R_grid] = meshgrid(z_grid, r_grid);
    
    F_omega = scatteredInterpolant(nodes_P2(:,2), nodes_P2(:,1), omega_i, 'linear', 'none');
    OMEGA_grid = F_omega(Z_grid, R_grid);
    
    contourf(Z_grid, R_grid, OMEGA_grid, 20, 'EdgeColor', 'none');
    title(sprintf('Vorticité ω_θ (t=%.1fs)', time_saves(i)));
    xlabel('z [m]'); ylabel('r [m]'); 
    colorbar; axis equal;
    
    % Sous-graphique 2 : Lignes de courant
    subplot(1,3,2);
    psi_i = streamlines_history(i, :)';
    
    F_psi = scatteredInterpolant(nodes_P2(:,2), nodes_P2(:,1), psi_i, 'linear', 'none');
    PSI_grid = F_psi(Z_grid, R_grid);
    
    contour(Z_grid, R_grid, PSI_grid, 20);
    title(sprintf('Lignes de courant ψ (t=%.1fs)', time_saves(i)));
    xlabel('z [m]'); ylabel('r [m]');
    axis equal; colorbar;
    
    % Sous-graphique 3 : Champ de vitesse
    subplot(1,3,3);
    u_r_i = velocity_history(i, 1:numP2)';
    u_z_i = velocity_history(i, numP2+1:2*numP2)';
    
    % Magnitude de vitesse
    u_mag = sqrt(u_r_i.^2 + u_z_i.^2);
    F_umag = scatteredInterpolant(nodes_P2(:,2), nodes_P2(:,1), u_mag, 'linear', 'none');
    UMAG_grid = F_umag(Z_grid, R_grid);
    
    contourf(Z_grid, R_grid, UMAG_grid, 20, 'EdgeColor', 'none');
    title(sprintf('|u| (t=%.1fs)', time_saves(i)));
    xlabel('z [m]'); ylabel('r [m]');
    axis equal; colorbar;
    
    % Superposer quelques vecteurs vitesse
    hold on;
    step = 5;  % Sous-échantillonnage
    quiver(nodes_P2(1:step:end,2), nodes_P2(1:step:end,1), ...
           u_z_i(1:step:end), u_r_i(1:step:end), 'k', 'AutoScale', 'on');
    
    sgtitle(sprintf('Évolution temporelle - t = %.1f s', time_saves(i)));
    pause(1);  % Pause pour animation
end

%% Vérifications théoriques


%% fonctions auxiliaires

function [C, res_c] = term_c(u_r, u_z, tri_P2, numP2, nodes_P2, numT, piG, gradP2, quad_pts, quad_wts, rho)
    % Calcul du terme convectif et sa dérivée de Fréchet
    % C : matrice jacobienne du terme convectif (dérivée de Fréchet)
    % res_c : résidu du terme convectif (vecteur)
    
    C = sparse(2*numP2, 2*numP2);
    res_c = zeros(2*numP2, 1);  % Résidu = vecteur
    
    for e = 1:numT
        nodes_e_P2 = tri_P2(e, :);
        coords_P2 = nodes_P2(nodes_e_P2, :);
        r = coords_P2(:, 1); 
        z = coords_P2(:, 2);
        
        % Transformation géométrique
        Bmat = [r(2)-r(1), r(3)-r(1); z(2)-z(1), z(3)-z(1)];
        detB = abs(det(Bmat));
        invB = inv(Bmat);
        
        % Matrices locales
        c_e = zeros(12, 12);    % 6 noeuds P2 × 2 composantes
        res_e = zeros(12, 1);   % Résidu local
        
        for q = 1:7  % Points de quadrature
            xi = quad_pts(q, 1); 
            eta = quad_pts(q, 2);
            
            % Coordonnée radiale et poids
            rq = sum(r .* piG(xi, eta));
            wq = quad_wts(q) * detB;
            fact = 2 * pi * rho * rq * wq;  % Facteur axisymétrique
            
            % Fonctions de forme et gradients
            phi = piG(xi, eta);           % 6×1 - fonctions P2
            dphi_ref = gradP2(xi, eta);   % 6×2 - gradients de référence
            dphi = (invB' * dphi_ref')';  % 6×2 - gradients physiques
            
            % Interpolation des vitesses au point de quadrature
            u_r_q = phi' * u_r(nodes_e_P2);
            u_z_q = phi' * u_z(nodes_e_P2);
            
            % Gradients des vitesses
            grad_u_r_r = dphi(:,1)' * u_r(nodes_e_P2);  % ?u_r/?r
            grad_u_r_z = dphi(:,2)' * u_r(nodes_e_P2);  % ?u_r/?z
            grad_u_z_r = dphi(:,1)' * u_z(nodes_e_P2);  % ?u_z/?r  
            grad_u_z_z = dphi(:,2)' * u_z(nodes_e_P2);  % ?u_z/?z
            
            % === DÉRIVÉE DE FRÉCHET DU TERME CONVECTIF ===
            % Terme convectif : (u·?)u = u_r ?u/?r + u_z ?u/?z
            % Dérivée : D[(u·?)u][v] = (v·?)u + (u·?)v
            
            dphi_r = dphi(:,1);  % 6×1
            dphi_z = dphi(:,2);  % 6×1
            
            % Blocs de la matrice jacobienne 12×12
            % Structure : [?F_r/?u_r, ?F_r/?u_z]
            %            [?F_z/?u_r, ?F_z/?u_z]
            
            % ?F_r/?u_r : dérivée de (u_r ?u_r/?r + u_z ?u_r/?z) par rapport à u_r
            C11 = grad_u_r_r * (phi * phi') + u_r_q * (phi * dphi_r') + u_z_q * (phi * dphi_z');
            
            % ?F_r/?u_z : dérivée de (u_r ?u_r/?r + u_z ?u_r/?z) par rapport à u_z  
            C12 = grad_u_r_z * (phi * phi');
            
            % ?F_z/?u_r : dérivée de (u_r ?u_z/?r + u_z ?u_z/?z) par rapport à u_r
            C21 = grad_u_z_r * (phi * phi');
            
            % ?F_z/?u_z : dérivée de (u_r ?u_z/?r + u_z ?u_z/?z) par rapport à u_z
            C22 =grad_u_z_z * (phi * phi') + u_r_q * (phi * dphi_r') + u_z_q * (phi * dphi_z');
            
            % Assemblage de la matrice locale 12×12
            c_e = c_e + fact * [C11, C12; C21, C22];
            
            % === RÉSIDU (TERME CONVECTIF ÉVALUÉ) ===
            % F_r = u_r ?u_r/?r + u_z ?u_r/?z
            % F_z = u_r ?u_z/?r + u_z ?u_z/?z
            
            conv_r = u_r_q * grad_u_r_r + u_z_q * grad_u_r_z;
            conv_z = u_r_q * grad_u_z_r + u_z_q * grad_u_z_z;
            
            % Contribution au résidu local
            res_e = res_e + fact * [phi * conv_r; phi * conv_z];
        end
        
        % Assemblage global
        idx = [nodes_e_P2, nodes_e_P2 + numP2];
        C(idx, idx) = C(idx, idx) + c_e;
        res_c(idx) = res_c(idx) + res_e;
    end
end

function [omega, psi] = compute_vorticity_streamlines(u_r, u_z, nodes_P2, tri_P2, numT, piG, gradP2, quad_pts, quad_wts, Ro)
    % Calcul de la vorticité et fonction de courant à partir du champ de vitesse
    
    numP2 = size(nodes_P2, 1);
    omega = zeros(numP2, 1);    % Vorticité aux nœuds P2
    psi = zeros(numP2, 1);      % Fonction de courant aux nœuds P2
    node_count = zeros(numP2, 1); % Compteur pour moyennage
    
    %% 1) CALCUL DE LA VORTICITÉ
    % En axisymétrique : ω_θ = ∂u_z/∂r - ∂u_r/∂z
    
    for e = 1:numT
        nodes_e_P2 = tri_P2(e, :);
        coords_P2 = nodes_P2(nodes_e_P2, :);
        r = coords_P2(:, 1); 
        z = coords_P2(:, 2);
        
        % Transformation géométrique
        Bmat = [r(2)-r(1), r(3)-r(1); z(2)-z(1), z(3)-z(1)];
        detB = abs(det(Bmat));
        invB = inv(Bmat);
        
        for q = 1:7
            xi = quad_pts(q, 1); 
            eta = quad_pts(q, 2);
            
            % Fonctions de forme et gradients
            phi = piG(xi, eta);
            dphi_ref = gradP2(xi, eta);
            dphi = (invB' * dphi_ref')';
            
            % Gradients des vitesses au point de quadrature
            du_z_dr = dphi(:,1)' * u_z(nodes_e_P2);  % ∂u_z/∂r
            du_r_dz = dphi(:,2)' * u_r(nodes_e_P2);  % ∂u_r/∂z
            
            % Vorticité au point de quadrature
            omega_q = du_z_dr - du_r_dz;
            
            % Projection sur les nœuds (moyennage)
            rq = sum(r .* phi);
            wq = quad_wts(q) * detB; % Poids ajusté par la taille de l'élément
            for i = 1:6
                node_idx = nodes_e_P2(i);
                omega(node_idx) = omega(node_idx) + phi(i) * omega_q ; % Pondération par r
                node_count(node_idx) = node_count(node_idx) + phi(i)  ;
            end
        end
    end
    
    % Normalisation pour obtenir la moyenne pondérée
    for i = 1:numP2
        if node_count(i) > 1e-10
            omega(i) = omega(i) / node_count(i);
        end
    end
    
    %% 2) CALCUL DE LA FONCTION DE COURANT
    % En axisymétrique : u_r = -1/r ∂ψ/∂z, u_z = 1/r ∂ψ/∂r
    % On résout : Δψ = -r·ω_θ avec ψ = 0 sur les frontières
    
    % Assemblage du Laplacien axisymétrique pour ψ
    L_psi = sparse(numP2, numP2);
    f_psi = zeros(numP2, 1);
    
    for e = 1:numT
        nodes_e_P2 = tri_P2(e, :);
        coords_P2 = nodes_P2(nodes_e_P2, :);
        r = coords_P2(:, 1); 
        z = coords_P2(:, 2);
        
        Bmat = [r(2)-r(1), r(3)-r(1); z(2)-z(1), z(3)-z(1)];
        detB = abs(det(Bmat));
        invB = inv(Bmat);
        
        L_e = zeros(6, 6);
        f_e = zeros(6, 1);
        
        for q = 1:7
            xi = quad_pts(q, 1); 
            eta = quad_pts(q, 2);
            
            rq = sum(r .* piG(xi, eta));
            rq = max(rq, 1e-10 * Ro); % Éviter division par 0
            wq = quad_wts(q) * detB;
            
            phi = piG(xi, eta);
            dphi_ref = gradP2(xi, eta);
            dphi = (invB' * dphi_ref')';
            
            % Laplacien axisymétrique : ∇²ψ = ∂²ψ/∂r² + (1/r)∂ψ/∂r + ∂²ψ/∂z²
            fact = 2 * pi * rq * wq;
            L_e = L_e + fact * (dphi * dphi');  % Terme ∇ψ·∇φ
            L_e = L_e + fact * (1/rq) * (dphi(:,1) * phi');  % Terme (1/r)∂ψ/∂r
            
            % Second membre : -r·ω
            omega_q = phi' * omega(nodes_e_P2);
            f_e = f_e - fact * rq * omega_q * phi;
        end
        
        L_psi(nodes_e_P2, nodes_e_P2) = L_psi(nodes_e_P2, nodes_e_P2) + L_e;
        f_psi(nodes_e_P2) = f_psi(nodes_e_P2) + f_e;
    end
    
    % Conditions aux limites : ψ = 0 sur toutes les frontières
    tol_bc = 1e-5;
    bc_all = find(abs(nodes_P2(:,1) - Ro) < tol_bc | ...  % Paroi r = Ro
                  abs(nodes_P2(:,1)) < tol_bc | ...       % Axe r = 0
                  abs(nodes_P2(:,2)) < tol_bc | ...       % Entrée z = 0
                  abs(nodes_P2(:,2) - 0.02) < tol_bc);    % Sortie z = L
    
    % Application des conditions de Dirichlet
    for i = 1:length(bc_all)
        idx = bc_all(i);
        L_psi(idx, :) = 0;
        L_psi(idx, idx) = 1;
        f_psi(idx) = 0;
    end
    
    % Résolution
    psi = L_psi \ f_psi;

end
