% Simulation de Poiseuille avec multiplicateur de Lagrange pour débit fixe
clc; clear; close all;

% Paramètres
rho = 1050; mu0 = 0.0056; Ro = 0.015; L = 0.05; U0 = 0.1; 
mu_inf = 0.0035; lambda = 3.43;  n = 0.3568; tol = 1e-12; n_iter = 10; a = 2.0;

carreau_params.eta_inf = mu_inf;   
carreau_params.eta_0 = mu0;        
carreau_params.lambda = lambda;       
carreau_params.n = n;
carreau_params.a = a; 

%% 2) Maillage P2-P1 avec PDE Toolbox
n_r = 8;   % Subdivisions radiales
n_z = 4;    % Subdivisions axiales
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

fprintf('Nombre d''éléments : %d\n', numT);
fprintf('Nombre de noeuds P1 : %d\n', numP1);
fprintf('Nombre de noeuds P2 : %d\n', numP2);

%% Fonctions de forme P2 sur triangle de référence
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

%% Initialisation des solutions
u_r = zeros(numP2,1);   
u_z = U0 * (1 - (nodes_P2(:,1)/Ro).^2);  % Profil parabolique initial
P = zeros(numP1,1);     
U = [u_r; u_z; P];

%% Quadrature de Gauss à 7 points
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


%% Boucle de Newton avec contrainte de débit
m =(n-1)/a; 
k = mu0 - mu_inf;
muk = @(dot_gamma) mu_inf + k * (1 + (lambda * dot_gamma)^a)^m;
dmuk =@(dot_gamma) k * m * a * (lambda^a) * (dot_gamma.^(a-1)) .* ...
                           (1 + (lambda * dot_gamma).^a).^(m - 1);


%% Conditions aux limites
tol_bc = 1e-5;
bc_r_wall = find(abs(nodes_P2(:,1) - Ro) < tol_bc);
bc_r_axis = find(abs(nodes_P2(:,1)) < tol_bc);
bc_z_in = find(abs(nodes_P2(:,2)) < tol_bc);
bc_z_out = find(abs(nodes_P2(:,2) - L) < tol_bc);
bc_p_out = find(abs(nodes_P1(:,2) - L) < tol_bc, 1);

% PROFIL D'ENTRÉE PROPREMENT CALCULÉ
r_inlet_sorted = sort(nodes_P2(bc_z_in,1));
% % Réorganiser selon l'ordre original des nœuds bc_z_in
r_inlet_original = nodes_P2(bc_z_in,1);
% u_in = zeros(size(r_inlet_original));


% Conditions de Dirichlet
dirichlet_nodes = [bc_r_wall; bc_r_wall + numP2; ...  % u_r = 0, u_z = 0 sur paroi
                  bc_r_axis; ...                     % u_r = 0 sur axe
                  bc_z_in;  ...                     % u_r = 0
                  bc_z_out;  ...                    % u_r = 0 à la sortie
                  bc_p_out + 2*numP2;...            % p = 0 sur un nœud
                  bc_z_in + numP2];                 % u_z profil à l'entrée  

dirichlet_values = [zeros(size(bc_r_wall)); ...       % u_r = 0 sur paroi
                   zeros(size(bc_r_wall)); ...       % u_z = 0 sur paroi
                   zeros(size(bc_r_axis)); ...       % u_r = 0 sur axe
                   zeros(size(bc_z_in)); ...         % u_r = 0 à l'entrée
                   zeros(size(bc_z_out)); ...        % u_r = 0 à la sortie
                   0                               % p = 0 sur un nœud
                   U0 * (1 - (nodes_P2(bc_z_in,1)/Ro).^2); ...         % u_r = 0 à l'entrée
                   ];                               

%% debit pulsatil


% Womersley : Débit pulsé (exemple)
f = 1.2;                    % Hz - fréquence cardiaque (72 bpm)
omega = 2*pi*f;             % rad/s - pulsation
T_cycle = 1/f;              % s - période cardiaque
Q_mean = 83e-6; % débit moyen
Q_amplitude = 40e-6; % amplitude pulsatile
Q_target = @(t) Q_mean + Q_amplitude * sin(omega * t);
%dt = T_cycle/50;            % Pas de temps (50 points par cycle)
dt = 0.1;
t_final = T_cycle;        % Simuler 2 cycles cardiaques
t_vec = 0:dt:t_final;
%%
u_in = zeros(size(bc_z_in)); % Vecteur pour les vitesses d'entrée
fprintf('\n=== SIMULATION TEMPORELLE AVEC WOMERSLEY ===\n');
fprintf('Période cardiaque: %.3f s, dt=%.4f s\n', T_cycle, dt);
fprintf('Nombre de Womersley: α = %.2f\n', Ro*sqrt(omega*rho/mu0));
solutions = sparse(length(t_vec), length(U));
%%
dt = 0.1;
U_old = U;  % Sauvegarde du pas de temps précédent
for i_time = 1:length(t_vec)
    t_current = t_vec(i_time);
    fprintf('\n--- t = %.4f s (%.1f%% du cycle) ---\n', ...
            t_current, 100*mod(t_current, T_cycle)/T_cycle);
    
    % 1) CALCUL DU PROFIL WOMERSLEY
    Q_instant = Q_target(t_current);
    u_womersley_sorted = womersley(r_inlet_sorted, Ro, Q_instant,t_current, carreau_params);
    
    % 2) RÉORGANISATION SELON L'ORDRE ORIGINAL DES NOEUDS
    for j = 1:length(bc_z_in)
        [~, pos] = min(abs(r_inlet_sorted - r_inlet_original(j)));
        u_in(j) = u_womersley_sorted(pos);
    end
    
    % 3) MISE À JOUR DES CONDITIONS DE DIRICHLET
    dirichlet_values(end-length(bc_z_in)+1:end) = u_in';
    
    fprintf('Débit cible: %.2f ml/s, u_max: %.3f m/s\n', ...
            Q_instant*1e6, max(u_in));
     iter = 1;
    while iter <= n_iter
        % Extraction des solutions courantes
        u_r_ = U(1:numP2);
        u_z_ = U(numP2+1:2*numP2);
        p_ = U(2*numP2+1:end);


        % Initialisation des matrices globales
        
        M = sparse(2*numP2, 2*numP2); 
        K = sparse(2*numP2, 2*numP2);  
        J = sparse(2*numP2, 2*numP2);
        B = sparse(numP1, 2*numP2);          
        f = zeros(2*numP2, 1);    
        %%
       
        %% Assemblage des matrices élémentaires
        for e = 1:numT
            nodes_e_P2 = tri_P2(e, :);
            nodes_e_P1 = tri_P1(e, :);
            coords_P2 = nodes_P2(nodes_e_P2, :);
            r = coords_P2(:, 1); z = coords_P2(:, 2);
            Bmat = [r(2)-r(1), r(3)-r(1); z(2)-z(1), z(3)-z(1)];
            detB = abs(det(Bmat));
            invB = inv(Bmat);

            % Matrices élémentaires
            M_e = zeros(12, 12);
            K_e = zeros(12, 12);  
            B_e = zeros(3, 12); 
            J_e = zeros(12, 12);
            f_e = zeros(12, 1);
            c_e = zeros(12, 12);
            % Solutions locales
            ur_e = u_r_(nodes_e_P2);
            uz_e = u_z_(nodes_e_P2);

            % Régularisation
            eps_r = 1e-10 * Ro;
            eps_gamma = 1e-6;
            h_elem = sqrt(detB);  % Taille d'élément
            

            for q = 1:7
                xi = quad_pts(q, 1); 
                eta = quad_pts(q, 2);

                % Coordonnées et poids
                rq = sum(r .* piG(xi, eta));
                rq_reg = max(rq, eps_r);
                wq = quad_wts(q) * detB;
                fact = 2 * pi * rq_reg * wq;

                % Fonctions de forme
                phi = piG(xi, eta);
                dphi_ref = gradP2(xi, eta);
                dphi = (invB' * dphi_ref')';
                psi = [psi1(xi, eta); psi2(xi, eta); psi3(xi, eta)];
                               
                % matrice de masse
                M_e(1:6, 1:6) = M_e(1:6, 1:6) + fact * (phi * phi');
                M_e(7:12, 7:12) = M_e(7:12, 7:12) + fact * (phi * phi');
                % Interpolation des vitesses
                ur_intp = phi' * ur_e;

                % Tenseur de déformation
                drr = dphi(:,1)' * ur_e;
                dzz = dphi(:,2)' * uz_e;
                doo = ur_intp / rq_reg;
                drz = 0.5 * (dphi(:,1)' * uz_e + dphi(:,2)' * ur_e);

                % Taux de cisaillement
                d = [drr, 0, drz; 0, doo, 0; drz, 0, dzz];
                dot_gamma = sqrt(2) * norm(d, 'fro');
                dot_gamma_reg = max(dot_gamma, eps_gamma);

                % Viscosité
                mu = muk(dot_gamma_reg);

                % ASSEMBLAGE PARTIE LINÉAIRE
                terme_rr = 2 * (dphi(:,1) * dphi(:,1)') + 2*(phi * phi') / (rq_reg^2) ...
                          + (dphi(:,2) * dphi(:,2)');
                K_e(1:6, 1:6) = K_e(1:6, 1:6) + fact * mu * terme_rr;

                terme_zz = 2 * (dphi(:,2) * dphi(:,2)') + (dphi(:,1) * dphi(:,1)');
                K_e(7:12, 7:12) = K_e(7:12, 7:12) + fact * mu * terme_zz;

                terme_rz = dphi(:,2) * dphi(:,1)';
                K_e(1:6, 7:12) = K_e(1:6, 7:12) + fact * mu * terme_rz;
                K_e(7:12, 1:6) = K_e(7:12, 1:6) + fact * mu * terme_rz';

                % ASSEMBLAGE PARTIE NON-LINÉAIRE

                dmu = dmuk(dot_gamma_reg);
                coeff_nl = 4 * (dmu / dot_gamma_reg) * fact;

                % Vecteur b : D : ?D
                b = zeros(12, 1);
                b(1:6) = drr * dphi(:,1) + doo * (phi / rq_reg) + drz * dphi(:,2);
                b(7:12) = dzz * dphi(:,2) + drz * dphi(:,1);

                % Contribution non-linéaire
                J_e = J_e + coeff_nl * (b * b');

                penalty_div = 1e-8;
                % MATRICE DE DIVERGENCE 
                div_v = dphi(:,1) + phi / rq_reg;

                for k = 1:3
                    B_e(k, 1:6) = B_e(k, 1:6) - fact * psi(k) * div_v';
                    B_e(k, 7:12) = B_e(k, 7:12) - fact * psi(k) * dphi(:,2)';
                end
                
            end

            % JACOBIEN COMPLET
            J_e = K_e + J_e ;

            % Assemblage global
            idx = [nodes_e_P2, nodes_e_P2 + numP2];
            K(idx, idx) = K(idx, idx) + K_e;
            J(idx, idx) = J(idx, idx) + J_e;
            B(nodes_e_P1, idx) = B(nodes_e_P1, idx) + B_e;
            M(idx, idx) = M(idx, idx) + M_e;
            
            
        end


        %% CONSTRUCTION DU SYSTÈME 
        u_previous = U_old(1:2*numP2);
        u_r_old = u_previous(1:numP2);
        u_z_old = u_previous((numP2+1):end);
        u_current = [u_r_; u_z_];
        % Résidu des équations de quantité de mouvement (SANS terme de contrainte ici)
        [C, res_c] = term_c(u_r_, u_z_, tri_P2, numP2, nodes_P2, numT, piG, gradP2, quad_pts, quad_wts, rho);

        % Résidu momentum complet (matrice*vecteur + termes sources)
        res_momentum = rho * M * (u_current - u_previous)/dt + K * u_current + res_c + B' * p_;

        % Résidu de l'équation de continuité
        res_continuity = B * u_current;
        % Résidu global
        res = [res_momentum; res_continuity];
    
        J = J + rho * M/dt + C;
        
        % ASSEMBLAGE CORRECT DU SYSTÈME 
         % Système jacobien standard
        % [J   B']
        % [B   0 ]
        J_ = [J, B'; B, zeros(numP1)];
        for i = 1:length(dirichlet_nodes)
            idx = dirichlet_nodes(i);

            % Ajuster le résidu AVANT d'annuler la colonne
            res = res - J_(:, idx) * (U(idx) - dirichlet_values(i));

            % Annuler ligne ET colonne (méthode stricte)
            J_(idx, :) = 0;
            J_(:, idx) = 0;
            J_(idx, idx) = 1;
            res(idx) = U(idx) - dirichlet_values(i);
        end


        %% Résolution et mise à jour
        delta_U = J_ \ (-res);
        U = U + delta_U;

        %% Test de convergence
        norme_residu = norm(res);
        norme_delta = norm(delta_U(1:2*numP2));
        norme_solution = norm([u_r_; u_z_]);

        tol_abs = 1e-13;  % Tolérance plus stricte
        tol_rel = 1e-12;  % Tolérance plus stricte

        converged = (norme_residu < tol_abs) && ...
                   (norme_delta < tol_rel * max(norme_solution, 1));

        fprintf('Iter %d: ||res|| = %.2e, ||delta_u|| = %.2e, div·u = %.2e\n', ...
                iter, norme_residu, norme_delta,...
                norm(B * [u_r_; u_z_])); % Norme de la divergence

        if converged
            fprintf('Convergence atteinte en %d itérations\n', iter);
            break;
        end

        iter = iter + 1;
    end
    % Mise à jour pour le pas suivant
    U_old = U;
    % Calcul du débit effectif
   % Calcul du débit effectif
    % Calcul du débit effectif
    u_z_current = U(numP2+1:2*numP2);
    r_inlet_nodes = nodes_P2(bc_z_in,1);
    u_z_inlet_nodes = u_z_current(bc_z_in);

    % Trier les noeuds pour l'intégration
    [r_sorted, idx_sort] = sort(r_inlet_nodes);
    u_sorted = u_z_inlet_nodes(idx_sort);

    % Intégration par trapèzes
    Q_effectif = 2*pi * trapz(r_sorted, u_sorted .* r_sorted);
    fprintf('Débit effectif: %.2f ml/s (écart: %.1f%%)\n', ...
            Q_effectif*1e6, 100*abs(Q_effectif-Q_instant)/Q_instant);
    
    % Optionnel: sauvegarde des solutions pour post-traitement
    if mod(i_time, 10) == 1  % Sauver 1 solution sur 10
        solutions(i_time,:) = U';  % Matrice de stockage à préallouer
    end
    %% Vérification de la divergence
    u_r = U(1:numP2);
    u_z = U(numP2+1:2*numP2);
    div_L2 = 0;
    for e = 1:numT
        nodes_e_P2 = tri_P2(e, :);
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

            div_integral = div_integral + 2 * pi * rq * wq * div_v^2;
        end
        div_L2 = div_L2 + div_integral;
    end
    div_L2 = sqrt(div_L2);
    fprintf('Norme L2 de la divergence : %.2e\n', div_L2);
end
%% Extraction des solutions finales
u_r = U(1:numP2);
u_z = U(numP2+1:2*numP2);
p = U(2*numP2+1:end);


%% Vérification de la divergence
div_L2 = 0;
for e = 1:numT
    nodes_e_P2 = tri_P2(e, :);
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
        
        div_integral = div_integral + 2 * pi * rq * wq * div_v^2;
    end
    div_L2 = div_L2 + div_integral;
end
div_L2 = sqrt(div_L2);
fprintf('Norme L2 de la divergence : %.2e\n', div_L2);

% --- Calcul divergence élémentaire et visualisation ---
div_elem = zeros(numT,1);
max_div_elem = 0;
for e = 1:numT
    nodes_e_P2 = tri_P2(e, :);
    coords_P2 = nodes_P2(nodes_e_P2, :);
    r = coords_P2(:,1); z = coords_P2(:,2);
    Bmat = [r(2)-r(1), r(3)-r(1); z(2)-z(1), z(3)-z(1)];
    detB = abs(det(Bmat));
    invB = inv(Bmat);
    div_integral = 0;
    for q = 1:7
        xi = quad_pts(q,1); eta = quad_pts(q,2);
        wq = quad_wts(q)*detB;
        phi = piG(xi,eta);
        dphi_ref = gradP2(xi,eta);
        dphi = (invB' * dphi_ref')';
        rq = sum(r .* phi);
        rq = max(rq, 1e-12); % régularisation locale pour le post-traitement
        dr_ur = dphi(:,1)' * u_r(nodes_e_P2);
        dz_uz = dphi(:,2)' * u_z(nodes_e_P2);
        ur = phi' * u_r(nodes_e_P2);
        div_v = dr_ur + ur / rq + dz_uz;
        div_integral = div_integral + wq * (div_v)^2;
    end
    div_elem(e) = sqrt(div_integral);
    max_div_elem = max(max_div_elem, div_elem(e));
end

figure;
patch('Faces', tri_P2(:,1:3), ...
      'Vertices', nodes_P2(:,[2 1]), ... % (z,r) pour cohérence
      'FaceVertexCData', div_elem, ...
      'FaceColor', 'flat', ...
      'EdgeColor', 'none');
colorbar;
axis equal tight;
xlabel('z [m]'); ylabel('r [m]');
title('Divergence RMS par élément (P2)');


%% Visualisation
figure;
trisurf(tri_P2(:,1:3), nodes_P2(:,2), nodes_P2(:,1), u_z, 'EdgeColor', 'none');
title('Vitesse axiale u_z');
xlabel('z [m]'); ylabel('r [m]'); zlabel('u_z [m/s]');
colorbar;

figure;
trisurf(tri_P2(:,1:3), nodes_P2(:,2), nodes_P2(:,1), u_r, 'EdgeColor', 'none');
title('Vitesse radiale u_r');
xlabel('z [m]'); ylabel('r [m]'); zlabel('u_r [m/s]');
colorbar;

figure;
trisurf(tri_P1, nodes_P1(:,2), nodes_P1(:,1), p,'EdgeColor', 'none');
title('Pression p');
xlabel('z [m]'); ylabel('r [m]'); zlabel('p [Pa]');
colorbar;

% % Profil à l'entrée
% figure;
% inlet_nodes = find(abs(nodes_P2(:,2)) < 1e-5);
% r_inlet = nodes_P2(inlet_nodes, 1);
% u_z_inlet = u_z(inlet_nodes);
% [r_sorted, idx] = sort(r_inlet);
% u_z_sorted = u_z_inlet(idx);
% 
% plot(r_sorted, u_z_sorted, 'bo-', 'LineWidth', 1.5);
% hold on;
% r_analytical = linspace(0, Ro, 100);
% u_z_analytical = fem_in(r_analytical, Ro, G, carreau_params);
% plot(r_analytical, u_z_analytical, 'r-', 'LineWidth', 2);
% xlabel('r [m]');
% ylabel('u_z [m/s]');
% title('Profil de vitesse à l''entrée');
% legend('Numérique', 'Analytique', 'Location', 'best');
% grid on;

%% Fonctions auxiliaires
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


function [u_z_inlet] = womersley(r_nodes, R, Q_target, t, carreau_params)
%--------------------------------------------------------------------------
% Profil de Womersley (axial) pour une entrée axisymétrique
%   r_nodes    : vecteur radial des noeuds [m]
%   R          : rayon de l’aorte [m]
%   Q_target   : débit instantané imposé [m^3/s]
%   t          : temps courant [s]
%   carreau_params : struct avec champs
%          eta_0 , eta_inf , lambda , n , a
% OUT :
%   u_z_inlet  : vitesse axiale aux noeuds [m/s]
%   u_r_inlet  : vitesse radiale aux noeuds (≈0 sauf ∂z) [m/s]
%--------------------------------------------------------------------------

    %% --- paramètres physiques
    rho    = 1050;                          % kg/m^3
    mu_inf = carreau_params.eta_inf;        % Pa·s
    mu0    = carreau_params.eta_0;          % Pa·s
    lam    = carreau_params.lambda;         % s
    n      = carreau_params.n;
    a      = carreau_params.a;

    %% --- fréquence cardiaque (exemple)
    f     = 1.2;                % Hz
    omega = 2*pi*f;             % rad/s

    %% --- viscosité effective (Carreau)
    Uavg  = abs(Q_target)/(pi*R^2);     % vitesse moyenne
    gdot  = max(4*Uavg/R, 1e-6);        % taux cisaillement moyen
    mu_eff = mu_inf + (mu0 - mu_inf) * (1 + (lam*gdot)^a)^((n-1)/a);
    nu = mu_eff / rho;

    %% --- nombre de Womersley
    alpha = R * sqrt(omega/nu);

    %% --- facteur complexe (solution Womersley)
    % i^(3/2) = exp(i*3π/4)
    icmplx = exp(1i*3*pi/4);
    J0a    = besselj(0, icmplx*alpha);
    phi_r  = 1 - besselj(0, icmplx*alpha*(r_nodes/R)) ./ J0a;

    %% --- phase temporelle
    phi_t  = phi_r .* exp(1i*omega*t);
    phi    = real(phi_t);      % on prend la partie réelle pour le champ axial

    %% --- normalisation pour respecter Q_target
    [r_sort, idx] = sort(r_nodes(:));
    phi_sort      = phi(idx);
    Int_phi       = 2*pi * trapz(r_sort, phi_sort .* r_sort);

    if abs(Int_phi) < 1e-12
        % repli parabolique si dégénérescence
        warning('Profil Womersley dégénéré – repli sur parabolique.');
        umax = 2*Q_target / (pi*R^2);
        u_z_inlet = umax * (1 - (r_nodes/R).^2);
        u_r_inlet = zeros(size(r_nodes));
        u_z_inlet(r_nodes>R) = 0;
        return;
    end

    A = Q_target / Int_phi;
    u_z_inlet = A * phi;
    u_z_inlet(r_nodes>R) = 0;

    %% --- débit de contrôle
    Q_check = 2*pi * trapz(r_sort, ...
                interp1(r_nodes,u_z_inlet,r_sort,'linear','extrap').*r_sort);
    if abs(Q_check - Q_target) > 0.05*abs(Q_target)
        warning('Écart débit: calc=%.3e  cible=%.3e', Q_check, Q_target);
    end

    %% --- petite composante radiale pour ∂z u_z ≠ 0
    Qamp   = 0.4*abs(Q_target);          % estimation amplitude pulsatile
    dQdt   =  Qamp * omega * cos(omega*t);
    duz_dz = -dQdt / (pi*R^2);           % approx uniforme

    u_r_inlet = zeros(size(r_nodes));
    for k = 1:numel(r_nodes)
        r = r_nodes(k);
        if r>1e-12
            u_r_inlet(k) = -duz_dz * r / 2;  % incompressibilité
        end
    end
    u_r_inlet(r_nodes>R) = 0;
end