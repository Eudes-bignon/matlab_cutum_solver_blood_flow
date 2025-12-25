function [solutions, simulation_params, geometry_params, t_vec] = compute_ns( carreau_params, U0, rho) 

    % Paramètres
    Ro = 0.015; 
    n_iter = 10;

    mu_inf = carreau_params.eta_inf;   
    mu0 = carreau_params.eta_0 ;        
    lambda =  carreau_params.lambda;       
    n = carreau_params.n ;
    a = carreau_params.a; 

    %% 2) Maillage P2-P1 avec PDE Toolbox
    L = 0.05;          % Longueur totale (m)
    R_inlet = Ro;   % Rayon sain (entrée/sortie) (m)
    R_max = 0.025;     % Rayon maximal de l'anévrisme (m)
    L_aneurysm = 0.015; % Longueur anévrisme (m)
    x_start_aneurysm = (L - L_aneurysm) / 2;
    x_end_aneurysm   = x_start_aneurysm + L_aneurysm;
    x_mid_aneurysm   = x_start_aneurysm + L_aneurysm / 2;

    % 2) Construction du contour sans duplication
    % Nombre de points pour la courbe supérieure
    n_curve = 1000;

    % Points de contrôle pour la spline
    x_points = [0, x_start_aneurysm, x_mid_aneurysm, x_end_aneurysm, L];
    y_points = [R_inlet, R_inlet, R_max, R_inlet, R_inlet];
    pp = pchip(x_points, y_points);

    % Courbe supérieure (sans les extrémités pour éviter duplication)
    x_curve = linspace(0, L, n_curve);
    y_curve = ppval(pp, x_curve);
    y_curve = max(y_curve, 0.9*R_inlet);
    y_curve([1,end]) = R_inlet;

    % 3') Construction du contour fermé propre
    % Méthode : construire le contour dans l'ordre, sans duplication
    contour_x = [];
    contour_y = [];

    % 1) Axe inférieur de gauche à droite (z=0, r=0)
    x_bottom = linspace(0, L, 100);
    y_bottom = zeros(size(x_bottom));
    contour_x = [contour_x, x_bottom];
    contour_y = [contour_y, y_bottom];

    % 2) Paroi droite (remontée verticale à z=L)
    y_right = linspace(0, R_inlet, 50);
    x_right = L * ones(size(y_right));
    % Retirer le premier point pour éviter duplication avec l'axe
    contour_x = [contour_x, x_right(2:end)];
    contour_y = [contour_y, y_right(2:end)];

    % 3) Courbe supérieure de droite à gauche
    x_top = fliplr(x_curve);
    y_top = fliplr(y_curve);
    % Retirer le premier point pour éviter duplication avec la paroi droite
    contour_x = [contour_x, x_top(2:end)];
    contour_y = [contour_y, y_top(2:end)];

    % 4) Paroi gauche (descente verticale à z=0)
    y_left = linspace(R_inlet, 0, 50);
    x_left = zeros(size(y_left));
    % Retirer les deux extrémités pour éviter duplication
    contour_x = [contour_x, x_left(2:end-1)];
    contour_y = [contour_y, y_left(2:end-1)];

    % 4') Nettoyage supplémentaire des points trop proches
    tol_ = 1e-10;
    dx = diff([contour_x, contour_x(1)]);
    dy = diff([contour_y, contour_y(1)]);
    distances = sqrt(dx.^2 + dy.^2);
    keep_idx = distances > tol_;
    contour_x = contour_x(keep_idx);
    contour_y = contour_y(keep_idx);

    % 5) Création polyshape propre
    pg = polyshape(contour_x, contour_y, 'Simplify', false, 'KeepCollinearPoints', false);

    % 6) Vérification de la validité du polyshape
    if ~isempty(pg.Vertices)
        fprintf('✅ Polyshape créé avec succès, %d sommets\n', size(pg.Vertices, 1));
    else
        error('❌ Échec de création du polyshape');
    end

    % 7) Géométrie PDE + maillage propre
    tr = triangulation(pg);
    model = createpde(1);
    geometryFromMesh(model, tr.Points', tr.ConnectivityList');

    % Maillage quadratique, pas trop dense
    mesh = generateMesh(model, 'Hmax', R_inlet/16, 'GeometricOrder', 'quadratic');

    nodes_P2 = mesh.Nodes';          
    tri_P2 = mesh.Elements(1:6,:)';  
    numT = size(tri_P2, 1);          
    tri_P1 = mesh.Elements(1:3,:)';  
    all_P1_nodes = unique(tri_P1(:));
    nodes_P1 = nodes_P2(all_P1_nodes, :);
    numP1 = size(nodes_P1, 1);
    numP2 = size(nodes_P2, 1);

    % fprintf('Nombre d''éléments : %d\n', numT);
    % fprintf('Nombre de noeuds P1 : %d\n', numP1);
    % fprintf('Nombre de noeuds P2 : %d\n', numP2);

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
    u_z = U0 * (1 - (nodes_P2(:,2)/Ro).^2);  % Profil parabolique initial
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
    tol_bc = 1e-4;

    bc_r_axis = find(abs(nodes_P2(:,2)) < tol_bc);                    % r = 0 (axe)
    bc_r_wall = find(abs(nodes_P2(:,2) - ppval(pp, nodes_P2(:,1))) < tol_bc & ...
                    nodes_P2(:,2) > tol_bc);                          % r = R(x) (paroi)
    bc_z_in   = find(abs(nodes_P2(:,1)) < tol_bc);                   % x = 0 (entrée)
    bc_z_out  = find(abs(nodes_P2(:,1) - L) < tol_bc);              % x = L (sortie)

    % Point de référence pour la pression
    bc_p_out = find(abs(nodes_P1(:,1) - L) < tol_bc, 1);
    % PROFIL D'ENTRÉE PROPREMENT CALCULÉ
    r_inlet_sorted = sort(nodes_P2(bc_z_in,2));
    % % Réorganiser selon l'ordre original des nœuds bc_z_in
    r_inlet_original = nodes_P2(bc_z_in,2);
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
                       zeros(size(bc_z_in)); ...         % u_z = 0 à l'entrée
                       ];                               

    %% debit pulsatil

    % Womersley : Débit pulsé (exemple)
    f = 1.2;                    % Hz - fréquence cardiaque (72 bpm)
    omega = 2*pi*f;             % rad/s - pulsation
    T_cycle = 1/f;              % s - période cardiaque
    CO = 6.5;                  % L/min
    Q_mean = (CO/60) * 1e-3; % m^3/s
    
    % Option 1 : heuristique (frac)
    amp_frac = 1.2;          % 0.3..0.8 typique au repos
    Q_amplitude = amp_frac * Q_mean;
       % amplitude pulsatile
    Q_target = @(t) Q_mean + Q_amplitude * sin(omega * t);
    dt = T_cycle/10;            % Pas de temps (50 points par cycle)
    t_final = T_cycle;        % Simuler 2 cycles cardiaques
    t_vec = 0:dt:t_final;
    temp_theorique = [0, T_cycle/4, T_cycle/2, 3*T_cycle/4];
    [~, temp_idx] = arrayfun(@(t) min(abs(t_vec - t)), temp_theorique);
    temp = t_vec(temp_idx);
    %%
    U_old = U;  % Solution du pas précédent
    u_in = zeros(size(bc_z_in)); % Vecteur pour les vitesses d'entrée
    fprintf('\n=== SIMULATION TEMPORELLE AVEC WOMERSLEY ===\n');
    fprintf('Période cardiaque: %.3f s, dt=%.4f s\n', T_cycle, dt);
    fprintf('Nombre de Womersley: α = %.2f\n', Ro*sqrt(omega*rho/mu0));
    % Stockage temporel (pour un cycle complet, e.g., le deuxième)
    num_steps_cycle = length(t_vec(t_vec > T_cycle & t_vec <= 2*T_cycle));  % Nb pas dans le 2e cycle
    u_r_time = zeros(numP2, num_steps_cycle);
    u_z_time = zeros(numP2, num_steps_cycle);
    p_time = zeros(numP1, num_steps_cycle);
    idx_cycle = 0;
    for i_time = 1:(length(t_vec) - 1)
        t_current = t_vec(i_time);
        fprintf('\n--- t = %.4f s (%.1f%% du cycle) ---\n', ...
                t_current, 100*mod(t_current, T_cycle)/T_cycle);

        % 1) CALCUL DU PROFIL WOMERSLEY
        Q_instant = Q_target(t_current);
        u_womersley_sorted = womersley(r_inlet_sorted, Ro, Q_instant, carreau_params);

        % 2) RÉORGANISATION SELON L'ORDRE ORIGINAL DES NOEUDS
        for j = 1:length(bc_z_in)
            [~, pos] = min(abs(r_inlet_sorted - r_inlet_original(j)));
            u_in(j) = u_womersley_sorted(pos);
        end

        % 3) MISE À JOUR DES CONDITIONS DE DIRICHLET
        dirichlet_values(end-length(bc_z_in)+1:end) = u_in';

    %     fprintf('Débit cible: %.2f ml/s, u_max: %.3f m/s\n', ...
    %             Q_instant*1e6, max(u_in));

        % 4) BOUCLE DE NEWTON (votre code existant)
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
                z = coords_P2(:, 1);  % Position axiale
                r = coords_P2(:, 2);  % Position radiale

                Bmat = [z(2)-z(1), z(3)-z(1); r(2)-r(1), r(3)-r(1)];
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
                    drr = dphi(:,2)' * ur_e;
                    dzz = dphi(:,1)' * uz_e;
                    doo = ur_intp / rq_reg;
                    drz = 0.5 * (dphi(:,2)' * uz_e + dphi(:,1)' * ur_e);

                    % Taux de cisaillement
                    d = [drr, 0, drz; 0, doo, 0; drz, 0, dzz];
                    dot_gamma = sqrt(2) * norm(d, 'fro');
                    dot_gamma_reg = max(dot_gamma, eps_gamma);

                    % Viscosité
                    mu = muk(dot_gamma_reg);

                    % ASSEMBLAGE PARTIE LINÉAIRE
                    terme_rr = 2 * (dphi(:,2) * dphi(:,2)') + 2*(phi * phi') / (rq_reg^2) ...
                              + (dphi(:,1) * dphi(:,1)');
                    K_e(1:6, 1:6) = K_e(1:6, 1:6) + fact * mu * terme_rr;

                    terme_zz = 2 * (dphi(:,1) * dphi(:,1)') + (dphi(:,2) * dphi(:,2)');
                    K_e(7:12, 7:12) = K_e(7:12, 7:12) + fact * mu * terme_zz;

                    terme_rz = dphi(:,1) * dphi(:,2)';
                    K_e(1:6, 7:12) = K_e(1:6, 7:12) + fact * mu * terme_rz;
                    K_e(7:12, 1:6) = K_e(7:12, 1:6) + fact * mu * terme_rz';

                    % ASSEMBLAGE PARTIE NON-LINÉAIRE

                    dmu = dmuk(dot_gamma_reg);
                    coeff_nl = 4 * (dmu / dot_gamma_reg) * fact;

                    % Vecteur b : D : ?D
                    b = zeros(12, 1);
                    b(1:6) = drr * dphi(:,2) + doo * (phi / rq_reg) + drz * dphi(:,1);
                    b(7:12) = dzz * dphi(:,1) + drz * dphi(:,2);

                    % Contribution non-linéaire
                    J_e = J_e + coeff_nl * (b * b');

                    penalty_div = 1e-8;
                    % MATRICE DE DIVERGENCE 
                    div_v = dphi(:,2) + phi / rq_reg;

                    for k = 1:3
                        B_e(k, 1:6) = B_e(k, 1:6) - fact * psi(k) * div_v';
                        B_e(k, 7:12) = B_e(k, 7:12) - fact * psi(k) * dphi(:,1)';
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
            norme_solution = norm(U(1:2*numP2));

            tol_abs = 1e-13;  % Tolérance plus stricte
            tol_rel = 1e-12;  % Tolérance plus stricte

            converged = (norme_residu < tol_abs) && ...
                       (norme_delta < tol_rel * max(norme_solution, 1));

%             fprintf('Iter %d: ||res|| = %.2e, ||delta_u|| = %.2e, div·u = %.2e\n', ...
%                     iter, norme_residu, norme_delta,...
%                     norm(B * U(1:2*numP2))); % Norme de la divergence

            if converged
                fprintf('Convergence atteinte en %d itérations\n', iter);
                break;
            end

            iter = iter + 1;
        end
        % Mise à jour pour le pas suivant
        U_old = U;
        if t_current > T_cycle
            idx_cycle = idx_cycle + 1;
            u_r_time(:, idx_cycle) =  U(1:numP2);
            u_z_time(:, idx_cycle) =U(numP2+1:2*numP2);
            p_time(:, idx_cycle) = U(2*numP2+1:end);
        end
        % Calcul du débit effectif
       % Calcul du débit effectif
        u_z_current = u_z_;
        [r_sorted, sort_idx] = sort(r_inlet_original);
        u_z_sorted = u_z_current(bc_z_in(sort_idx));
        Q_effectif = 2*pi * trapz(r_sorted, u_z_sorted .* r_sorted);
        fprintf('Débit effectif: %.2f ml/s (écart: %.1f%%)\n', ...
                Q_effectif*1e6, 100*abs(Q_effectif-Q_instant)/Q_instant);

        % Optionnel: sauvegarde des solutions pour post-traitement
        solutions(i_time,:) = U';  % Matrice de stockage à préallouer
      
        
        %% Vérification de la divergence
        u_r = U(1:numP2);
        u_z = U(numP2+1:2*numP2);
        div_L2 = 0;
        for e = 1:numT
            nodes_e_P2 = tri_P2(e, :);
            coords_P2 = nodes_P2(nodes_e_P2, :);
            z = coords_P2(:, 1);  % Position axiale
            r = coords_P2(:, 2);  % Position radiale

            Bmat = [z(2)-z(1), z(3)-z(1); r(2)-r(1), r(3)-r(1)];
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

                dr_ur = dphi(:,2)' * u_r(nodes_e_P2);
                dz_uz = dphi(:,1)' * u_z(nodes_e_P2);
                ur = phi' * u_r(nodes_e_P2);
                div_v = dr_ur + ur / rq + dz_uz;

                div_integral = div_integral + 2 * pi * rq * wq * div_v^2;
            end
            div_L2 = div_L2 + div_integral;
        end
        div_L2 = sqrt(div_L2);
        fprintf('Norme L2 de la divergence : %.2e\n', div_L2);
        
        u_r = U(1:numP2);
        u_z = U(numP2+1:2*numP2);
        p = U(2*numP2+1:end);
        if ismember(t_current, temp)
            figure;
            trisurf(tri_P2(:,1:3), nodes_P2(:,1), nodes_P2(:,2), u_z);
            title(['Vitesse axiale u_z à t = ', num2str(t_current)]);
            xlabel('z [m]'); ylabel('r [m]'); zlabel('u_z [m/s]');
            colorbar;

            figure;
            trisurf(tri_P2(:,1:3), nodes_P2(:,1), nodes_P2(:,2), u_r);
            title(['Vitesse radiale u_r à t = ', num2str(t_current)]);
            xlabel('z [m]'); ylabel('r [m]'); zlabel('u_r [m/s]');
            colorbar;

            figure;
            trisurf(tri_P1, nodes_P1(:,1), nodes_P1(:,2), p);
            title(['Pression p  à t = ', num2str(t_current)]);
            xlabel('z [m]'); ylabel('r [m]'); zlabel('p [Pa]');
            colorbar;
        end
        [WSS_instant, z_wall, gamma_dot] = compute_WSS(u_z, nodes_P2, tri_P2, ...
                                                 pp, carreau_params, ...
                                                 gradP2, quad_pts, quad_wts);
        if i_time == 1
            z_wall_ref = z_wall;
            WSS_history = zeros(length(t_vec), length(z_wall));
        end

        WSS_history(i_time, :) = WSS_instant';

        fprintf('WSS max: %.2f Pa, min: %.2f Pa\n', max(WSS_instant), min(WSS_instant));

    end
    %% Extraction des solutions finales
    u_r = U(1:numP2);
    u_z = U(numP2+1:2*numP2);
    p = U(2*numP2+1:end);
    %%
    % figure;
    % subplot(2,1,1);
    % Q_vec = arrayfun(Q_target, t_vec) * 1e6;  % Conversion en ml/s
    % plot(t_vec, Q_vec, 'b-', 'LineWidth', 2);
    % grid on; xlabel('Temps (s)'); ylabel('Débit (ml/s)');
    % title('Débit pulsatile imposé');
    % xlim([0 2*T_cycle]);
    % 
    % subplot(2,1,2);
    % % Évolution de la vitesse maximale à l'entrée
    % u_max_vec = zeros(size(t_vec));
    % for i = 1:length(t_vec)
    %     Q_inst = Q_target(t_vec(i));
    %     u_wom = womersley(r_inlet_sorted, Ro, Q_inst, carreau_params);
    %     u_max_vec(i) = max(u_wom);
    % end
    % plot(t_vec, u_max_vec, 'r-', 'LineWidth', 2);
    % grid on; xlabel('Temps (s)'); ylabel('u_{max} (m/s)');
    % title('Vitesse maximale d''entrée (Womersley)');
    % xlim([0 2*T_cycle]);

    %% 4) DIAGNOSTIC FINAL
    fprintf('\n=== RÉSUMÉ SIMULATION ===\n');
    fprintf('Nombre de pas de temps: %d\n', length(t_vec));
    fprintf('Débit moyen: %.1f ml/s\n', Q_mean*1e6);
    fprintf('Débit min/max: %.1f/%.1f ml/s\n', ...
            min(arrayfun(Q_target, t_vec))*1e6, ...
            max(arrayfun(Q_target, t_vec))*1e6);
    fprintf('Nombre de Womersley: α = %.2f\n', Ro*sqrt(omega*rho/mu0));
    
    geometry_params.L = L;
    geometry_params.R_inlet = R_inlet;
    geometry_params.R_max = R_max;
    geometry_params.x_start_aneurysm = x_start_aneurysm;
    geometry_params.x_end_aneurysm = x_end_aneurysm;
    geometry_params.pp = pp;

    simulation_params.nodes_P2 = nodes_P2;
    simulation_params.tri_P2 = tri_P2;
    simulation_params.numP2 = numP2;
    simulation_params.piG = piG;
    simulation_params.gradP2 = gradP2;
    simulation_params.quad_pts = quad_pts;
    simulation_params.quad_wts = quad_wts;
    
    %% Vérification de la divergence
    div_L2 = 0;
    for e = 1:numT
        nodes_e_P2 = tri_P2(e, :);
        coords_P2 = nodes_P2(nodes_e_P2, :);
        z = coords_P2(:, 1);  % Position axiale
        r = coords_P2(:, 2);  % Position radiale

        Bmat = [z(2)-z(1), z(3)-z(1); r(2)-r(1), r(3)-r(1)];
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

            dr_ur = dphi(:,2)' * u_r(nodes_e_P2);
            dz_uz = dphi(:,1)' * u_z(nodes_e_P2);
            ur = phi' * u_r(nodes_e_P2);
            div_v = dr_ur + ur / rq + dz_uz;

            div_integral = div_integral + 2 * pi * rq * wq * div_v^2;
        end
        div_L2 = div_L2 + div_integral;
    end
    div_L2 = sqrt(div_L2);
    % fprintf('Norme L2 de la divergence : %.2e\n', div_L2);

    % --- Calcul divergence élémentaire et visualisation ---
    div_elem = zeros(numT,1);
    max_div_elem = 0;
    for e = 1:numT
        nodes_e_P2 = tri_P2(e, :);
        coords_P2 = nodes_P2(nodes_e_P2, :);
        z = coords_P2(:, 1);  % Position axiale
        r = coords_P2(:, 2);  % Position radiale

        Bmat = [z(2)-z(1), z(3)-z(1); r(2)-r(1), r(3)-r(1)];
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
            dr_ur = dphi(:,2)' * u_r(nodes_e_P2);
            dz_uz = dphi(:,1)' * u_z(nodes_e_P2);
            ur = phi' * u_r(nodes_e_P2);
            div_v = dr_ur + ur / rq + dz_uz;
            div_integral = div_integral + wq * (div_v)^2;
        end
        div_elem(e) = sqrt(div_integral);
        max_div_elem = max(max_div_elem, div_elem(e));
    end

    % figure;
    % patch('Faces', tri_P2(:,1:3), ...
    %       'Vertices', nodes_P2(:,[2 1]), ... % (z,r) pour cohérence
    %       'FaceVertexCData', div_elem, ...
    %       'FaceColor', 'flat', ...
    %       'EdgeColor', 'none');
    % colorbar;
    % axis equal tight;
    % xlabel('z [m]'); ylabel('r [m]');
    % title('Divergence RMS par élément (P2)');
    % Calcul du débit effectif
      %% Visualisation
    figure;
    trisurf(tri_P2(:,1:3), nodes_P2(:,1), nodes_P2(:,2), u_z);
    title('Vitesse axiale u_z');
    xlabel('z [m]'); ylabel('r [m]'); zlabel('u_z [m/s]');
    colorbar;
    
    figure;
    trisurf(tri_P2(:,1:3), nodes_P2(:,1), nodes_P2(:,2), u_r);
    title('Vitesse radiale u_r');
    xlabel('z [m]'); ylabel('r [m]'); zlabel('u_r [m/s]');
    colorbar;
    
    figure;
    trisurf(tri_P1, nodes_P1(:,1), nodes_P1(:,2), p);
    title('Pression p');
    xlabel('z [m]'); ylabel('r [m]'); zlabel('p [Pa]');
    colorbar;
    
   % Calcul de l'OSI
    [OSI, z_osi, TAWSS, RRT] = compute_OSI(WSS_history, z_wall_ref, t_vec);

    % Visualisation
    figure;
    subplot(2,1,1);
    plot(z_osi*1000, TAWSS, 'b-', 'LineWidth', 2);
    xlabel('Position axiale z [mm]');
    ylabel('TAWSS [Pa]');
    title('Time-Averaged Wall Shear Stress');
    grid on;

    subplot(2,1,2);
    plot(z_osi*1000, OSI, 'r-', 'LineWidth', 2);
    xlabel('Position axiale z [mm]');
    ylabel('OSI [-]');
    title('Oscillatory Shear Index');
    ylim([0 0.5]);
    grid on;

    % Carte spatio-temporelle du WSS
    figure;
    imagesc(z_osi*1000, t_vec, WSS_history);
    colorbar;
    xlabel('Position axiale z [mm]');
    ylabel('Temps [s]');
    title('Évolution spatio-temporelle du WSS');
    set(gca, 'YDir', 'normal');
    % Identification des zones à risque
    low_WSS = find(TAWSS < 0.4);  % Faible cisaillement
    high_OSI = find(OSI > 0.2);    % Forte oscillation

    risk_zones = intersect(low_WSS, high_OSI);

    fprintf('\n=== ANALYSE DE RISQUE ===\n');
    fprintf('Zones à faible WSS (<0.4 Pa): %.1f%% du vaisseau\n', ...
            100*length(low_WSS)/length(TAWSS));
    fprintf('Zones à fort OSI (>0.2): %.1f%% du vaisseau\n', ...
            100*length(high_OSI)/length(OSI));
    fprintf('Zones à risque combiné: %.1f%% du vaisseau\n', ...
            100*length(risk_zones)/length(TAWSS));
    % Votre définition initiale
    R_z = @(z) Ro + (R_max - Ro) * sqrt(max(0, 1 - ((z - x_mid_aneurysm )/(L_aneurysm/2)).^2));
    plot_OSI_analysis(z_osi, OSI, TAWSS, RRT, R_z)

    compute_vorticity_streamlines(solutions, simulation_params, geometry_params, t_vec, temp)
end

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
        z = coords_P2(:, 1);  % Position axiale
        r = coords_P2(:, 2);  % Position radiale
        
        % Transformation géométrique
        Bmat = [z(2)-z(1), z(3)-z(1); r(2)-r(1), r(3)-r(1)];
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
            grad_u_r_r = dphi(:,2)' * u_r(nodes_e_P2);  % ?u_r/?r
            grad_u_r_z = dphi(:,1)' * u_r(nodes_e_P2);  % ?u_r/?z
            grad_u_z_r = dphi(:,2)' * u_z(nodes_e_P2);  % ?u_z/?r  
            grad_u_z_z = dphi(:,1)' * u_z(nodes_e_P2);  % ?u_z/?z
            
            % === DÉRIVÉE DE FRÉCHET DU TERME CONVECTIF ===
            % Terme convectif : (u·?)u = u_r ?u/?r + u_z ?u/?z
            % Dérivée : D[(u·?)u][v] = (v·?)u + (u·?)v
            
            dphi_r = dphi(:,2);  % 6×1
            dphi_z = dphi(:,1);  % 6×1
            
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

function u_inlet = womersley(r_nodes, R, Q_target, carreau_params)
    % Calcule le profil de Womersley pour un débit donné Q_target
    % Adapté pour fluide de Carreau
    % INPUTS:
    %   r_nodes    : positions radiales des noeuds d'entrée [m]
    %   R          : rayon du tube [m] 
    %   Q_target   : débit cible à imposer [m^3/s]
    %   carreau_params : paramètres du modèle de Carreau
    % OUTPUT:
    %   u_inlet    : profil de vitesse axiale aux noeuds [m/s]
    
    % Paramètres physiques
    rho = 1050;         % kg/m^3
    mu = carreau_params.eta_inf;  % Viscosité à cisaillement nul
    nu = mu/rho;
    
    % Paramètres temporels (fréquence cardiaque typique)
    f = 1.2;            % Hz - fréquence cardiaque 
    omega = 2*pi*f;     % rad/s
    
    % Nombre de Womersley
    alpha = R * sqrt(omega/nu);
    
    % Coefficient complexe pour Womersley
    i_complex = 1i^(3/2);
    J0_alpha = besselj(0, i_complex * alpha);
    
    % Profil radial complexe φ(r)
    phi_r = 1 - besselj(0, i_complex * alpha * (r_nodes/R)) / J0_alpha;
    phi_real = real(phi_r);  % Partie réelle (même phase que Q)
    
    % Normalisation pour respecter le débit cible
    % Intégration du profil : ∫₀ᴿ φ(r) * r dr
    r_sorted = sort(r_nodes);
    phi_sorted = interp1(r_nodes, phi_real, r_sorted, 'linear', 'extrap');
    
    % Intégrale numérique avec trapz
    Int_phi = 2*pi * trapz(r_sorted, phi_sorted .* r_sorted);
    
    % Amplitude pour respecter Q_target
    if abs(Int_phi) > 1e-12
        A = Q_target / Int_phi;
    else
        % Profil de secours (parabolique) si problème numérique
        warning('Profil Womersley dégénéré, utilisation profil parabolique');
        u_max = 2 * Q_target / (pi * R^2);
        u_inlet = u_max * (1 - (r_nodes/R).^2);
        return;
    end
    
    % Profil final
    u_inlet = A * phi_real;
    
    % Vérification du débit calculé
    Q_check = 2*pi * trapz(r_sorted, interp1(r_nodes, u_inlet, r_sorted) .* r_sorted);
    if abs(Q_check - Q_target) / abs(Q_target) > 0.05
        warning('Écart débit: calculé=%.2e, cible=%.2e', Q_check, Q_target);
    end
end
%% CALCUL DU WSS (Wall Shear Stress)
function [WSS, z_wall, gamma_dot_wall] = compute_WSS(u_z, nodes_P2, tri_P2, pp, carreau_params, gradP2, quad_pts, quad_wts)
    % Calcul robuste du WSS avec diagnostic des anomalies
    
    tol_wall = 5e-4;
    wall_elements = [];
    
    % Identification des éléments pariétaux
    for e = 1:size(tri_P2, 1)
        nodes_e = tri_P2(e, :);
        coords = nodes_P2(nodes_e, :);
        z_nodes = coords(:, 1);
        r_nodes = coords(:, 2);
        R_local = ppval(pp, z_nodes);
        
        if sum(abs(r_nodes - R_local) < tol_wall) >= 2
            wall_elements = [wall_elements; e];
        end
    end
    
    if isempty(wall_elements)
        error('Aucun élément pariétal trouvé !');
    end
    
    n_wall = length(wall_elements);
    z_wall = zeros(n_wall, 1);
    gamma_dot_wall = zeros(n_wall, 1);
    WSS = zeros(n_wall, 1);
    
    % Paramètres Carreau
    mu_inf = carreau_params.eta_inf;
    mu0 = carreau_params.eta_0;
    lambda = carreau_params.lambda;
    n = carreau_params.n;
    a = carreau_params.a;
    m = (n-1)/a;
    
    % Diagnostic
    n_anomalies = 0;
    
    for idx = 1:n_wall
        e = wall_elements(idx);
        nodes_e = tri_P2(e, :);
        coords = nodes_P2(nodes_e, :);
        z = coords(:, 1);
        r = coords(:, 2);
        
        z_elem_center = mean(z);
        R_elem = ppval(pp, z_elem_center);
        
        Bmat = [z(2)-z(1), z(3)-z(1); r(2)-r(1), r(3)-r(1)];
        detB = abs(det(Bmat));
        invB = inv(Bmat);
        
        % ⚠️ AMÉLIORATION : Chercher le point le PLUS PROCHE de la paroi
        % au lieu du plus grand r
        max_r = -inf;
        best_q = 1;
        min_dist = inf;
        
        for q = 1:length(quad_wts)
            xi = quad_pts(q, 1);
            eta = quad_pts(q, 2);
            psi1 = 1 - xi - eta;
            psi2 = xi;
            psi3 = eta;
            phi = [psi1*(2*psi1-1); psi2*(2*psi2-1); psi3*(2*psi3-1);
                   4*psi1*psi2; 4*psi2*psi3; 4*psi3*psi1];
            
            zq = phi' * z;
            rq = phi' * r;
            R_wall_q = ppval(pp, zq);
            
            dist_to_wall = abs(rq - R_wall_q);
            
            % Critère : point le plus proche de la paroi ET r élevé
            if dist_to_wall < min_dist && rq > 0.8*R_wall_q
                min_dist = dist_to_wall;
                best_q = q;
                max_r = rq;
            end
        end
        
        % Calcul au point optimal
        xi = quad_pts(best_q, 1);
        eta = quad_pts(best_q, 2);
        psi1 = 1 - xi - eta;
        psi2 = xi;
        psi3 = eta;
        phi = [psi1*(2*psi1-1); psi2*(2*psi2-1); psi3*(2*psi3-1);
               4*psi1*psi2; 4*psi2*psi3; 4*psi3*psi1];
        
        dphi_ref = gradP2(xi, eta);
        dphi = (invB' * dphi_ref')';
        z_wall(idx) = phi' * z;
        
        % Gradient ∂u_z/∂r (SIGNÉ)
        duz_dr = dphi(:,2)' * u_z(nodes_e);
        
        % ⚠️ DÉTECTION D'ANOMALIES
        if abs(duz_dr) > 1000  % Gradient anormalement élevé
            n_anomalies = n_anomalies + 1;
            if n_anomalies <= 5  % Afficher seulement les 5 premières
                warning('WSS anomalie #%d: z=%.4f m, duz_dr=%.1e s^-1 (élément %d)', ...
                    n_anomalies, z_wall(idx), duz_dr, e);
            end
            
            % ⚠️ CORRECTION : Limiter les gradients aberrants
            duz_dr = sign(duz_dr) * min(abs(duz_dr), 800);
        end
        
        % Magnitude pour la viscosité
        gamma_dot_wall(idx) = abs(duz_dr);
        gamma_reg = max(gamma_dot_wall(idx), 1e-12);
        mu_wall = mu_inf + (mu0 - mu_inf) * (1 + (lambda * gamma_reg)^a)^m;
        
        % ✅ WSS SIGNÉ (garde le sens du flux)
        WSS(idx) = mu_wall * duz_dr;
    end
    
    % Tri par position axiale
    [z_wall, idx_sort] = sort(z_wall);
    gamma_dot_wall = gamma_dot_wall(idx_sort);
    WSS = WSS(idx_sort);
    
    % Diagnostic final
    if n_anomalies > 0
        fprintf('⚠️  %d gradients anormaux détectés et corrigés\n', n_anomalies);
    end
    
    % Statistiques WSS
    WSS_positive = WSS(WSS > 0);
    WSS_negative = WSS(WSS < 0);
    
    fprintf('WSS positif : max=%.2f Pa, moy=%.2f Pa (%d points)\n', ...
        max(WSS_positive), mean(WSS_positive), length(WSS_positive));
    
    if ~isempty(WSS_negative)
        fprintf('WSS négatif : min=%.2f Pa, moy=%.2f Pa (%d points)\n', ...
            min(WSS_negative), mean(WSS_negative), length(WSS_negative));
        fprintf('  → Flux inverse dans %.1f%% de la paroi\n', ...
            100*length(WSS_negative)/length(WSS));
    end
    
    fprintf('Gradient vitesse max: %.1f s^-1\n', max(gamma_dot_wall));
end

%% CALCUL DE L'OSI (Oscillatory Shear Index)
function [OSI, z_osi, TAWSS, RRT] = compute_OSI(WSS_history, z_wall_history, t_vec)
% Calcul correct de l'OSI (Oscillatory Shear Index)
%
% Inputs:
%   WSS_history : matrice (n_time × n_wall_points) - WSS VECTORIEL (avec signe!)
%   z_wall_history : vecteur des positions z à la paroi
%   t_vec : vecteur temps
%
% Outputs:
%   OSI : Oscillatory Shear Index [0-0.5]
%   z_osi : positions z correspondantes
%   TAWSS : Time-Averaged Wall Shear Stress [Pa]
%   RRT : Relative Residence Time [Pa^-1]

n_time = size(WSS_history, 1);
n_wall = size(WSS_history, 2);
T = t_vec(end) - t_vec(1);  % Période du cycle

% ===================================================================
% 1. TAWSS : Time-Averaged Wall Shear Stress
% ===================================================================
% On intègre la MAGNITUDE du WSS (toujours positif)
TAWSS = trapz(t_vec, abs(WSS_history), 1) / T;

% ===================================================================
% 2. Intégrale VECTORIELLE du WSS (garde le signe!)
% ===================================================================
% Ceci capture les inversions de direction
WSS_vector_integrated = trapz(t_vec, WSS_history, 1);  % Pas de abs() ici!

% ===================================================================
% 3. OSI : Oscillatory Shear Index
% ===================================================================
% Formule correcte selon Ku et al. (1985)
% OSI = 0.5 * (1 - |∫WSS_vec dt| / ∫|WSS_vec| dt)

denominator = trapz(t_vec, abs(WSS_history), 1);

% Éviter division par zéro
denominator(denominator < 1e-12) = 1e-12;

OSI = 0.5 * (1 - abs(WSS_vector_integrated) ./ denominator);

% Limiter OSI entre 0 et 0.5 (valeurs théoriques)
OSI = max(0, min(0.5, OSI));

% ===================================================================
% 4. RRT : Relative Residence Time (bonus)
% ===================================================================
% Mesure le temps de résidence des particules
% RRT = 1 / ((1 - 2*OSI) * TAWSS)
% Zones à fort RRT = zones de stagnation

RRT = 1 ./ ((1 - 2*OSI) .* TAWSS + 1e-12);

% Limiter RRT pour éviter les valeurs infinies
RRT = min(RRT, 100);  % Cap à 100 Pa^-1

z_osi = z_wall_history;

% ===================================================================
% Affichage diagnostique
% ===================================================================
fprintf('\n=== ANALYSE OSI ===\n');
fprintf('OSI moyen    : %.3f\n', mean(OSI));
fprintf('OSI max      : %.3f (z = %.4f m)\n', max(OSI), z_osi(OSI == max(OSI)));
fprintf('OSI min      : %.3f\n', min(OSI));
fprintf('TAWSS moyen  : %.3f Pa\n', mean(TAWSS));
fprintf('TAWSS max    : %.3f Pa\n', max(TAWSS));
fprintf('TAWSS min    : %.3f Pa\n', min(TAWSS));

% Zones à risque
zone_low_wss = sum(TAWSS < 0.4) / n_wall * 100;
zone_high_osi = sum(OSI > 0.2) / n_wall * 100;
zone_high_rrt = sum(RRT > 10) / n_wall * 100;

fprintf('\n=== ZONES À RISQUE ===\n');
fprintf('TAWSS < 0.4 Pa   : %.1f%% du vaisseau\n', zone_low_wss);
fprintf('OSI > 0.2        : %.1f%% du vaisseau\n', zone_high_osi);
fprintf('RRT > 10 Pa^-1   : %.1f%% du vaisseau\n', zone_high_rrt);

% Interprétation physique
if mean(OSI) < 0.05
    fprintf('\n⚠️  OSI très faible : écoulement quasi-unidirectionnel\n');
    fprintf('    → Vérifier si les inversions de flux sont capturées\n');
elseif mean(OSI) > 0.3
    fprintf('\n⚠️  OSI élevé : fortes oscillations et recirculations\n');
    fprintf('    → Zone à très haut risque de dégradation pariétale\n');
end

end


% ===================================================================
% FONCTION BONUS : Visualisation OSI
% ===================================================================
function plot_OSI_analysis(z_osi, OSI, TAWSS, RRT, R_z)
% Visualisation complète de l'analyse hémodynamique

figure('Position', [100 100 1200 800]);

% Subplot 1 : OSI le long du vaisseau
subplot(3,1,1);
plot(z_osi, OSI, 'b-', 'LineWidth', 2);
hold on;
yline(0.2, 'r--', 'Seuil risque (OSI=0.2)', 'LineWidth', 1.5);
xlabel('Position axiale z [m]');
ylabel('OSI');
title('Oscillatory Shear Index (OSI)');
grid on;
ylim([0 0.5]);

% Subplot 2 : TAWSS
subplot(3,1,2);
plot(z_osi, TAWSS, 'r-', 'LineWidth', 2);
hold on;
yline(0.4, 'b--', 'Seuil risque (0.4 Pa)', 'LineWidth', 1.5);
xlabel('Position axiale z [m]');
ylabel('TAWSS [Pa]');
title('Time-Averaged Wall Shear Stress');
grid on;

% Subplot 3 : Géométrie avec zones à risque
subplot(3,1,3);
r_plot = R_z(z_osi);
% Colormap selon le risque combiné
risk_score = (OSI > 0.2) + (TAWSS < 0.4);
scatter(z_osi, r_plot, 50, risk_score, 'filled');
hold on;
plot(z_osi, r_plot, 'k-', 'LineWidth', 1);
xlabel('Position axiale z [m]');
ylabel('Rayon [m]');
title('Zones à risque (0=sain, 1=risque modéré, 2=risque élevé)');
colorbar;
colormap(jet);
grid on;
axis equal;

end


function compute_vorticity_streamlines(solutions, simulation_params, geometry_params, t_vec, plot_times)
% Calcule et visualise la vorticité et les lignes de courant à des instants donnés.

% Extraction des paramètres
nodes_P2 = simulation_params.nodes_P2;
tri_P2 = simulation_params.tri_P2;
numP2 = simulation_params.numP2;
piG = simulation_params.piG;
gradP2 = simulation_params.gradP2;
quad_pts = simulation_params.quad_pts;
quad_wts = simulation_params.quad_wts;
L = geometry_params.L;
pp = geometry_params.pp;

% Nombre d'éléments
numT = size(tri_P2, 1);

% Vérification des temps demandés
plot_times = plot_times(ismember(plot_times, t_vec));
if isempty(plot_times)
    error('Aucun temps demandé n''est présent dans t_vec.');
end

% Boucle sur les temps demandés
for t = plot_times
    % Trouver l'indice du temps le plus proche
    [~, idx_time] = min(abs(t_vec - t));
    U = solutions(idx_time, :)';
    
    % Extraction des composantes
    u_r = U(1:numP2);
    u_z = U(numP2+1:2*numP2);
    
    % --- Calcul de la vorticité aux nœuds ---
    vorticity_nodes = zeros(numP2, 1); % Vorticité aux nœuds
    node_count = zeros(numP2, 1);      % Compteur pour la moyenne
    
    for e = 1:numT
        nodes_e_P2 = tri_P2(e, :);
        coords_P2 = nodes_P2(nodes_e_P2, :);
        z = coords_P2(:, 1);  % Position axiale
        r = coords_P2(:, 2);  % Position radiale
        
        % Matrice de transformation
        Bmat = [z(2)-z(1), z(3)-z(1); r(2)-r(1), r(3)-r(1)];
        detB = abs(det(Bmat));
        invB = inv(Bmat);
        
        % Calcul de la vorticité à chaque nœud de l'élément
        for i = 1:6
            node_idx = nodes_e_P2(i);
            
            % Coordonnées locales du nœud i
            if i == 1
                xi = 0; eta = 0;
            elseif i == 2
                xi = 1; eta = 0;
            elseif i == 3
                xi = 0; eta = 1;
            elseif i == 4
                xi = 0.5; eta = 0;
            elseif i == 5
                xi = 0.5; eta = 0.5;
            elseif i == 6
                xi = 0; eta = 0.5;
            end
            
            % Gradients au nœud
            dphi_ref = gradP2(xi, eta);
            dphi = (invB' * dphi_ref')';
            
            % Gradients des vitesses
            duz_dr = dphi(:,2)' * u_z(nodes_e_P2); % ∂u_z/∂r
            dur_dz = dphi(:,1)' * u_r(nodes_e_P2); % ∂u_r/∂z
            
            % Vorticité (composante θ en axisymétrique)
            omega_theta = duz_dr - dur_dz;
            
            % Accumulation pour la moyenne
            vorticity_nodes(node_idx) = vorticity_nodes(node_idx) + omega_theta;
            node_count(node_idx) = node_count(node_idx) + 1;
        end
    end
    
    % Moyenne de la vorticité aux nœuds
    vorticity_nodes = vorticity_nodes ./ max(node_count, 1);
    
    % --- Calcul de la fonction de courant ---
    psi = zeros(numP2, 1); % Fonction de courant aux nœuds P2
    M_psi = sparse(numP2, numP2); % Matrice de masse
    f_psi = zeros(numP2, 1); % Second membre
    
    for e = 1:numT
        nodes_e_P2 = tri_P2(e, :);
        coords_P2 = nodes_P2(nodes_e_P2, :);
        z = coords_P2(:, 1);
        r = coords_P2(:, 2);
        
        % Matrice de transformation
        Bmat = [z(2)-z(1), z(3)-z(1); r(2)-r(1), r(3)-r(1)];
        detB = abs(det(Bmat));
        invB = inv(Bmat);
        
        % Matrices locales
        M_e = zeros(6, 6);
        f_e = zeros(6, 1);
        
        for q = 1:length(quad_wts)
            xi = quad_pts(q, 1);
            eta = quad_pts(q, 2);
            wq = quad_wts(q) * detB;
            
            phi = piG(xi, eta);
            dphi_ref = gradP2(xi, eta);
            dphi = (invB' * dphi_ref')';
            
            % Coordonnée radiale
            rq = max(sum(r .* phi), 1e-12);
            
            % Interpolation des vitesses
            ur_q = phi' * u_r(nodes_e_P2);
            uz_q = phi' * u_z(nodes_e_P2);
            
            % Contribution à la matrice de masse
            M_e = M_e + 2 * pi * rq * wq * (phi * phi');
            
            % Second membre
            f_e = f_e + 2 * pi * wq * (rq * uz_q * dphi(:,2) - rq * ur_q * dphi(:,1));
        end
        
        % Assemblage global
        M_psi(nodes_e_P2, nodes_e_P2) = M_psi(nodes_e_P2, nodes_e_P2) + M_e;
        f_psi(nodes_e_P2) = f_psi(nodes_e_P2) + f_e;
    end
    
    % Condition aux limites : ψ = 0 sur l'axe (r=0)
    tol_bc = 1e-4;
    bc_axis = find(abs(nodes_P2(:,2)) < tol_bc);
    for idx = bc_axis'
        M_psi(idx, :) = 0;
        M_psi(idx, idx) = 1;
        f_psi(idx) = 0;
    end
    
    % Résolution du système pour ψ
    psi = M_psi \ f_psi;
    
    % --- Visualisation ---
    figure('Position', [100, 100, 1200, 400]);
    
    % Vorticité
    subplot(1, 2, 1);
    trisurf(tri_P2(:,1:3), nodes_P2(:,1), nodes_P2(:,2), vorticity_nodes);
    title(['Vorticité \omega_\theta à t = ', num2str(t), ' s']);
    xlabel('z [m]'); ylabel('r [m]'); zlabel('\omega_\theta [s^{-1}]');
    colorbar;
    shading interp;
    axis equal tight;
    view(2);
    
    % Lignes de courant
    subplot(1, 2, 2);
    % Interpolation sur une grille régulière
    z_grid = linspace(0, L, 100);
    r_max = max(ppval(pp, z_grid));
    r_grid = linspace(0, r_max, 50);
    [Z, R] = meshgrid(z_grid, r_grid);
    
    % Interpolation de ψ sur la grille
    PSI_grid = griddata(nodes_P2(:,1), nodes_P2(:,2), psi, Z, R);
    
    % Masquer les points hors du domaine
    R_wall = ppval(pp, Z);
    mask = R > R_wall;
    PSI_grid(mask) = NaN;
    
    % Trace des lignes de courant
    contour(Z, R, PSI_grid, 20, 'LineWidth', 1.5);
    title(['Lignes de courant à t = ', num2str(t), ' s']);
    xlabel('z [m]'); ylabel('r [m]');
    colorbar;
    axis equal tight;
    hold on;
    % Trace de la paroi
    plot(z_grid, ppval(pp, z_grid), 'k-', 'LineWidth', 2);
    plot(z_grid, zeros(size(z_grid)), 'k-', 'LineWidth', 2);
    hold off;
end

end