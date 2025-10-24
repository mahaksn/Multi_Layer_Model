function F = ReactKineInt1D(p, edges, u1, u2, eta, f, ord)
% ReactKineInt1D computes ∫ f(u(x)) * φ_i(x) dx over 1D domain
% Inputs:
%   p     - 1 x N node coordinates
%   edges - Ne x 2 connectivity of nodes (each row = [i j])
%   u     - N x 1 nodal values
%   f     - function handle for reaction kinetics (e.g., f = @(u) a - u)
%   ord   - quadrature order

N = length(p);          % Number of nodes
Ne = size(edges, 1);    % Number of elements (edges)
F = zeros(N, 1);        % Global load vector

% Quadrature rule and basis functions
[iw, ip] = RefEdgeQuad(ord);             % weights and points on [0,1]
[basis_ip,~] = basis_linear_1D(ip);          % 2 x n_quad_pts

% Loop over each element (edge)
for k = 1:Ne
    node_ids = edges(k, :);             % Node indices for this element
    x = p(node_ids);                    % Coordinates of the two nodes
    h = abs(x(2) - x(1));               % Element length

    u1_local = u1(node_ids);             % Local nodal values
    u2_local = u2(node_ids);
    % u3_local = u3(node_ids);
    u1_ip = u1_local' * basis_ip;         % u at quadrature points (1 x n_quad)
    u2_ip = u2_local' * basis_ip;
    % u3_ip = u3_local' * basis_ip;

    f_vals = f(u1_ip,u2_ip,eta);
    f_vals(1);
if any(isnan(u1_ip)) || any(isnan(u2_ip)) %|| any(isnan(u3_ip))
    error('NaN detected in u_ip values at element %d', k);
end

if any(isnan(f_vals))
    disp('Inputs to f:');
    disp([u1_ip; u2_ip]);
    error('f returned NaN at element %d', k);
end
    % Compute local load vector
    local_load = zeros(2, 1);
    for i = 1:2
        local_load(i) = sum(f_vals .* basis_ip(i,:) .* iw) * h;
    end

    % Assemble into global load vector
    F(node_ids) = F(node_ids) + local_load;
end
end
