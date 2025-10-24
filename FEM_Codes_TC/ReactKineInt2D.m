function F = ReactKineInt2D(p1, t1, p2, t2, u1, u2, eta, f, ord)
% ReactKineInt2D computes ∫ f(u1(x,y),u2(x,y),eta) * φ_i(x,y) dxdy over 2D domain
%
% Inputs:
%   p     - 2 x N node coordinates
%   t     - 3 x Ne element connectivity (triangles)
%   u1    - N x 1 nodal values (field 1)
%   u2    - N x 1 nodal values (field 2)
%   eta   - parameter(s) passed to f
%   f     - function handle f(u1,u2,eta)
%   ord   - quadrature order
%
% Output:
%   F     - N x 1 global load vector

N  = size(p1,2);        % number of nodes
Ne = size(t1,2);        % number of elements
F  = zeros(N,1);       % global load vector

% Quadrature rule on reference triangle (xi,eta)
[iw, ip] = RefElemQuad(ord);         % ip: nq x 2, iw: nq x 1
[basis_ip, ~] = basis_linear_2D(ip); % 3 x nq

% Loop over triangles
for k = 1:Ne
    nodes1 = t1(:,k);             % indices of nodes in this element
    coords1 = p1(:,nodes1);        % 2x3 coordinates of triangle
    nodes2 = t2(:,k);
    coords2 = p2(:,nodes2);

    % Jacobian mapping from reference triangle
    J1 = [coords1(:,2)-coords1(:,1), coords1(:,3)-coords1(:,1)]; % 2x2
    detJ1 = abs(det(J1));
    J2 = [coords2(:,2)-coords2(:,1), coords2(:,3)-coords2(:,1)]; % 2x2
    detJ2 = abs(det(J2));

    % Local nodal values
    u1_local = u1(nodes1);
    u2_local = u2(nodes2);

    % Interpolate u1,u2 at quadrature points
    u1_ip = u1_local' * basis_ip;   % 1 x nq
    u2_ip = u2_local' * basis_ip;

    % Reaction evaluation
    f_vals = f(u1_ip, u2_ip, eta);  % 1 x nq

    if any(isnan(f_vals))
        error('f returned NaN at element %d', k);
    end

    % Local load vector
    local_load = zeros(3,1);
    for i = 1:3
        local_load(i) = sum(f_vals .* basis_ip(i,:) .* iw) * detJ1;
    end

    % Assemble into global load vector
    F(nodes1) = F(nodes1) + local_load;
end
end
