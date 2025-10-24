function [St,Ma]=AssembleGlobalMatrices(p,t,ord)
N=length(p); % No of nodes
M=length(t); % No of elements
St=sparse(N,N);
Ma=sparse(N,N);

%% Reference Quadrature Rule and Basis Functions
[iw,ip]=RefElemQuad(ord);
[basis_ip,grad_basis_ip]=basis_linear_2D(ip);

Local_M_matr=zeros(3,3);
Local_S_matr=zeros(3,3);

%% Loop Over All Triangles
for k=1:M
    nodes=p(:,t(1:3,k));
    C=[nodes(:,2)-nodes(:,1), nodes(:,3)-nodes(:,1)];
    cv=nodes(:,1);
    Cdet=C(1,1)*C(2,2)-C(1,2)*C(2,1);
    CinvT=(1/Cdet)*[C(2,2),-C(2,1);-C(1,2),C(1,1)];

    node_ids = t(1:3,k);
    % Compute Local Matrices
    for i=1:3
        for j=1:3
            Local_M_matr(i,j) = (basis_ip(i,:) .* basis_ip(j,:)) * iw' * Cdet;
            Local_S_matr(i,j)=(CinvT*grad_basis_ip(:,1,i))'*...
                               (CinvT*grad_basis_ip(:,1,j))*Cdet;
        end
    end

    % Assemble into global matrices
    St(node_ids, node_ids) = St(node_ids, node_ids) + Local_S_matr;
    Ma(node_ids, node_ids) = Ma(node_ids, node_ids) + Local_M_matr;
end
end