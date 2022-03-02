function [Doms] = subDoms(NMatSDD,CMatSDD,Dvvs,subInd,subOL,F,uSoln)
% uSoln should include exact boundary data

% Initialize full adjacency matrix, including the boundary nodes
% [N,m]= size(NMatSDD);
M = max(max(NMatSDD));
% adj = zeros(N,M); 

% number of sub indices
K = numel(subInd);

% % Fill in the adjacency matrix using the neighbor matrix
% for i = 1:N
%     for j = 1:m
%         adj(i,NMatSDD(i,j))=1;
%     end
% end
% 
% % Adjacency matrix only with interior nodes
% Adj = adj(1:N,1:N);

% Build the local subdomains
d = length(Dvvs);
Doms = [];
for k = 1:K

    % global and local index arrays
    % Want all index arrays to be column vectors 
    
    % Subdomain with no (exterior) overlap
    nOL_g = subInd{k};
    nOL_l = (1:length(nOL_g))';

    % Only exterior overlap
    OL_g = subOL{k};
    OL_l = (1:length(OL_g))'+length(nOL_g);

    % Interior nodes
    I_g = [nOL_g; OL_g];
    I_l = [nOL_l;OL_l];
    
    % Boundary Nodes
    B_g = setdiff(NMatSDD(I_g,:),I_g);
    B_l = ((1:length(B_g))+length(I_g))';
    
    % Local to global and global to local maps
    L2G = [I_g;B_g];
    G2L = zeros(M,1);
    G2L(L2G) = [I_l;B_l];

    % Local NMat and CMat
    NMat_l = G2L(NMatSDD(I_g,:));
    CMat_l = CMatSDD(I_g,:);

    % Local Dvvs
    Dvvs_l = cell(size(Dvvs));
    for j = 1:d
        Dvv = Dvvs{j};
        Dvvs_l{j} = Dvv(I_g,L2G);
    end
    % Local F, Interior Soln, and boundary
    F_l = F(I_g);
    uInt_l = uSoln(I_g);
    uBdry_l = uSoln(B_g);


    % Build domain object
    Dom.nOL_g = nOL_g;
    Dom.nOL_l = nOL_l;
    Dom.OL_g = OL_g;
    Dom.OL_l = OL_l;
    Dom.I_g = I_g;
    Dom.I_l = I_l;
    Dom.B_g = B_g;
    Dom.B_l = B_l;
    Dom.L2G = L2G;
    Dom.G2L = G2L;
    Dom.NMat_l = NMat_l;
    Dom.CMat_l = CMat_l;
    Dom.Dvvs_l = Dvvs_l;
    Dom.F_l = F_l;
    Dom.uInt_l = uInt_l;
    Dom.uBdry_l = uBdry_l;

    Doms = [Doms Dom];
end


end