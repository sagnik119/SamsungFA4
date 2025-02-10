function [FEAS_FLAG, bu_a, info] = admMAC_rate_region(H, Lxu, bu, Eu)
%admMAC_rate_region Determine whether the target rate vector bu can be
%achieved under channel H and power constraint Eu via rate region
%   Input arguments:
%       - H: Ly-by-Lx-by-N channel matrix. H(:,:,n) denotes the channel for
%           the n-th tone.
%       - Lxu: number of transmit antennas of each user. It can be either a
%           scalar or a length-U vector. If it is a scalar, every user has
%           Lxu transmit antennas; otherwise user u has Lxu(u) transmit
%           antennas.
%       - bu: target rate of each user, length-U vector.
%       - Eu: Power constraint on each user, length-U vector.
%   Outputs:
%       - FEAS_FLAG: indicator of achievability. FEAS_FLAG=0 if the target
%           is not achievable; FEAS_FLAG=1 if the target is achievable by a
%           single ordering; FEAS_FLAG=2 if the target is achievable by
%           time-sharing
%       - bu_a: U-by-1 vector showing achieved sum rate of each user. 
%       - info: various length output depending on FEAS_FLAG
%           --if FEAS_FLAG=0: empty
%           --if FEAS_FLAG=1: 1-by-5 cell array containing
%               {Rxxs, Eun, bun, theta, w} corresponds to the single vertex
%           --if FEAS_FLAG=2: v-by-6 cell array, with each row representing
%               a time-shared vertex {fraction, Rxxs, Eun, bun, theta, w}
%
%       - Rxxs: U-by-N cell array containing Rxx(u,n)'s if Lxu is a
%           length-U vector; or Lxu-by-Lxu-by-U-by-N tensor if Lxu is a
%           scalar. If the rate target is infeasible, output 0.
%       - Eun:  U-by-N matrix showing users' transmit power on each tone.
%           If infeasible, output 0.
%       - bun: U-by-N matrix showing users' rate on each tone. If
%           infeasible, output 0.
%       - theta: U-by-1 Lagrangian multiplier w.r.t. target rates
%       - w: U-by-1 Lagrangian multiplier w.r.t. power constraints

[Ly, ~, N] = size(H);
if N == 1
    H = reshape(H,Ly,[],1);
end
Eu = reshape(Eu,[],1);
bu = reshape(bu,[],1);
U = length(Eu);
theta = ones(U,1);
A = eye(U)*U;
count = 1
% error tolerance
err = 1e-6;
% tolerance of convex hull boundary
conv_tol = 1e-6;


if max(Lxu) == 1
    SCALAR_FLAG = 1;
    [Eun, w, bun] = maxRMAC_cvx(H, Eu, theta);
    Rxxs = Eun;
else
    SCALAR_FLAG = 0;
    [Rxxs, Eun, w, bun] = maxRMAC_vector_cvx(H, Lxu, Eu, theta);
end

% achieved users' rates
bu_v = sum(bun,2)';
g = bu_v' - bu;

if theta'*g < -err % infeasible criterion
    FEAS_FLAG = 0;
    bu_a = bu_v';
    info=[];
    return
end

if sqrt(g'*A*g) <= err || min(g) >= 0 % achievable by current vertex
    FEAS_FLAG = 1;
    bu_a = bu_v';
    info = table(bu_v, {Rxxs}, {Eun}, {bun}, {theta}, {w}); % detailed info of boundary vertices
	info.Properties.VariableNames(2:end) = {'Rxxs' 'Eun' 'bun' 'theta' 'w'};
    return
end

% add the first vertex
known_vertices = table(bu_v, {Rxxs}, {Eun}, {bun}, {theta}, {w}); % detailed info of boundary vertices
known_vertices.Properties.VariableNames(2:end) = {'Rxxs' 'Eun' 'bun' 'theta' 'w'};
bd_vertices = bu_v; % track critical boundary vertices
bd_V = 1; % track number of critical boundary vertices
vertices = de2bi(0:2^U-1).*bu_v; % all boundary vertices
tess = convhulln(vertices);

while 1
    count = count + 1
    
    % update ellipsoid
    tmp = A*g/sqrt(g'*A*g);
    theta = theta - 1/(U + 1)*tmp;
    A = U^2/(U^2 - 1) * (A - 2/(U + 1)*(tmp*tmp'));
    ind = find(theta < zeros(U,1));
    while ~isempty(ind)
        g = zeros(U,1);
        g(ind(1)) = -1;
        tmp = A*g/sqrt(g'*A*g);
        theta = theta - 1/(U + 1)*tmp;
        A = U^2/(U^2 - 1) * (A - 2/(U + 1)*(tmp*tmp'));
        ind = find(theta < zeros(U,1));
    end   
    
    % get the next vertex
    if SCALAR_FLAG
        [Eun, w, bun] = maxRMAC_cvx(H, Eu, theta);
        Rxxs = Eun;
    else
        [Rxxs, Eun, w, bun] = maxRMAC_vector_cvx(H, Lxu, Eu, theta);
    end
    % achieved users' rates
    bu_v = sum(bun,2)';
    g = bu_v' - bu;
    
    if theta'*g < -err % infeasible criteria
        FEAS_FLAG = 0;
        bu_a = bu_v';
        info = {};
        break
    end
    
    if sqrt(g'*A*g) <= err || min(g) >= 0 % achievable by current vertex
        FEAS_FLAG = 1;
        bu_a = bu_v';
        info = {Rxxs, Eun, bun, theta, w};
        break
    end
    
    if ~inhull(bu_v, vertices, tess, conv_tol)
        % add new vertices and build new convex hull
        known_vertices = [known_vertices; {bu_v, {Rxxs}, {Eun}, {bun}, {theta}, {w}}];
        bd_vertices_extend = [bd_vertices; bu_v];
        bd_V = bd_V + 1;
        new_vs = de2bi(0:2^U-1).*bu_v;
        vertices = [vertices;new_vs];
        tess = convhulln(vertices);
        % keep only the vertices at the boundary of convex hull
        vertices = vertices(unique(tess),:);
        bd_vertices = intersect(bd_vertices_extend, vertices, 'rows');
        tess = convhulln(vertices);
        % delete inner vertices
        if size(bd_vertices, 1) < bd_V
            to_remove = setdiff(bd_vertices_extend, bd_vertices, 'rows');
            known_vertices(ismember(known_vertices.bu_v, to_remove, 'rows'),:) = [];
            bd_V = size(bd_vertices, 1);
        end
        if inhull(bu', vertices, tess, conv_tol) % achievable by time-share
            FEAS_FLAG = 2;
            cvx_begin quiet  % compute time-share combination
                variable frac(bd_V) nonnegative
                minimize norm(frac,1)
                subject to
                    sum(frac) <= 1;
                    bd_vertices'*frac == bu;
            cvx_end
            active_idx = find(frac>=err);
            a_v = bd_vertices(active_idx,:);
            a_frac = frac(active_idx);
            [tab_pos,frac_loc] = ismember(known_vertices.bu_v, a_v,'rows');
            info = known_vertices(tab_pos,:);
            info.frac = a_frac(frac_loc);
            bu_a = bu;
            break
        end
    end
end

