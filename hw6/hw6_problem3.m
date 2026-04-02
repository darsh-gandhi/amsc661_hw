close all; clear;

%note <- need to use distmesh2d package with this code to generate mesh for
%the plate

%% base params & problem geometry
L = 3;
cx = 1.5; cy = 1.5; % circle center
r = 1; % radius
h0 = 0.12;

% defining geometry of the plate
N_circ = round(2*pi*r/h0);
t_circ = linspace(0,2*pi*(N_circ-1)/N_circ, N_circ)';
circle_pts = [cx + r*cos(t_circ), cy + r*sin(t_circ)];
corners = [0,0; L,0; L,L; 0,L];
pfix = [corners; circle_pts];

% dist function
fd = @(p) drectangle(p,0,L,0,L);
fh = @(p) huniform(p);

[pts, tri] = distmesh2d(fd, fh, h0, [0,0; L,L], pfix); %distmesh for mesh generation

Npts = size(pts,1);
Ntri = size(tri,1);

%helper function
function M = genStiffMat(verts,a)
    G = [ones(1,3); verts']\[zeros(1,2); eye(2)];
    M = a*0.5*abs(det([ones(1,3); verts']))*2*(G*G');
end

%% fem assembly, solution, & plotting

tol = 1e-10;
left_idx = find(abs(pts(:,1)) < tol);
right_idx = find(abs(pts(:,1) - L) < tol);
dirichlet = [left_idx; right_idx];
FreeNodes = setdiff(1:Npts, dirichlet')';

centroids = (pts(tri(:,1),:) + pts(tri(:,2),:) + pts(tri(:,3),:)) / 3;

cases = {1.2, 1.0, '(a) a_1=1.2, a_2=1'; 0.8, 1.0, '(b) a_1=0.8, a_2=1'};

fig = figure(1);
theta_plot = linspace(0, 2*pi, 300);
cx_c = cx + r*cos(theta_plot);
cy_c = cy + r*sin(theta_plot);

for row = 1:2
    a1 = cases{row,1};
    a2 = cases{row,2};
    label = cases{row,3};
    inside = (centroids(:,1)-cx).^2 + (centroids(:,2)-cy).^2 <= r^2;
    a_tri = a2*ones(Ntri,1);
    a_tri(inside) = a1;

    %assembly stiffness matrix
    A = sparse(Npts, Npts);
    for j = 1:Ntri
        verts = pts(tri(j,:),:);
        M = genStiffMat(verts, a_tri(j));
        A(tri(j,:),tri(j,:)) = A(tri(j,:),tri(j,:)) + M;
    end

    %boundary
    u = sparse(Npts, 1);
    u(left_idx) = 0;
    u(right_idx) = 1;

    b = -A*u;

    %solve
    u(FreeNodes) = A(FreeNodes, FreeNodes)\b(FreeNodes);
    u = full(u);

    grad_u = zeros(Ntri, 2);
    current_mag = zeros(Ntri, 1);
    for j=1:Ntri
        verts = pts(tri(j,:),:);
        G = [ones(1,3); verts']\[zeros(1,2); eye(2)];
        grad_u(j,:) = G'*u(tri(j,:));
        j_vec = -a_tri(j)*grad_u(j,:);
        current_mag(j) = norm(j_vec);
    end

    abs_current_verts = zeros(Npts,1);
    count_tri= zeros(Npts,1);
    for j=1:Ntri
        abs_current_verts(tri(j,:)) = abs_current_verts(tri(j,:)) + current_mag(j);
        count_tri(tri(j,:)) = count_tri(tri(j,:)) + 1;
    end
    abs_current_verts = abs_current_verts./count_tri;

    % voltage
    subplot(2,2,2*(row-1)+1);
    trisurf(tri, pts(:,1), pts(:,2), u, 'EdgeColor','none', 'FaceColor','interp');
    view(2); axis equal tight;
    colorbar;
    hold on;
    plot3(cx_c, cy_c, ones(size(cx_c)), 'k--', 'LineWidth', 1.5);
    title(['Voltage ' label], 'FontSize', 11);
    xlabel('x'); ylabel('y');
    clim([0 1]);

    % density
    subplot(2,2,2*(row-1)+2);
    trisurf(tri, pts(:,1), pts(:,2), abs_current_verts, 'EdgeColor','none', 'FaceColor','interp');
    view(2); axis equal tight;
    colorbar;
    hold on;
    plot3(cx_c, cy_c, max(abs_current_verts)*ones(size(cx_c)), 'w--', 'LineWidth', 1.5);
    title(['Density ' label], 'FontSize', 11);
    xlabel('x'); ylabel('y');
end