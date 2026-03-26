clear; close all

%% L-shape
pv = [0 0; 0 3; 1 3; 1 2; 2 2; 2 0; 0 0];

fd = @(p) dpoly(p, pv);
[p,t] = distmesh2d(fd, @huniform, 0.2, [0,0; 2,3], pv);

%% pentagon w/ pentagonal hole in center and inscribed star
clear; close all

t5 = linspace(pi/2, pi/2+2*pi, 6);
pv = [cos(t5)', sin(t5)'];

fd = @(p) dpoly(p, pv);
[p,t] = distmesh2d(fd, @huniform, 0.1, [-1,-1; 1,1], pv);

%% semicircle w/ two circular holes
clear; close all

t5 = linspace(pi/2, pi/2+2*pi, 6);
pv_outer = [cos(t5)', sin(t5)'];
pv_inner = [0.4*cos(t5)', 0.4*sin(t5)'];  % same shape, scaled down

fd = @(p) ddiff(dpoly(p, pv_outer), dpoly(p, pv_inner));
[p, t] = distmesh2d(fd, @huniform, 0.1, [-1,-1;1,1], [pv_outer(1:5,:); pv_inner(1:5,:)]);
fd = @(p) ddiff(dpoly(p, pv_outer), dpoly(p, pv_inner));
[p, t] = distmesh2d(fd, @huniform, 0.1, [-1,-1;1,1], pv_outer(1:5,:));