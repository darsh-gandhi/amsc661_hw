clear; close all

%% L-shape

clear; close all

initmsh();
node = [0 0; 2 0; 2 2; 1 2; 1 3; 0 3];
edge = [1 2; 2 3; 3 4; 4 5; 5 6; 6 1];
[vert,etri,tria,tnum] = refine2(node, edge, [], [], 0.2);
patch('faces',tria(:,1:3),'vertices',vert,'facecolor',[0.8 0.9 1],'edgecolor',[.2 .2 .2]);
axis image off;

%% pentagon

clear; close all

initmsh();
t5 = linspace(pi/2, pi/2+2*pi, 6);
outer = [cos(t5(1:5))', sin(t5(1:5))'];
inner = 0.5*outer;
node = [outer; inner];
edge = [1 2; 2 3; 3 4; 4 5; 5 1; 6 7; 7 8; 8 9; 9 10; 10 6];
[vert,etri,tria,tnum] = refine2(node, edge, [], [], 0.1);
patch('faces',tria(:,1:3),'vertices',vert,'facecolor',[0.8 0.9 1],'edgecolor',[.2 .2 .2]);
axis image off;

%% semi-circle with holes

clear; close all

initmsh();
th = linspace(pi,2*pi,40)';
arc = [cos(th), sin(th)];
n = size(arc,1);
tc = linspace(0,2*pi,25)';
tc = tc(1:end-1);
c1 = [0.45+0.2*cos(tc), -0.45+0.2*sin(tc)];
c2 = [-0.45+0.2*cos(tc), -0.45+0.2*sin(tc)];
nc = size(c1,1);
node = [arc; c1; c2];
e_arc = [(1:n-1)' (2:n)'; n 1];
e_c1  = n + [(1:nc-1)' (2:nc)'; nc 1];
e_c2  = n + nc + [(1:nc-1)' (2:nc)'; nc 1];
edge = [e_arc; e_c1; e_c2];
[vert,etri,tria,tnum] = refine2(node, edge, [], [], 0.08);
patch('faces',tria(:,1:3),'vertices',vert,'facecolor',[0.8 0.9 1],'edgecolor',[.2 .2 .2]);
axis image off;
