% generate random network on periodic plane
% written by Thomas Fai (tfai@seas.harvard.edu)
% and Steven Petteruti, September 2016

clear all
format long
close all

load('/Users/Steven/Downloads/erData/myvertices.mat');
Edges = xlsread('/Users/Steven/Downloads/erData/myedges1.xls');
Vertices = Expression1;
max_x = max(Vertices(:,1))
max_y = max(Vertices(:,2))
max_value = max([max_x max_y])
Vertices = Vertices/max_value;

numb_points = 100;
N = numb_points;
R = 1;

%%%%%%%%%%%%%%%%%%%%%% using random points in plane

% regular triangulation
% ruler=linspace(0,1-1/sqrt(N),sqrt(N));
% 
% V=repmat(ruler,1,sqrt(N));
% U=reshape(V,sqrt(N),sqrt(N)); x=reshape(U',1,N);
% y=V;

% random graph
x = rand(1,N);
y = rand(1,N);
%x = Expression2(:,1)';
%y = Expression2(:,2)';
max_val = max([x y]);
x = x/max_val;
y = y/max_val;
X=[x' y'];
dt=DelaunayTri(x',y');
 
 % Make graph on periodic plane
 % Place periodic copies
 X=[x' y'; x'-1 y'; x'+1 y';
    x'-1 y'-1; x' y'-1; x'+1 y'-1;
    x'-1 y'+1; x' y'+1; x'+1 y'+1];

%perform Delaunay triangulation
tr=DelaunayTri(X);
[szT scratch]=size(tr);
TV=tr(1:szT,:);

% Build dual graph

%get dual graph
[szXX scratch]=size(X);
[A,X_dual] = compute_dual_graph(TV,[X zeros(szXX,1)]');
X_dual = X_dual(1:2,:)';
[szXX_dual scratch]=size(X_dual);
EV_dual = zeros(szXX_dual,3); %assumes 3 edges per vertex
e_dual = [];
for i=1:szXX_dual
    inds = find(A(i,:) == 1);
    EV_dual(i,1:numel(inds)) = inds;
end

%find dual graph edges in unit square
X11 = (X_dual(:,1) <= 1 & X_dual(:,1) >= 0);
X12 = (X_dual(:,2) <= 1 & X_dual(:,2) >= 0);
indd=logical(X11.*X12);
inds = find(indd == 1);

N_dual = sum(indd);
indd2 = zeros(numel(indd),1);
indd2(indd) = 1;
for i=1:N_dual
   indd2(EV_dual(inds(i),:)) = 1; 
end

%identify dual graph ghost points
iXst_dual = sum(indd2)-N_dual;
indd2 = logical(indd2);
% X_dual = X_dual(indd2,:);

inds_ghost = find(indd2 == 1 & indd == 0);

%rearrange so that unit square points are first
X_dual = [X_dual(inds,:); X_dual(inds_ghost,:)];
EV_dual_tmp = EV_dual;
for i=1:N_dual
    EV_dual(EV_dual_tmp == inds(i)) = i;
end
for i=1:iXst_dual
    EV_dual(EV_dual_tmp == inds_ghost(i)) = N_dual+i;
end

EV_dual_tmp = EV_dual;
for i=1:N_dual
    EV_dual(i,:) = EV_dual_tmp(inds(i),:);
end
for i=1:iXst_dual
    EV_dual(N_dual+i,:) = EV_dual_tmp(inds_ghost(i),:);
end
EV_dual = EV_dual(1:(N_dual+iXst_dual),:);

% e_dual = zeros(3*N_dual/2,2);
% eind = 1;
e_dual = [];
for i=1:N_dual
    for j=1:3
        if (EV_dual(i,j) > i)
            e_dual = [e_dual; i EV_dual(i,j)];
%             e_dual(eind,:) = [i,EV_dual(i,j)];
%             eind = eind+1;
        end
    end
end

%find indices of ghost points periodic images in unit square
X_unit = X_dual(1:N_dual,:);
pim_dual = zeros(iXst_dual,1);
pxind_dual = zeros(numel(iXst_dual),1);
pyind_dual = zeros(numel(iXst_dual),1);

% find coordinates of ghost points in (3 by 3) supercell
% (ghost point position is given by
% X(ghost,:) = X(pim(ghost),:)+L*[pxind(ghost); pyind(ghost)];)

tol = 1e-9;
for i=1:iXst_dual
    inds = find(abs(X_unit(:,1) - mod(X_dual(N_dual+i,1),1)) < tol & ...
        abs(X_unit(:,2) - mod(X_dual(N_dual+i,2),1)) < tol);
    if (numel(inds) ~= 1)
        'couldnt find periodic image'
        i
        inds
    end
    pim_dual(i) = inds;
    if (X_dual(N_dual+i,1) < 0)
        pxind_dual(i) = -1;
    elseif (X_dual(N_dual+i,1) > 1)
        pxind_dual(i) = 1;
        
    else
        pxind_dual(i) = 0;
    end
    if (X_dual(N_dual+i,2) < 0)
        pyind_dual(i) = -1;
        
    elseif (X_dual(N_dual+i,2) > 1)
        pyind_dual(i) = 1;
        
    else
        pyind_dual(i) = 0;
    end
    
end

X_dual = R*X_dual;
Xplt_dual = X_dual;
[szXX_dual scratch]=size(Xplt_dual);

% Continue with original graph

% Find vertices in the unit square
X11=(TV(:,1)>N);
X12=(TV(:,2)>N);
X13=(TV(:,3)>N);
indd=logical(1-X11.*X12.*X13);
TV=TV(indd,:);

% find number of triangles in the unit square
[szT scratch]=size(TV);

% identify ghost points
iXst = TV(TV>N);
Xplt = [X(1:N,:); X(iXst,:)];
szet = numel(iXst);

%find indices of ghost points periodic images in unit square
pim = mod(iXst-1,N)+1;

% find coordinates of ghost points in (3 by 3) supercell
% (ghost point position is given by  
% X(ghost,:) = X(pim(ghost),:)+L*[pxind(ghost); pyind(ghost)];)

pind = floor((iXst-1)/N);
pxind = zeros(numel(iXst),1);
pyind = zeros(numel(iXst),1);

pxind(pind == 0 | pind == 4 | pind == 7) = 0;
pxind(pind == 1 | pind == 3 | pind == 6) = -1;
pxind(pind == 2 | pind == 5 | pind == 8) = 1;

pyind(pind == 0 | pind == 1 | pind == 2) = 0;
pyind(pind == 3 | pind == 4 | pind == 5) = -1;
pyind(pind == 6 | pind == 7 | pind == 8) = 1;

TVplt = TV;
%%% DELETE this potentially
%Xplt = Expression2;
%reindex triangles
TVplt(TV>N) = N+(1:szet)';
X = Xplt;
TV = TVplt;

X = R*X;
Xplt = X;
[szXX scratch]=size(Xplt);

tr = triangulation(TV,X(:,1),X(:,2));
e = tr.edges;
eplt = e;

% plot triangulation
figure(1); clf;
triplot(tr);
axis equal

% plot edges inside of unit square
title('Delauney Triangulation')
Xplt = [Xplt; 0 0; 0 R; R 0; R R];
eplt = [eplt; szXX+1 szXX+2; szXX+1 szXX+3; szXX+2 szXX+4; szXX+3 szXX+4];
figure(2); clf;
patch('faces', eplt, 'vertices', Xplt, 'edgecolor', 'b','Marker','o');
axis equal;
% plot ghost points
hold on
scatter(X(N+1:end,1),X(N+1:end,2),50,'r','filled');

%%%%%%%%%% plot on torus to see periodicity
Xplt=2*pi/R*Xplt;
Xplt=[cos(Xplt(:,1)).*(2+cos(Xplt(:,2))) sin(Xplt(:,1)).*(2+cos(Xplt(:,2))) sin(Xplt(:,2))];
figure(3); clf;
trisurf(TV, Xplt(:,1), Xplt(:,2), Xplt(:,3), 'FaceColor', 'cyan','FaceAlpha', 0.8);

hold on
Xplt_dual=2*pi/R*X_dual;
Xplt_dual=[cos(Xplt_dual(:,1)).*(2+cos(Xplt_dual(:,2)))...
    sin(Xplt_dual(:,1)).*(2+cos(Xplt_dual(:,2))) sin(Xplt_dual(:,2))];
patch('faces', e_dual, 'vertices', Xplt_dual, 'edgecolor', 'b','Marker','o');

figure(4);
axis equal;
% plot ghost points
hold on
N = numb_points;
title('voronoi')
patch('faces', e_dual, 'vertices', X_dual, 'edgecolor', 'b','Marker','o');

figure(5);
axis equal;
hold on
patch('faces', Edges, 'vertices', Vertices, 'edgecolor', [0 0 0],'Marker','.', 'MarkerSize', 18, 'MarkerEdgeColor', [.5 .5 .5]);







%% Set up data structures and find parabolic curcves

%swap original and dual vertices
%uncomment this to revert to original
%N = N_dual;
%numb_points = N;
pim = pim_dual;
%e = e_dual;
%X = X_dual;


Edges( isnan(Edges) )=-1; % fix edges
N = max(size(Vertices)); % number of points
numb_points = N;
e = Edges;
X = Vertices;
%%%% build a table CV with all vertex connections

%%% find max number of connections
N_vertex = N;
N_ghost = max(size(pim));
N = N_ghost + N_vertex;

num_cons = zeros(N,1);% N = number of vertices
ne = max(size(e));

for i=1:ne
    num_cons(e(i,1)) = num_cons(e(i,1)) + 1;
    num_cons(e(i,2)) = num_cons(e(i,2)) + 1;
end

max_cons = max(num_cons);

%%% build CV

CV = -1 * ones(N, max_cons);
num_cons = zeros(N,1);% N = number of vertices

for i = 1:ne

   num_cons(e(i,1)) = num_cons(e(i,1)) + 1;
   num_cons(e(i,2)) = num_cons(e(i,2)) + 1;  
   
   CV(e(i,1), num_cons(e(i,1))) = e(i,2);
   CV(e(i,2), num_cons(e(i,2))) = e(i,1);

end

%%% build EV - used for keeping track of edge indices

EV = zeros(ne,2);
[n_total, max_connections] = size(CV);
k=1;
for i = 1:n_total
    for j = 1:max_connections
        edge = CV(i,j);

        if edge > -1

             EV(k,1) = i;
             EV(k,2) = edge;
             k = k+1;
        end
    end
end

   
%% Random Walk simulation - basic approach
starting_index = 30;
delta_t = .1;
vertex_index = starting_index;
N_total = N;
vertex_frequency = zeros(1,N_total);
iters = 1E7
for i = 1:iters
    vertex_next_indices = CV(vertex_index,:);
    vertex_next_indices=vertex_next_indices(vertex_next_indices~=-1);
    % flip coin to select next vertex
    r = rand();
    coin = randi([1, length(vertex_next_indices)]);

    vertex_next_index = vertex_next_indices(coin);
  
    x0 = X(vertex_index,:);
    x1 = X(vertex_next_index,:); 

    p_return = probability_return(x0,x1,delta_t);
    p_transition(i) = p_return;
    r2 = rand();
    
    if r2 < p_return % vertex should return
       %do nothing   
       
    else
        % vertex should move on
        vertex_index = vertex_next_index;
    end
    vertex_frequency(vertex_next_index) = vertex_frequency(vertex_next_index)+1;
end

vfa = (vertex_frequency)/max(vertex_frequency);

%%


fvc2 = zeros(max(size(Vertices)),3); % blue
for i=1:length(fvc2)
   freq = 1-vfa(i);   
   fvc2(i,:) = 1.*[1-freq 0 freq];
   if freq == 1
       fvc2(i,:) = .9*[1 1 1];
   end
end

figure(18)
clf;
patch('vertices', X, 'FaceVertexCData', fvc2, ...
      'faces', EV, ...
      'EdgeColor', 'interp', ...
      'Marker','o', 'MarkerFaceColor', 'flat','LineWidth',4, 'MarkerSize', 4);
  
  figure(7)
clf;
patch('vertices', X, 'FaceVertexCData', fvc2, ...
      'faces', EV, ...
      'EdgeColor', 'interp', ...
      'Marker','o', 'MarkerFaceColor', 'flat','LineWidth',4, 'MarkerSize', 6);
  xlim([.35 .65])
  ylim([.35 .65])
  %%
for i = 1:max(size(CV))
    n_im = sum(CV(i,:) == -1);
    if n_im == 4
    n_connections(i) = 0;
    end
    if n_im == 3
    n_connections(i) = 1;
    end
    if n_im == 2
    n_connections(i) = 2;
    end
    if n_im == 1
    n_connections(i) = 3;
    end
    if n_im == 0
    n_connections(i) = 4;
    end
end
vfam = [];
ncm = [];
for i = 1:length(vfa)
    if vfa(i)~=0
        
        vfam=[vfam vfa(i)];
        ncm = [ncm n_connections(i)];
    end

end
ncm_char = [];
for i = 1:length(ncm)
    if ncm(i) == 1
    ncm_char =  [ncm_char; '1'];
    end
    if ncm(i) == 2
    ncm_char =  [ncm_char; '2'];
    end
    if ncm(i) == 3
    ncm_char =  [ncm_char; '3'];
    end
    if ncm(i) == 4
    ncm_char =  [ncm_char; '4'];
    end
end

t = ncm_char(1);
ncm_char(1) = ncm_char(4)
ncm_char(4) = t;

t = vfam(1);
vfam(1) = vfam(4)
vfam(4) = t;


  close all;
  figure(19)
boxplot(vfam, ncm_char)
xlabel('number of connections')
ylabel('relative frequency')

%% 
starting_index = 30;
delta_t = .1;
vertex_index = starting_index;
N_total = N;
vertex_frequency = zeros(1,N_total);
iters = 1E4
for i = 1:iters
    vertex_next_indices = CV(vertex_index,:);
    
    if sum(vertex_next_indices==-1) > 0
        % hit an edge
        vertex_next_indices
    end
    vertex_next_indices=vertex_next_indices(vertex_next_indices~=-1);
    
    % flip coin to select next vertex
    r = rand();
    coin = randi([1, length(vertex_next_indices)]);

    vertex_next_index = vertex_next_indices(coin);
  
    x0 = X(vertex_index,:);
    x1 = X(vertex_next_index,:); 

    p_return = probability_return(x0,x1,delta_t);
    p_transition(i) = p_return;
    r2 = rand();
    
    if r2 < p_return % vertex should return
       %do nothing   
       
    else
        % vertex should move on
        vertex_index = vertex_next_index;
    end
    vertex_frequency(vertex_next_index) = vertex_frequency(vertex_next_index)+1;
end

vfa = (vertex_frequency)/max(vertex_frequency);
