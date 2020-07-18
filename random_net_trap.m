% generate random network on periodic plane
% written by Thomas Fai (tfai@seas.harvard.edu)
% and Steven Petteruti, September 2016

clear all
format long
close all

numb_points = 100;
N = numb_points;
R = 1;

%%%%%%%%%%%%%%%%%%%%%% using random points in plane

%% regular triangulation
% ruler=linspace(0,1-1/sqrt(N),sqrt(N));
% 
% V=repmat(ruler,1,sqrt(N));
% U=reshape(V,sqrt(N),sqrt(N)); x=reshape(U',1,N);
% y=V;

%% random graph
  x = rand(1,N);
  y = rand(1,N);
  
 dt=DelaunayTri(x',y');
 X=[x' y'];
 
 %% Make graph on periodic plane
 % Place periodic copies
 X=[x' y'; x'-1 y'; x'+1 y';
    x'-1 y'-1; x' y'-1; x'+1 y'-1;
    x'-1 y'+1; x' y'+1; x'+1 y'+1];

%perform Delaunay triangulation
tr=DelaunayTri(X);
[szT scratch]=size(tr);
TV=tr(1:szT,:);

%% Build dual graph

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

%% Continue with original graph

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

%% Set up data structures and find parabolic curcves

%swap original and dual vertices
%uncomment this to revert to original
N = N_dual;
numb_points = N;
pim = pim_dual;
e = e_dual;
X = X_dual;

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

%%% compute S = [a0 a1 a2; a0' a1' a2'];
 EV_sort = sort(EV,2);
[EV_condensed, IA, IC] = unique(EV_sort, 'rows');


[n_total, m] = size(EV_condensed);
P = zeros(1,6);
for i=1:n_total
   x_ind_1 = EV_condensed(i,1);
   x_ind_2 = EV_condensed(i,2);
   
   x1 = X(x_ind_1,:)';
   x2 = X(x_ind_2,:)';

  
   [a0 a1 a2] = quadratic_curve(x1, x2, pi/2, 0);
   
   V = [a0 a1 a2];
   v = reshape(V,[1,6]);
   P = [P;v];
   

end

%% plotting the curves
[n_total, m] = size(EV_condensed);

for i=1:n_total

    a0 = [P(i,1); P(i,2)];
    a1 = [P(i,3); P(i,4)];
    a2 = [P(i,5); P(i,6)];

    s_array = linspace(0,1,500);
    curve = [];
    for i=1:length(s_array)
        s = s_array(i);
        c = [a0(1) + a1(1)*s + a2(1)*s^2 ; a0(2) + a1(2)*s + a2(2)*s^2];
        curve = [curve c];
    end
      
    figure(4)
    hold on
    plot(curve(1,:), curve(2,:),'b');
    
end

figure(4);
axis equal;
% plot ghost points
hold on
N = numb_points;
scatter(X(1:N,1), X(1:N,2),50,'b', 'filled');
scatter(X(N+1:end,1),X(N+1:end,2),50,'r','filled');

%%
xp = X(1,:)';
x1 = X(2,:)';
x2 = X(3,:)';
x3 = X(4,:)';
   
theta1 = 0; 
theta2 = 0;
theta3 = 0;
   
theta = theta_optimal(xp, x1, x2, x3, theta1, theta2, theta3);
   
   


%% Random Walk simulation: advanced

starting_index = 30;
delta_t = .001;
v_ind_current = starting_index;
edge_num = zeros(1,max(size(EV)));

for i = 1:200000
    vertex_next_indices = CV(v_ind_current,:);
    % flip coin to select next vertex
    r_random = rand();
    
    x0 = X(v_ind_current,:);
    x1 = X(vertex_next_indices(1),:);
    x2 = X(vertex_next_indices(2),:);
    x3 = X(vertex_next_indices(3),:);
    
    [p1, p2, p3] = weighted_coin(x0,x1,x2,x3);
    
    if r_random < p1
        coin = 1;
    elseif r_random < (p1+p2)
        coin = 2;   
    else
        coin = 3;
    end

    v_ind_next = vertex_next_indices(coin);
    x_next = X(vertex_next_index,:);
    
    p_tot = 0;
    for j=1:3
        v_ind_j = vertex_next_indices(j);
        x_end = X(v_ind_j,:);
        p_tot = p_tot + probability_return(x0, x_end, delta_t);
    end
    r = 1/3 * p_tot;
    
    for j=1:3
        
        v_ind_j = vertex_next_indices(j);
        [~,indx]=ismember([v_ind_current v_ind_j],EV,'rows');
        
        x_end = X(v_ind_j,:);
        
        E_nr = expected_steps_given_ruin(x0, x_end, delta_t);
        E_r = expected_steps_given_win(x0, x_end, delta_t);
        p_return = probability_return(x0, x_end, delta_t);
        
        edge_num(indx) = edge_num(indx) + ((p_return*E_r)/(3*(1-r))); %fill this in later

        if v_ind_j==v_ind_next
             edge_num(indx) = edge_num(indx) + E_nr; %fill this in later
        end
    end
    
    v_ind_current = v_ind_next;
    
    if v_ind_current > N
       % we are at a ghost point, go to a non-ghost point
       ghost_index = mod(v_ind_current,N);
       real_index = pim(ghost_index);
       v_ind_current  = real_index;
    end
end


EV_sort = sort(EV,2);
[EV_condensed, IA, IC] = unique(EV_sort, 'rows');

edge_num_condensed = zeros(1,max(size(EV_condensed))); %row vector

for i = 1:length(edge_num)
    
   edge_num_condensed(IC(i)) = edge_num_condensed(IC(i)) + edge_num(i);
    
end

vertex_num = zeros(N+N_ghost,1);
for i = 1:ne
    vert_index = EV_condensed(i,1);
    if (vert_index > N)
        vert_index = pim(vert_index-N);
    end
    vertex_num(vert_index) = vertex_num(vert_index) + edge_num_condensed(i)*.5;
    
    vert_index = EV_condensed(i,2);
    if (vert_index > N)
        vert_index = pim(vert_index-N);
    end
    vertex_num(vert_index) = vertex_num(vert_index) + edge_num_condensed(i)*.5;
end

for i = 1:N_ghost
    vertex_num(N+i) = vertex_num(pim(i));
end

fvc = zeros(N+N_ghost,3);
vertex_num_freq = vertex_num/max(vertex_num);

for i=1:max(size(fvc))
   freq = 1-vertex_num_freq(i);   
   fvc(i,:) = 1.0*[freq freq freq];
end

vertex_frequency = vertex_frequency + 2.3*vertex_num';
vertex_frequency_adjusted = (vertex_frequency)/max(vertex_frequency);

fvc2 = zeros(max(size(X_dual)),3); % blue
for i=1:max(size(fvc2))
   freq = 1-vertex_frequency_adjusted(i);   
   fvc2(i,:) = 1.*[freq freq freq];
end

figure(17)
title('Edge Method')
subplot(2,2,1)
plot(vertex_num_freq,'o')
xlabel('vertex number (200+ = ghost)')
subplot(2,2,2)
title('advanced method')
patch('vertices', X_dual, 'FaceVertexCData', fvc, ...
      'faces', e_dual, ...
      'EdgeColor', 'interp', ...
      'Marker','o', 'MarkerFaceColor', 'flat','LineWidth',3);
  
%% Traps
starting_index = 30;
delta_t = .001;
v_ind_current = starting_index;
edge_num = zeros(1,max(size(EV)));
traps = [18 40];
bias = 1;
count = zeros(1,N+N_ghost);
for i = 1:100000
    vertex_next_indices = CV(v_ind_current,:);
    % flip coin to select next vertex
    r_random = rand();
    
    x0 = X(v_ind_current,:);
    x1 = X(vertex_next_indices(1),:);
    x2 = X(vertex_next_indices(2),:);
    x3 = X(vertex_next_indices(3),:);

    for k=1:length(traps)
        
        trap = traps(k);
        
    if sum(find(trap==vertex_next_indices))~=0 
       
        index = find(vertex_next_indices==trap);
       
        
        [p1, p2, p3] = weighted_coin_trap(x0,x1,x2,x3,bias,delta_t, index);
        
        M = [p1 p2 p3 index trap];

    else
         [p1, p2, p3] = weighted_coin(x0,x1,x2,x3);

    end
    end

    
    if r_random < p1
        coin = 1;
    elseif r_random < (p1+p2)
        coin = 2;   
    else
        coin = 3;
    end

    v_ind_next = vertex_next_indices(coin);
    x_next = X(v_ind_next,:);
    
    p_tot = 0;
    for j=1:3
        v_ind_j = vertex_next_indices(j);
        x_end = X(v_ind_j,:);
        p_tot = p_tot + probability_return(x0, x_end, delta_t);
    end
    r = 1/3 * p_tot;
    
    for j=1:3
        
        v_ind_j = vertex_next_indices(j);
        [~,indx]=ismember([v_ind_current v_ind_j],EV,'rows');
        
        x_end = X(v_ind_j,:);
        
        E_nr = expected_steps_given_ruin(x0, x_end, delta_t);
        E_r = expected_steps_given_win(x0, x_end, delta_t);
        p_return = probability_return(x0, x_end, delta_t);
        
        edge_num(indx) = edge_num(indx) + ((p_return*E_r)/(3*(1-r))); %fill this in later

        if v_ind_j==v_ind_next
             edge_num(indx) = edge_num(indx) + E_nr; %fill this in later
        end
    end
    
    v_ind_current = v_ind_next;
    
    if v_ind_current > N
       % we are at a ghost point, go to a non-ghost point
       ghost_index = mod(v_ind_current,N);
       real_index = pim(ghost_index);
       v_ind_current  = real_index;
    end
    
    count(v_ind_current) = count(v_ind_current) + 1;
end


EV_sort = sort(EV,2);
[EV_condensed, IA, IC] = unique(EV_sort, 'rows');

edge_num_condensed = zeros(1,max(size(EV_condensed))); %row vector

for i = 1:length(edge_num)
    
   edge_num_condensed(IC(i)) = edge_num_condensed(IC(i)) + edge_num(i);
    
end

vertex_num = zeros(N+N_ghost,1);
for i = 1:ne
    vert_index = EV_condensed(i,1);
    if (vert_index > N)
        vert_index = pim(vert_index-N);
    end
    vertex_num(vert_index) = vertex_num(vert_index) + edge_num_condensed(i)*.5;
    
    vert_index = EV_condensed(i,2);
    if (vert_index > N)
        vert_index = pim(vert_index-N);
    end
    vertex_num(vert_index) = vertex_num(vert_index) + edge_num_condensed(i)*.5;
end

for i = 1:N_ghost
    vertex_num(N+i) = vertex_num(pim(i));
end

fvc = zeros(N+N_ghost,3);
vertex_num_freq = vertex_num/max(vertex_num);
count = count/max(count);

for i=1:max(size(fvc))
   freq = 1-vertex_num_freq(i);   
   fvc(i,:) = 1.0*[freq freq freq];
   fvc(18,:) = [1 0 1];
   fvc(40,:) = [.8 0 .8];
end

for i=1:max(size(fvc))
   freq = 1-count(i);   
   fvc(i,:) = 1.0*[freq freq freq];
   
end


figure(13)

title('advanced method')
patch('vertices', X_dual, 'FaceVertexCData', fvc, ...
      'faces', e_dual, ...
      'EdgeColor', 'interp', ...
      'Marker','o', 'MarkerFaceColor', 'flat','LineWidth',3);
  
  %% Traps Simple
  
  starting_index = 30;
delta_t = .001;
v_ind_current = starting_index;
edge_num = zeros(1,max(size(EV)));
traps = [13];
bias = 1;
count = zeros(1,N+N_ghost);
for i = 1:1000000
    
            count(v_ind_current) = count(v_ind_current) + 1;

            
    vertex_next_indices = CV(v_ind_current,:);
    % flip coin to select next vertex
    r_random = rand();
    
    x0 = X(v_ind_current,:);
    x1 = X(vertex_next_indices(1),:);
    x2 = X(vertex_next_indices(2),:);
    x3 = X(vertex_next_indices(3),:);

    for k=1:length(traps)
        
        trap = traps(k);
        
    if sum(find(trap==vertex_next_indices))~=0 
       
        index = find(vertex_next_indices==trap);
       
        
        [p1, p2, p3] = weighted_coin_trap(x0,x1,x2,x3,bias,delta_t, index);
        
    else
         [p1, p2, p3] = weighted_coin(x0,x1,x2,x3);

    end
    end

    
    if r_random < p1
        coin = 1;
    elseif r_random < (p1+p2)
        coin = 2;   
    else
        coin = 3;
    end

    v_ind_next = vertex_next_indices(coin);
    x_next = X(v_ind_next,:);
    

  
    v_ind_current = v_ind_next;
    
    if v_ind_current > N
       % we are at a ghost point, go to a non-ghost point
       ghost_index = mod(v_ind_current,N);
       real_index = pim(ghost_index);
       v_ind_current  = real_index;
    end
    
end


fvc = zeros(N+N_ghost,3);
count = count/max(count);

%%
for i=1:max(size(fvc))
   freq = 1-count(i);   
   fvc(i,:) = 1.0*[freq freq freq];
   if freq > .1
          fvc(i,:) = .96*[freq freq freq];

   end
   
end

close all
figure(6)

plot(count(1:200),'o')
hax = axes
SP = 13;
%line([SP SP], [0 1], 'color', [1 0 0], 'LineStyle', ':')
ylim([0 1])
xlim([0 200])

figure(9)
patch('vertices', X_dual, 'FaceVertexCData', fvc, ...
      'faces', e_dual, ...
      'EdgeColor', 'interp', ...
      'Marker','o', 'MarkerFaceColor', 'flat','LineWidth',3);