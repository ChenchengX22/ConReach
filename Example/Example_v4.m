clc; clear all;

%% Attitude dynamics
% dotR = R hat{w}
% dotw = Aw
A = [-2,0,0;0,-1,0;0,0,-3];
% A = [-2,0,0;0,-4,0;0,0,-3];
% A = [-2,-0.2,-0.1;-0.2,-1,0;-0.1,0,-3];

% Reachability settings
T = 4; % terminal time T
% initial condition
R_0 = eye(3);
R_r = 0.1; Q0 = eye(3);
[~,~,expr_0] = get_exp_coord(R_0,0,false);
% w: center = w_c0 radius = 0.2 P = I
w_c0 = [0.65,0.54,0.61];
w_r = 0.1; P0 = eye(3);

yalmip('clear');
%% Simulation
sim_step = 0.001; sim_N = T/sim_step;

% Nominal trajectory: starts from the center of the initial set.
w_rec = zeros(sim_N+1,3); w_rec(1,:) = w_c0(:);
dw_rec = zeros(sim_N+1,3); 
R_rec = zeros(sim_N+1,3,3); R_rec(1,:,:) = R_0(:,:);
expr_rec = zeros(sim_N+1,3); expr_rec(1,:) = expr_0(:);
Cur_r = zeros(3,1); Cur_w = zeros(3,1);
dotw = zeros(3,1);
for i = 1:sim_N
    dotw(:) = A*w_rec(i,:)'; dw_rec(i,:) = dotw(:);
    w_rec(i+1,:) = w_rec(i,:) + dotw(:)'*sim_step;

    Cur_r(:) = expr_rec(i,:); Cur_w(:) = w_rec(i,:);
    [Next_r,Next_Ang_Vel] = Example_Dynamic(Cur_r,Cur_w,sim_step,A);
    expr_rec(i+1,:) = Next_r(:);
    R_rec(i+1,:,:) = exp_coord2rot(Next_r,0);
end
dotw = A*w_rec(end,:)'; dw_rec(end,:) = dotw(:);

% some nearby trajectories
N_traj = 4;

% tmpa = (1-2*rand(1,3)); a = tmpa/norm(tmpa); 
% a1 = R_r*a; R_1 = R_0*expfun(skewfun(a1)); 
tmpa = [0.004,0.071,0.071]; a = tmpa/norm(tmpa); 
a1 = R_r*a; R_1 = R_0*expfun(skewfun(a1));
[~,~,expr_1] = get_exp_coord(R_1,0,false);

% tmpa = (1-2*rand(1,3)); a = tmpa/norm(tmpa); 
% a2 = R_r*a; R_2 = R_0*expfun(skewfun(a2)); 
tmpa = [-0.066,0.048,-0.058]; a = tmpa/norm(tmpa); 
a2 = R_r*a; R_2 = R_0*expfun(skewfun(a2)); 
[~,~,expr_2] = get_exp_coord(R_2,0,false);

% tmpa = (1-2*rand(1,3)); a = tmpa/norm(tmpa); 
% a3 = R_r*a; R_3 = R_0*expfun(skewfun(a3));
tmpa = [0.071,0.043,0.055]; a = tmpa/norm(tmpa); 
a3 = R_r*a; R_3 = R_0*expfun(skewfun(a3));
[~,~,expr_3] = get_exp_coord(R_3,0,false);

% tmpa = (1-2*rand(1,3)); a = tmpa/norm(tmpa); 
% a4 = R_r*a; R_4 = R_0*expfun(skewfun(a4));
tmpa = [0.063,-0.042,0.065]; a = tmpa/norm(tmpa); 
a4 = R_r*a; R_4 = R_0*expfun(skewfun(a4));
[~,~,expr_4] = get_exp_coord(R_4,0,false);

% tmpa = (1-2*rand(1,3)); a = tmpa/norm(tmpa); 
% w_c1 = w_c0 + w_r*a;
% tmpa = (1-2*rand(1,3)); a = tmpa/norm(tmpa); 
% w_c2 = w_c0 + w_r*a;
% tmpa = (1-2*rand(1,3)); a = tmpa/norm(tmpa); 
% w_c3 = w_c0 + w_r*a;
% tmpa = (1-2*rand(1,3)); a = tmpa/norm(tmpa); 
% w_c4 = w_c0 + w_r*a;
w_c1 = [0.623,0.487,0.530];
w_c2 = [0.659,0.449,0.650]; 
w_c3 = [0.684,0.634,0.621];
w_c4 = [0.710,0.568,0.535];

w_rec_n = zeros(N_traj,sim_N+1,3); 
w_rec_n(1,1,:) = w_c1(:); 
w_rec_n(2,1,:) = w_c2(:);
w_rec_n(3,1,:) = w_c3(:); 
w_rec_n(4,1,:) = w_c4(:);
dw_rec_n = zeros(N_traj,sim_N+1,3);

R_rec_n = zeros(N_traj,sim_N+1,3,3); 
R_rec_n(1,1,:,:) = R_1(:,:); 
R_rec_n(2,1,:,:) = R_2(:,:);
R_rec_n(3,1,:,:) = R_3(:,:); 
R_rec_n(4,1,:,:) = R_4(:,:);

expr_rec_n = zeros(N_traj,sim_N+1,3); 
expr_rec_n(1,1,:) = expr_1(:); 
expr_rec_n(2,1,:) = expr_2(:);
expr_rec_n(3,1,:) = expr_3(:); 
expr_rec_n(4,1,:) = expr_4(:);

tmpw = zeros(3,1);
for j = 1:N_traj
    for i = 1:sim_N
        tmpw(:) = w_rec_n(j,i,:);
        dotw = A*tmpw;
        dw_rec_n(j,i,:) = dotw(:);
        w_rec_n(j,i+1,:) = w_rec_n(j,i,:) + reshape(dotw,1,1,3)*sim_step;

        Cur_r(:,1) = expr_rec_n(j,i,:);
        Cur_w(:,1) = w_rec_n(j,i,:);
        [Next_r,Next_Ang_Vel] = Example_Dynamic(Cur_r,Cur_w,sim_step,A);
        expr_rec_n(j,i+1,:) = Next_r(:);
        R_rec_n(j,i+1,:,:) = exp_coord2rot(Next_r,0);
    end
    tmpw(:) = w_rec_n(j,end,:);
    dotw = A*tmpw; dw_rec_n(j,end,:) = dotw(:);
end

%% Reachtube computation: for angular velocity only
% First_iteration: reachtube of w 
% get standard contraction rate, i.e., 
% metric = I, | | is computed by 2-norm 
% |w_0-w_c0|< w_r  -> |w(t)-w_c(t)|< w_r*e^{k*(t-t0)} 
% 
% Reference: W. Lohmiller and J.-J. E. Slotine, “On Contraction
% Analysis for Non-linear Systems,” Automatica, vol. 34, no. 6, pp.
% 683–696, Jun. 1998, doi: 10.1016/S0005-1098(98)00019-3.

tic
sdpvar gamma;
Const1 = [A'+A <= 2*gamma*eye(3)];
Objective1 = gamma;
ops = sdpsettings('solver','sdpt3','verbose',0);
optimize(Const1,Objective1,ops);
w_gamma_max = value(gamma);

% Compute and plot ... at every interval
split_t_int = 0.1; 
split_step = split_t_int/sim_step;
split_Num = sim_N/split_step;

% Compute the reachtube in the form of boxes
x = zeros(2,3); center = zeros(3,1);
w_lb = 10*ones(split_Num,3); 
w_ub = -10*ones(split_Num,3);
for j = 1:split_Num
    for l = 1:split_step+1
        radius = w_r*exp(w_gamma_max*((j-1)*split_step+(l-1))*sim_step);
        center(:) = w_rec((j-1)*split_step+l,:);
        for i = 1:3
            x(1,i) = center(i) - radius;
            x(2,i) = center(i) + radius;
            if x(1,i) < w_lb(j,i)
                w_lb(j,i) = x(1,i);
            end
            if x(2,i) > w_ub(j,i)
                w_ub(j,i) = x(2,i);
            end
        end
    end
end
toc

%% ConReach: Reach set computation for R and w 
% Note: Computing a global contraction rate is not ideal. It's better to
% approximate the reachable set over small time intervals and iterate.
% (Here, we start with the identity metric and perform two iterations.)

% contraction rate
c = zeros(split_Num,1);
sdpvar c1;

% optimization setting
opt_c = zeros(split_Num,1); 
opt_R_Q = zeros(split_Num,3,3);
opt_w_P = zeros(split_Num,3,3);
c_opt_step = 0.01; 
N_linesearch = 3;

Q_R = sdpvar(3,3); 
P_w = sdpvar(3,3);

R_Q_init = eye(3);
w_P_init = eye(3);

tic
for i = 1:split_Num
    Objective2 = c1;
    Const2 = const_W(w_lb(i,:),w_ub(i,:),eye(3),eye(3),A,c1);
    sol1_I = optimize(Const2,Objective2,ops);
    c(i) = value(c1);
    tmp_c = c(i);

    % search for a smaller contraction rate
    opt_c(i) = tmp_c;
    opt_R_Q(i,:,:) = eye(3);
    opt_w_P(i,:,:) = eye(3);
    for j = 1:N_linesearch
        tmp_c = tmp_c - c_opt_step;

        % Constraint: cover the initial set
        Const_opt = [Q_R>=0,P_w>=0,Q_R<=R_Q_init,P_w<=w_P_init];
        Const_PQ = const_W(w_lb(i,:),w_ub(i,:),Q_R,P_w,A,tmp_c);
        Const_PQ = [Const_opt,Const_PQ];
        % obj3: minimize the initial ball (on SO(3))
        Objective3 = -trace(Q_R);
        sol2_Q = optimize(Const_PQ,Objective3,ops);
        tmp_Q = value(Q_R); tmp_P = value(P_w);
        if sol2_Q.problem == 0 && Check_Matrix_Definiteness(tmp_Q) ...
                               && abs(tmp_Q(1,1))> 0.01 && abs(tmp_Q(2,2))> 0.01 && abs(tmp_Q(3,3))> 0.01...
                               && Check_Matrix_Definiteness(tmp_P) ...
                               && abs(tmp_P(1,1))> 0.01 && abs(tmp_P(2,2))> 0.01 && abs(tmp_P(3,3))> 0.01
            opt_c(i) = tmp_c;
            opt_R_Q(i,:,:) = tmp_Q(:,:);
            opt_w_P(i,:,:) = tmp_P(:,:);
        end
    end
    R_Q_init(:,:) = opt_R_Q(i,:,:);
    w_P_init(:,:) = opt_w_P(i,:,:);
end
toc



%% visualization
figure1 = figure('Name','Angular velocity trajectory and reachtube: [1st iteration]');

% plot angular velocity trajectory
h1 = plot3(w_rec(:,1),w_rec(:,2),w_rec(:,3),'k','LineWidth',1.2);
hold on; axis on; grid on;
h2 = plot3(w_rec(1,1),w_rec(1,2),w_rec(1,3),'b*','MarkerSize',8);
h3 = plot3(w_rec(end,1),w_rec(end,2),w_rec(end,3),'rpentagram','MarkerSize',8);
for j = 1:N_traj
    h4(j) = plot3(w_rec_n(j,:,1),w_rec_n(j,:,2),w_rec_n(j,:,3));
    plot3(w_rec_n(j,1,1),w_rec_n(j,1,2),w_rec_n(j,1,3),'b*','MarkerSize',8);
    plot3(w_rec_n(j,end,1),w_rec_n(j,end,2),w_rec_n(j,end,3),'rpentagram','MarkerSize',8);
end

% plot the reachtube
plot_step = 20; plot_Num = sim_N/plot_step;
radius_rec = zeros(plot_Num+1,1);
for i = 1:plot_Num+1
    radius = w_r*exp(w_gamma_max*((i-1)*plot_step*sim_step)); 
    radius_rec(i) = radius;
    center(:) = w_rec((i-1)*plot_step+1,:);
    [x, y, z] = ellipsoid(center(1),center(2),center(3),radius,radius,radius,30);
    surf(x, y, z,'LineStyle','none','FaceColor',"b",'Facealpha',0.03)
end

% plot the reach set
% selected time interval
show_rec_Num = [1,7];

faces = [1 2 3 4;5 6 7 8;4 3 6 5;3 2 7 6;2 1 8 7;1 4 5 8];
verts = zeros(3,8);
for j = show_rec_Num
    % verts(:,1) = [w_ub(j,1),w_lb(j,2),w_ub(j,3)];
    % verts(:,2) = [w_ub(j,1),w_ub(j,2),w_ub(j,3)];
    % verts(:,3) = [w_lb(j,1),w_ub(j,2),w_ub(j,3)];
    % verts(:,4) = [w_lb(j,1),w_lb(j,2),w_ub(j,3)];
    % verts(:,5) = [w_lb(j,1),w_lb(j,2),w_lb(j,3)];
    % verts(:,6) = [w_lb(j,1),w_ub(j,2),w_lb(j,3)];
    % verts(:,7) = [w_ub(j,1),w_ub(j,2),w_lb(j,3)];
    % verts(:,8) = [w_ub(j,1),w_lb(j,2),w_lb(j,3)];
    % patch('Faces',faces,'Vertices',verts','FaceColor','blue','FaceAlpha',0.02)
    plot3([w_lb(j,1),w_ub(j,1)],[w_lb(j,2),w_lb(j,2)],[w_lb(j,3),w_lb(j,3)],'b--');
    plot3([w_lb(j,1),w_ub(j,1)],[w_lb(j,2),w_lb(j,2)],[w_ub(j,3),w_ub(j,3)],'b--');
    plot3([w_lb(j,1),w_ub(j,1)],[w_ub(j,2),w_ub(j,2)],[w_lb(j,3),w_lb(j,3)],'b--');
    plot3([w_lb(j,1),w_ub(j,1)],[w_ub(j,2),w_ub(j,2)],[w_ub(j,3),w_ub(j,3)],'b--');
    plot3([w_lb(j,1),w_lb(j,1)],[w_lb(j,2),w_ub(j,2)],[w_lb(j,3),w_lb(j,3)],'b--');
    plot3([w_ub(j,1),w_ub(j,1)],[w_lb(j,2),w_ub(j,2)],[w_lb(j,3),w_lb(j,3)],'b--');
    plot3([w_lb(j,1),w_lb(j,1)],[w_lb(j,2),w_ub(j,2)],[w_ub(j,3),w_ub(j,3)],'b--');
    plot3([w_ub(j,1),w_ub(j,1)],[w_lb(j,2),w_ub(j,2)],[w_ub(j,3),w_ub(j,3)],'b--');
    plot3([w_lb(j,1),w_lb(j,1)],[w_lb(j,2),w_lb(j,2)],[w_lb(j,3),w_ub(j,3)],'b--');
    plot3([w_ub(j,1),w_ub(j,1)],[w_lb(j,2),w_lb(j,2)],[w_lb(j,3),w_ub(j,3)],'b--');
    plot3([w_lb(j,1),w_lb(j,1)],[w_ub(j,2),w_ub(j,2)],[w_lb(j,3),w_ub(j,3)],'b--');
    plot3([w_ub(j,1),w_ub(j,1)],[w_ub(j,2),w_ub(j,2)],[w_lb(j,3),w_ub(j,3)],'b--');
    
end

zlabel({'\omega_1'}); ylabel({'\omega_2'}); xlabel({'\omega_3'});
% axis equal;
% legend1 = legend([h1 h2 h3 h4(1) h4(2) h4(3) h4(4)],{'main trajectory','initial states','terminal states','1','2','3','4'});
legend1 = legend([h2 h3],{'sampled initial states','terminal states'});

set(legend1,'Location','northeast');
view([-81,11.4]);




figure3 = figure('Name','Rotation trajectory and reach set in exp charts: [1st iteration]');
h1 = plot3(expr_rec(:,1),expr_rec(:,2),expr_rec(:,3),'k','LineWidth',1.2);
hold on; axis on; grid on;
h2 = plot3(expr_rec(1,1),expr_rec(1,2),expr_rec(1,3),'b*','MarkerSize',8);
h3 = plot3(expr_rec(end,1),expr_rec(end,2),expr_rec(end,3),'rpentagram','MarkerSize',8);
for j = 1:N_traj
    h4(j) = plot3(expr_rec_n(j,:,1),expr_rec_n(j,:,2),expr_rec_n(j,:,3));
    plot3(expr_rec_n(j,1,1),expr_rec_n(j,1,2),expr_rec_n(j,1,3),'b*','MarkerSize',8);
    plot3(expr_rec_n(j,end,1),expr_rec_n(j,end,2),expr_rec_n(j,end,3),'rpentagram','MarkerSize',8);
end

% plot reach set for the 1st iteration [identity metric]
exp_plot_Num = split_step/plot_step;
% plot_List = [1,4];
% split_List = [1,25];

% for test
% exp_radius_rec = zeros(plot_Num+1,1);
% exp_radius_rec_n = zeros(N_traj,plot_Num+1,1);

Sample_Num = 200;
Sub_sample_Num = 10;
terminal_t = 1; 
geo_sim_N = terminal_t/sim_step;

x_traj = zeros(Sample_Num,6,geo_sim_N + 1);
x_R = zeros(3,1); x_R0 = zeros(3,1);
radius_tmp = 1;

num_mb = 80;
x = zeros(num_mb,num_mb);
y = zeros(num_mb,num_mb);
z = zeros(num_mb,num_mb);

radius_0 = sqrt(R_r^2 + w_r^2);
tmp_Q(:,:) = eye(3); tmp_P(:,:) = eye(3);
% N_center = zeros(3,1);
R_center = zeros(3,3);
NR_center = zeros(3,3);
last_plot_center = zeros(3,1);
last_plot_center(:) = expr_rec(1,:);
exp_Q_R = sdpvar(3,3);
dis_plot = 0.04;
for i = 1:split_Num
% for i = split_List
    for j = 1:exp_plot_Num
    % for j = plot_List
        radius = radius_0*exp(c(i)*(j-1)*plot_step*sim_step);
        % exp_radius_rec((i-1)*exp_plot_Num+j) = radius;
        center(:) = expr_rec((i-1)*split_step+(j-1)*plot_step+1,:);

        % Choose the time point for reach set visualization.
        dist = norm(last_plot_center - center);
        if dist < dis_plot
            continue;
        else
            last_plot_center(:) = center(:);
        end
        R_center(:,:) = R_rec((i-1)*split_step+(j-1)*plot_step+1,:,:);
        for k = 1:N_traj
            % N_center(:) = expr_rec_n(k,(i-1)*split_step+(j-1)*plot_step+1,:);
            NR_center(:,:) = R_rec_n(k,(i-1)*split_step+(j-1)*plot_step+1,:,:);
            % exp_radius_rec_n(k,(i-1)*exp_plot_Num+j) = norm(veemap(NR_center'*R_center));
        end

        % Plot the reach set 
        % We densely sampled the surface of the reach ball and approximated it by a Q-ball in R3
        const_rch = [exp_Q_R>=0];
        x_R0(:) = center(:);
        for l = 1:Sample_Num
            w_tmp = (1-2*rand(1,3)); w_tmp = w_tmp/norm(w_tmp);
            w_tmp = radius*w_tmp;
            x_traj(l,:,1) = [center(1),center(2),center(3),w_tmp(1),w_tmp(2),w_tmp(3)];
            for n = 1:geo_sim_N
                Cur_r(:) = x_traj(l,1:3,n);
                Cur_w(:) = x_traj(l,4:6,n);
                [Next_r,Next_Ang_Vel] = Geodesic_Dynamic_Exp(Cur_r,Cur_w,tmp_Q,sim_step);
                x_traj(l,1:3,n+1) = Next_r(:);
                x_traj(l,4:6,n+1) = Next_Ang_Vel(:);
            end
            x_R(:) = x_traj(l,1:3,end);
            const_rch = [const_rch,(x_R-x_R0)'*exp_Q_R*(x_R-x_R0)<=radius_tmp*radius_tmp];
        end
        Objective_exp_Q = -trace(exp_Q_R);
        sol_Q_exp = optimize(const_rch,Objective_exp_Q,ops);
        opt_Exp_Q = value(exp_Q_R);

        if sol_Q_exp.problem == 0 && Check_Matrix_Definiteness(opt_Exp_Q) && abs(opt_Exp_Q(1,1))> 0.01 && abs(opt_Exp_Q(2,2))> 0.01 && abs(opt_Exp_Q(3,3))> 0.01
            [V,D] = eig(opt_Exp_Q);
            A_plot = sqrt(D)*V';
            for l = 1:num_mb
                init_x = radius_tmp - (l-1)*2*radius_tmp/(num_mb-1);
                r = real(sqrt(radius_tmp^2-init_x^2));
                for k = 1:num_mb
                    theta = 2*pi*(k-1)/(num_mb-1);
                    init_z = r*sin(theta);
                    init_y = r*cos(theta);
                    init_r = [init_x;...
                              init_y;...
                              init_z];
                    r1 = (A_plot\init_r) + center;
                    % if ~isreal(r1(1)) || ~isreal(r1(2)) || ~isreal(r1(3))
                    %     shdoa = 1;
                    % end
                    x(l,k) = r1(1);
                    y(l,k) = r1(2);
                    z(l,k) = r1(3);
                end
            end
            surf(x, y, z,'LineStyle','none','FaceColor',"b",'Facealpha',0.03)
        end
    end
    radius_0 = radius_0*exp(c(i)*split_step*sim_step);
end

% Plot the reach set at terminal time
radius = radius_0;
% exp_radius_rec(end) = radius;
center(:) = expr_rec(end,:);
x_R0(:) = center(:);
const_rch = [exp_Q_R>=0];
for l = 1:Sample_Num
    w_tmp = (1-2*rand(1,3)); w_tmp = w_tmp/norm(w_tmp);
    w_tmp = radius*w_tmp;
    x_traj(l,:,1) = [center(1),center(2),center(3),w_tmp(1),w_tmp(2),w_tmp(3)];
    for n = 1:geo_sim_N
        Cur_r(:) = x_traj(l,1:3,n);
        Cur_w(:) = x_traj(l,4:6,n);
        [Next_r,Next_Ang_Vel] = Geodesic_Dynamic_Exp(Cur_r,Cur_w,tmp_Q,sim_step);
        x_traj(l,1:3,n+1) = Next_r(:);
        x_traj(l,4:6,n+1) = Next_Ang_Vel(:);
    end
    x_R(:) = x_traj(l,1:3,end);
    const_rch = [const_rch,(x_R-x_R0)'*exp_Q_R*(x_R-x_R0)<=radius_tmp*radius_tmp];
end
Objective_exp_Q = -trace(exp_Q_R);
sol_Q_exp = optimize(const_rch,Objective_exp_Q,ops);
opt_Exp_Q = value(exp_Q_R);
if sol_Q_exp.problem == 0 && Check_Matrix_Definiteness(opt_Exp_Q) && abs(opt_Exp_Q(1,1))> 0.01 && abs(opt_Exp_Q(2,2))> 0.01 && abs(opt_Exp_Q(3,3))> 0.01
    [V,D] = eig(opt_Exp_Q);
    A_plot = sqrt(D)*V';
    for l = 1:num_mb
        init_x = radius_tmp - (l-1)*2*radius_tmp/(num_mb-1);
        r = real(sqrt(radius_tmp^2-init_x^2));
        for k = 1:num_mb
            theta = 2*pi*(k-1)/(num_mb-1);
            init_z = r*sin(theta);
            init_y = r*cos(theta);
            init_r = [init_x;...
                      init_y;...
                      init_z];
            r1 = (A_plot\init_r) + center;
            % if ~isreal(r1(1)) || ~isreal(r1(2)) || ~isreal(r1(3))
            %     shdoa = 1;
            % end
            x(l,k) = r1(1);
            y(l,k) = r1(2);
            z(l,k) = r1(3);
        end
    end
    surf(x, y, z,'LineStyle','none','FaceColor',"b",'Facealpha',0.03)
end

zlabel({'r_{03}'}); ylabel({'r_{02}'}); xlabel({'r_{01}'});
legend3 = legend([h2 h3],{'sampled initial states','terminal states'});
set(legend3,'Location','northeast');
view([-81,11.4]);




figure4 = figure('Name','Rotation trajectory and reach set in exp charts: [final iteration]');
h1 = plot3(expr_rec(:,1),expr_rec(:,2),expr_rec(:,3),'k','LineWidth',1.2);
hold on; axis on; grid on;
h2 = plot3(expr_rec(1,1),expr_rec(1,2),expr_rec(1,3),'b*','MarkerSize',8);
h3 = plot3(expr_rec(end,1),expr_rec(end,2),expr_rec(end,3),'rpentagram','MarkerSize',8);
for j = 1:N_traj
    h4(j) = plot3(expr_rec_n(j,:,1),expr_rec_n(j,:,2),expr_rec_n(j,:,3));
    plot3(expr_rec_n(j,1,1),expr_rec_n(j,1,2),expr_rec_n(j,1,3),'b*','MarkerSize',8);
    plot3(expr_rec_n(j,end,1),expr_rec_n(j,end,2),expr_rec_n(j,end,3),'rpentagram','MarkerSize',8);
end

plot_Num = split_step/plot_step;
plot_List = [1,4];

% for test
% exp_radius_rec_1 = zeros(plot_Num+1,1);

radius_0 = sqrt(R_r^2 + w_r^2);
for i = 1:split_Num
    tmp_Q(:,:) = opt_R_Q(i,:,:);
    tmp_P(:,:) = opt_w_P(i,:,:);
    % for j = 1:plot_Num
    for j = plot_List
        radius = radius_0*exp(opt_c(i)*(j-1)*plot_step*sim_step);
        exp_radius_rec_1((i-1)*plot_Num+j) = radius;
        center(:) = expr_rec((i-1)*split_step+(j-1)*plot_step+1,:);

        % Choose the time point for reach set visualization.
        dist = norm(last_plot_center - center);
        if dist < dis_plot
            continue;
        else
            last_plot_center(:) = center(:);
        end 
        const = [];
        x_R0(:) = center(:);
        for l = 1:Sub_sample_Num
            for m = 1:Sub_sample_Num
                cur_num = Sub_sample_Num*(l-1)+m;
                w_tmp = (1-2*rand(1,3)); w_tmp = w_tmp/norm(w_tmp);
                w_tmp = radius*w_tmp;
                x_traj(cur_num,:,1) = [center(1),center(2),center(3),w_tmp(1),w_tmp(2),w_tmp(3)];
                for n = 1:geo_sim_N
                    Cur_r(:) = x_traj(cur_num,1:3,n);
                    Cur_w(:) = x_traj(cur_num,4:6,n);
                    [Next_r,Next_Ang_Vel] = Geodesic_Dynamic_Exp(Cur_r,Cur_w,tmp_Q,sim_step);
                    x_traj(cur_num,1:3,n+1) = Next_r(:);
                    x_traj(cur_num,4:6,n+1) = Next_Ang_Vel(:);
                end
                % plot3(squeeze(x_traj(cur_num,1,:)),squeeze(x_traj(cur_num,2,:)),squeeze(x_traj(cur_num,3,:)),"b");
                x_R(:) = x_traj(cur_num,1:3,end);
                const = [const,(x_R-x_R0)'*exp_Q_R*(x_R-x_R0)<=radius_tmp*radius_tmp];
            end
        end
        Objective_exp_Q = -trace(exp_Q_R);
        sol_Q_exp = optimize(const,Objective_exp_Q,ops);
        opt_Exp_Q = value(exp_Q_R);

        if sol_Q_exp.problem == 0 && Check_Matrix_Definiteness(opt_Exp_Q) && abs(opt_Exp_Q(1,1))> 0.01 && abs(opt_Exp_Q(2,2))> 0.01 && abs(opt_Exp_Q(3,3))> 0.01
            [V,D] = eig(opt_Exp_Q);
            A_plot = sqrt(D)*V';
            for l = 1:num_mb
                init_x = radius_tmp - (l-1)*2*radius_tmp/(num_mb-1);
                r = real(sqrt(radius_tmp^2-init_x^2));
                for k = 1:num_mb
                    theta = 2*pi*(k-1)/(num_mb-1);
                    init_z = r*sin(theta);
                    init_y = r*cos(theta);
                    init_r = [init_x;...
                        init_y;...
                        init_z];
                    r1 = (A_plot\init_r) + center;
                    if ~isreal(r1(1)) || ~isreal(r1(2)) || ~isreal(r1(3))
                        shdoa = 1;
                    end
                    x(l,k) = r1(1);
                    y(l,k) = r1(2);
                    z(l,k) = r1(3);
                end
            end
            surf(x, y, z,'LineStyle','none','FaceColor',"b",'Facealpha',0.03)
        end
    end
    radius_0 = radius_0*exp(c(i)*split_step*sim_step);
end

% plot reach set at terminal time
radius = radius_0;
exp_radius_rec_1(end) = radius;
center(:) = expr_rec(end,:);
const = [];
x_R0(:) = center(:);
for l = 1:Sub_sample_Num
    for m = 1:Sub_sample_Num
        cur_num = Sub_sample_Num*(l-1)+m;
        w_tmp = (1-2*rand(1,3)); w_tmp = w_tmp/norm(w_tmp);
        w_tmp = radius*w_tmp;
        x_traj(cur_num,:,1) = [center(1),center(2),center(3),w_tmp(1),w_tmp(2),w_tmp(3)];
        for n = 1:geo_sim_N
            Cur_r(:) = x_traj(cur_num,1:3,n);
            Cur_w(:) = x_traj(cur_num,4:6,n);
            [Next_r,Next_Ang_Vel] = Geodesic_Dynamic_Exp(Cur_r,Cur_w,tmp_Q,sim_step);
            x_traj(cur_num,1:3,n+1) = Next_r(:);
            x_traj(cur_num,4:6,n+1) = Next_Ang_Vel(:);
        end
        % plot3(squeeze(x_traj(cur_num,1,:)),squeeze(x_traj(cur_num,2,:)),squeeze(x_traj(cur_num,3,:)),"b");
        x_R(:) = x_traj(cur_num,1:3,end);
        const = [const,(x_R-x_R0)'*exp_Q_R*(x_R-x_R0)<=radius_tmp*radius_tmp];
    end
end
Objective_exp_Q = -trace(exp_Q_R);
sol_Q_exp = optimize(const,Objective_exp_Q,ops);
opt_Exp_Q = value(exp_Q_R);
if sol_Q_exp.problem == 0 && Check_Matrix_Definiteness(opt_Exp_Q) && abs(opt_Exp_Q(1,1))> 0.01 && abs(opt_Exp_Q(2,2))> 0.01 && abs(opt_Exp_Q(3,3))> 0.01
    [V,D] = eig(opt_Exp_Q);
    A_plot = sqrt(D)*V';
    for l = 1:num_mb
        init_x = radius_tmp - (l-1)*2*radius_tmp/(num_mb-1);
        r = real(sqrt(radius_tmp^2-init_x^2));
        for k = 1:num_mb
            theta = 2*pi*(k-1)/(num_mb-1);
            init_z = r*sin(theta);
            init_y = r*cos(theta);
            init_r = [init_x;...
                init_y;...
                init_z];
            r1 = (A_plot\init_r) + center;
            if ~isreal(r1(1)) || ~isreal(r1(2)) || ~isreal(r1(3))
                shdoa = 1;
            end
            x(l,k) = r1(1);
            y(l,k) = r1(2);
            z(l,k) = r1(3);
        end
    end
    surf(x, y, z,'LineStyle','none','FaceColor',"b",'Facealpha',0.03)
end


zlabel({'r_{03}'}); ylabel({'r_{02}'}); xlabel({'r_{01}'});
legend3 = legend([h2 h3],{'sampled initial states','terminal states'});
set(legend3,'Location','northeast');
view([-81,11.4]);

%% get the constraints of R
function const = const_W(x_lb,x_ub,Q,P,J,c1)
% Reference: C. Fan, J. Kapinski, X. Jin, and S. Mitra, "Locally optimal
% reach set over-approximation for nonlinear systems,” in 2016
% International Conference on Embedded Software (EMSOFT), Oct. 2016, pp.
% 1–10. doi: 10.1145/2968478.2968482.

const = [];
for k1 = 1:2
    if k1 == 1
        w1 = x_lb(1);
    else
        w1 = x_ub(1);
    end
    for k2 = 1:2
        if k2 == 1
            w2 = x_lb(2);
        else
            w2 = x_ub(2);
        end
        for k3 = 1:2
            if k3 == 1
                w3 = x_lb(3);
            else
                w3 = x_ub(3);
            end
            % hat(w)
            F1 = [0,-w3,w2;...
                  w3,0,-w1;...
                  -w2,w1,0];

            W = [F1*Q-Q*F1-2*c1*Q,Q;...
                 Q,J'*P+P*J-2*c1*P];
            const = [const,W <= 0];
        end
    end

end
end



