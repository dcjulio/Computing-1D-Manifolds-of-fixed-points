

%% Adding the path
    clear all
    close all
    addpath('./GrowFundCurv1D_functions');

%% Options

    %--- Information of the system
    opts.thesystem=StdHenon3D; % What is the name of the system file
    opts.par=struct('a', 4.2,'b', 0.3, 'xi', 1.2); % The parameter values and names (has to match with the names defined in StdHenon3D)
    opts.user_arclength=100; % What is the approximate arclength of the manifold
    
    %--- Number of iterations used to compute the manifold
    opts.max_funditer=30; % how many times (max) the algorithm iterates the fundamental domain

    %--- Accuracy parameters (default)
    %opts.accpar.alphamax=0.3;
    %opts.accpar.deltalphamax=0.001; 
    %opts.accpar.deltamin=0.000001;
    %opts.accpar.deltamax=0.01;  

    %--- Initial step (default)
    %opts.accpar.init_step=10^-7; % pos value is positive branch and viceversa


%% Computing the manifold: Ws(pmin) orientation-preserving
    opts.name_fixpoint='pmin'; % name of the fixed point associated to the manifold (has to match with the names defined in StdHenon3D)
    opts.stability='Smanifold';

    opts.branch='pos'; %which branch (Optional. Could be done with the sign of the init_step instead)
    pmin=GrowFundCurv1D(opts);

%% Computing intersection points
    angle=pi/2; %the angle of the plane from [-pi, pi]. (angle=pi/2: x==0 (y>0), angle=0: y==0 (x>0))
    pmin=inter_plane(pmin,angle);

%% Epsilon pseudo orbit (orientation-preserving)
    k=pmin.inter.points.idx(3); % pseudo orbit that finish at the point k
    orbit_pmin=eps_pseudo_orbit(pmin,k);

    % plot the epsilon orbit. Starting and end point are colored in solid red.
    f1=manifplot(pmin);
    hold on
    plot3(orbit_pmin.x([1,end]),orbit_pmin.y([1,end]),orbit_pmin.z([1,end]),'r.','MarkerSize',27) %epsilon orbit
    plot3(orbit_pmin.x,orbit_pmin.y,orbit_pmin.z,'ko--','LineWidth',1.2) %epsilon orbit
    
 %% Computing the manifold: Ws(pplu) orientation-reversing
    opts.user_arclength=opts.user_arclength/2; % each branch has half the total arclength
    opts.name_fixpoint='pplu'; % name of the fixed point associated to the manifold (has to match with the names defined in StdHenon3D)
    opts.stability='Smanifold';

    opts.branch='pos'; %which branch (Optional. Could be done with the sign of the init_step instead)
    pplu=GrowFundCurv1D(opts);

%% Computing intersection points
    angle=pi/2; %the angle of the plane from [-pi, pi]. (angle=pi/2: x==0 (y>0), angle=0: y==0 (x>0))
    pplu=inter_plane(pplu,angle);

%% Epsilon pseudo orbit (orientation-preserving)
    k=pplu.inter.pointspos.idx(3); % pseudo orbit that finish at the point k
    orbit_pplu=eps_pseudo_orbit(pplu,k);

    % plot the epsilon orbit. Starting and end point are colored in solid red.
    f2=manifplot(pplu);
    hold on
    hold on
    plot3(orbit_pplu.x([1,end]),orbit_pplu.y([1,end]),orbit_pplu.z([1,end]),'r.','MarkerSize',27) %epsilon orbit
    plot3(orbit_pplu.x,orbit_pplu.y,orbit_pplu.z,'ko--','LineWidth',1.2) %epsilon orbit
    