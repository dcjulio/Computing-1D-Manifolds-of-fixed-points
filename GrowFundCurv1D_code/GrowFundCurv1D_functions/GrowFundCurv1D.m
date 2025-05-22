function Manif=GrowFundCurv1D(opts,eigenval_idx)
% Input:
% - opts: all the options needed to compute the manifold
% - eigenval_idx: optional argument. It is which eigenvector to consider for 2D manifolds
if nargin < 2
    eigenval_idx = [];
else
    if eigenval_idx~=1 && eigenval_idx~=2
        fprintf('\n ERROR! The optional second argument eigenval_idx can be either 1 or 2\n')
        return;
    end
end
    
    %% Initializing the manifold structure
    manif = init_manif(opts);

    tic
    %%
    %saves the information of the system stored in the structure 'manif' and 'opts'
    system_info=manif.system_info;
    thesystem=opts.thesystem;

    %% Warnings in case the specifications are not correct for the computation of this manifold

    if strcmp(manif.orientability,'orientation-reversing')
        fprintf('\n----  Manifold is non-orientable. We will proceed computing both branches\n');
        mapiter=2;
    else
        mapiter=1;
    end



    %% Warnings if the manifold we want to compute is not 1D
    
    if manif.grow_info.dimension>1
        fprintf('\n----  The dimension of the manifold is %i and this algorithm is designed to compute one-dimensional manifolds. -----\n',manif.grow_info.dimension)
        fprintf(' To choose automatically which of the two eigenvalues to use, include a second argument to GrowFundCurv1D(opts,#) where #: 1 or 2\n\n') 
        fprintf('The eigenvectors are:\n')

        for i=1:manif.grow_info.dimension
            fprintf('(%i) Eigenvalue: %.3f, Eigenvector: (%.3f, %.3f, %.3f)\n',i, manif.grow_info.eigval(i), manif.grow_info.eigvec(i,1),manif.grow_info.eigvec(i,2),manif.grow_info.eigvec(i,3))
        end

        if ~isempty(eigenval_idx)
            x = eigenval_idx;
            fprintf('\nChoice: eigenval_idx = %i\n\n', eigenval_idx)
        else
            fprintf('\n Which eigenvector do you want to use?') 
            prompt = "\n... Press 1 or 2\n... Press 0 in case you want to exit \n\n";
            x = input(prompt);
        end

        if x == 0
            return
        else
            manif.grow_info.eigval=manif.grow_info.eigval(x);
            manif.grow_info.eigvec=manif.grow_info.eigvec(x,:);
        end

    end


%% Printing the general information of this run

    names=fieldnames(system_info.par);
    
    fprintf('\n----');
    for k=1:length(names)
        fprintf(' %s:%0.2f ',names{k},system_info.par.(names{k}));
    end
    fprintf('----');

    if strcmp(manif.orientability,'orientation-preserving')
        fprintf('\n----  Computing W%s(%s), %s branch up to arclength %.0f ---- ',lower(opts.stability(1)),opts.name_fixpoint,manif.name(end-2:end), opts.user_arclength);
    else
        fprintf('\n----  Computing W%s(%s), both branches up to arclength %.0f each ---- ',lower(opts.stability(1)),opts.name_fixpoint,opts.user_arclength);
    end
    
    fprintf('\n----  Total maximum iterations of the fundamental domain: %i ---- ',manif.grow_info.max_funditer);
    fprintf('\n\n')
    fprintf('----- Acc. Conditions ----- ')
    fprintf('\n| AlphaMax:%.2e',manif.grow_info.alphamax)
    fprintf('\n| DeltaAlphaMax:%.2e ',manif.grow_info.deltalphamax)
    fprintf('\n| DeltaMin:%.2e ' ,manif.grow_info.deltamin)
    fprintf('\n| DeltaMax:%.2e \n' ,manif.grow_info.deltamax)
    fprintf('---------- \n')
    
    
    %-------------------------------------------------------------------
    %-------------------------------------------------------------------
    % initializing the fields for the initial information
    manif.grow_info.runinf.rem_deltamin=0; %how many points are removed because of deltamin
    manif.grow_info.runinf.rem_nan=0; %how many points are removed because of NaN values
    manif.grow_info.runinf.rem_inf=0; %how many points are removed at infinity because of duplication
    manif.grow_info.runinf.add_alphamax=0; %how many points are added because of alpha_max
    manif.grow_info.runinf.add_deltamax=0;%how many points are added because of delta_max
    manif.grow_info.runinf.add_deltalphamax=0;  %how many points are added because of (delta alpha)_max

%% computes first fundamental domain

% by default, we choose the eigenvector going to positive x
if manif.grow_info.eigvec(1)<0 % is it is going to negative x, then consider the other side
    manif.grow_info.eigvec=-manif.grow_info.eigvec;
end

% Choosing the initial distance correctly:
% at what distance from the fix point the fundamental domain starts so 
% the end of the fundamental domain is approx at distance: opts.accpar.init_step
% we assume the linearlize system to do the calculation

% if orientation preserving we need just the first iteration of the map and
% it is scaled by the eigenvalue
if strcmp(manif.orientability,'orientation-preserving')
    if strcmp(manif.stability,'Smanifold')
        dist = manif.grow_info.init_step*abs(manif.grow_info.eigval);
    elseif strcmp(manif.stability,'Umanifold')
        dist = manif.grow_info.init_step/abs(manif.grow_info.eigval);
    end
end

% if orientation-reversing it is scaled by the eigenval^(mapiter+1)
  if strcmp(manif.orientability,'orientation-reversing') 
    if strcmp(manif.stability,'Smanifold')
        dist = manif.grow_info.init_step*abs(manif.grow_info.eigval^(mapiter+1));
    elseif strcmp(manif.stability,'Umanifold')
        dist = manif.grow_info.init_step/abs(manif.grow_info.eigval^(mapiter+1));
    end

  end

% starting point of the fundamental domain
fund_init.x=manif.fixp.x+dist*manif.grow_info.eigvec(1);
fund_init.y=manif.fixp.y+dist*manif.grow_info.eigvec(2);
fund_init.z=manif.fixp.z+dist*manif.grow_info.eigvec(3);

% mapping the fundamental domain

fund_end = thesystem.mapping(fund_init,manif.stability,opts,mapiter);


% interpolate N=4 points linearly
fund.points.x = linspace(fund_init.x, fund_end.x, 6);
fund.points.y = linspace(fund_init.y, fund_end.y, 6);
fund.points.z = linspace(fund_init.z, fund_end.z, 6);


fund.points.arc = arclength(fund.points);
fund.points.idx_fund_dom = [1 numel(fund.points.x)];

if strcmp(manif.orientability,'orientation-preserving') 

    manif.points=fund.points;

    total_arc = fund.points.arc(end);

elseif strcmp(manif.orientability,'orientation-reversing') 

    manif.pointspos=fund.points;

    total_arc_pos = fund.points.arc(end);

    manif.pointsneg.x=[]; manif.pointsneg.y=[]; manif.pointsneg.z=[]; manif.pointsneg.arc=[]; manif.pointsneg.idx_fund_dom = []; 
    total_arc_neg = 0;
end


fund_dom_arc = fund.points.arc(end);

%%
fprintf(' Starting fundamental domain arclength  %f\n',fund.points.arc(end))

Manif=manif;

%----------- Starting the loop
iter=0;
stop_arc=0; % flag to stop iterating if we reached the desired arclength
stop_arc_neg = 0; % flag to stop iterating if we reached the desired arclength (neg branch)
stop_arc_pos = 0; % flag to stop iterating if we reached the desired arclength (pos branch)

% iter: number of iterations of the fundamental domain
while iter < manif.grow_info.max_funditer && stop_arc ~= 1 && (stop_arc_neg + stop_arc_pos) ~= 2
      % we haven't reached the max number of iterations
      % and, we haven't reached the arclength (orientation-preserving case)
      % and, we haven't reached the arclength for the two branches of the manifolds (orientation-reversing case)

        iter=iter+1;

        %mapping the points
        mappoints = thesystem.mapping(fund.points,manif.stability,opts);


    %% STARTING THE ALGORITHM

    %----------- Removing points at infinity or NaN values. This is only because issues from compactification


    %Delete NaN values
    nan_idx=union(union(find(isnan(mappoints.x)),find(isnan(mappoints.y))),find(isnan(mappoints.z)));
    Manif.grow_info.runinf.rem_nan=Manif.grow_info.runinf.rem_nan+numel(nan_idx);

      if numel(nan_idx)==numel(mappoints.x)
        fprintf('\n\n ALL THE POINTS HAVE BEEN MAPPED TO NAN VALUES. LAST POINT IS AT (%f,%f,%f)\n\n',fund.points.x(end),fund.points.y(end),fund.points.z(end))

        break;
       end
    
    mappoints.x(nan_idx)=[]; 
    mappoints.y(nan_idx)=[]; 
    mappoints.z(nan_idx)=[]; 
    
    fund.points.x(nan_idx)=[]; 
    fund.points.y(nan_idx)=[]; 
    fund.points.z(nan_idx)=[]; 

    %Delete points at infinity
    inf_idx=union(find((mappoints.x.^2+mappoints.y.^2)==1),find(abs(mappoints.z)==1));
    Manif.grow_info.runinf.rem_inf=Manif.grow_info.runinf.rem_inf+numel(inf_idx);



    if numel(inf_idx)==numel(mappoints.x)
        fprintf('\n\n THE MANIFOLD HAS REACHED INFINITY AT (%f,%f,%f)\n\n',fund.points.x(end),fund.points.y(end),fund.points.z(end))
        break;
    end



    mappoints.x(inf_idx)=[]; 
    mappoints.y(inf_idx)=[]; 
    mappoints.z(inf_idx)=[]; 
    
    fund.points.x(inf_idx)=[]; 
    fund.points.y(inf_idx)=[]; 
    fund.points.z(inf_idx)=[]; 

    
    fprintf('\n ITERATION OF THE FUNDAMENTAL DOMAIN %i\n',iter)

    
    arc_mappoints = arclength(mappoints); %arclength of the mapped fundamental domain

    % Check how much arclength is left to compute
    if strcmp(manif.orientability,'orientation-preserving')
        needed_arc = opts.user_arclength - total_arc;
    elseif strcmp(manif.orientability,'orientation-reversing')
        if mod(iter,2)==0 
            needed_arc = opts.user_arclength - total_arc_pos;
        else
            needed_arc = opts.user_arclength - total_arc_neg;
        end
    end
    
    
    % if in this iteration we exceed the desired arclength, then we chop
    % the fundamental domain up to the desired arclength
    if arc_mappoints(end) > needed_arc 
        idx_arc = find(arc_mappoints > needed_arc,1); %where we exceed the extra needed arc
        % chop the fund domain and the mappoints up to there
        fund.points.x=fund.points.x(1:idx_arc);
        fund.points.y=fund.points.y(1:idx_arc);
        fund.points.z=fund.points.z(1:idx_arc);

        mappoints.x=mappoints.x(1:idx_arc);
        mappoints.y=mappoints.y(1:idx_arc);
        mappoints.z=mappoints.z(1:idx_arc);
        
        
        if strcmp(manif.orientability,'orientation-preserving')
            stop_arc=1; % this is the last iteration
        elseif strcmp(manif.orientability,'orientation-reversing')
            if mod(iter,2)==0 % we are in the positive branch
                stop_arc_pos=1;
            else
                stop_arc_neg=1;
            end
        end
    end

%------------------------------------------------------- 
%---%----------- Adding points depending on Acc. Cond.  
    
    %initializing the structures to add points
    add_acc=struct(); 
    newpoints=struct();
    mapnewpoints=struct();
    
    
	add_acc.iter=0;
    add_acc.failed=[]; % points that failed acc cond in last loop
    add_acc.loop=true; % still doing the while loop 


    % Interpolation mesh
    fund_initial=fund.points; %starting mesh of the fundamental domain
    t_initial=0:1/(numel(fund_initial.x)-1):1; % parametrization for meshpoints

%---%-------------- Loop of the same mesh checking acc cond (this adds points)
    while add_acc.loop 
%---%--------------
        add_acc.loop=false; %to stop while loop % if at least one point is added this turns true
        add_acc.iter=add_acc.iter+1; %counter of iterations


        % Interpolating points from previous fundamental domain.
         if add_acc.iter==1
            tt=t_initial;
            t_interp=tt(1:end-1)+(tt(2:end)-tt(1:end-1))/2; % parametrization of interpolated points
            interp = makima3D(fund_initial,t_initial,t_interp); % compute interpolated preimage
            mapinterp = thesystem.mapping(interp,manif.stability,opts); % interpolated image
         else
             tt=sort([tt t_interp(add_acc.add)]); %parametrization of (new) mesh points
             t_interp=tt(1:end-1)+(tt(2:end)-tt(1:end-1))/2; % parametrization of (new) interpolated points
             interp = makima3D(fund_initial,t_initial,t_interp); % compute interpolated preimage
             mapinterp = thesystem.mapping(interp,manif.stability,opts); % interpolated image
         end

%%
        add_acc.add=[]; %points we are going to add

        % idx of points to check acc cond
        if add_acc.iter==1
            for_idx=2:(numel(mappoints.x)-1);
        else
            for_idx=add_acc.failed;
        end      

        fprintf('  loop number %i (points to check %i...',add_acc.iter,numel(for_idx));

%-------%---------- Going through the points that failed
        million=0;
        for k=for_idx
%-------%----------      
            % a flag for when # million points have been checked
            million=million+1;

            if floor(million/1000000)==ceil(million/1000000) %is integer?
                fprintf(' -checkpoint %i million points checked-...',floor(million/1000000));
            end

                
            % coordinates of mapped points
            add_acc.p0=[mappoints.x(k-1), mappoints.y(k-1), mappoints.z(k-1)];
            add_acc.p1=[mappoints.x(k), mappoints.y(k), mappoints.z(k)]; % the point we are actually looking at
            add_acc.p2=[mappoints.x(k+1), mappoints.y(k+1), mappoints.z(k+1)];
            
            % Distance btw points
            add_acc.delta0=norm(add_acc.p1-add_acc.p0); % before
            add_acc.delta2=norm(add_acc.p1-add_acc.p2); % after
            add_acc.alpha = angles(add_acc.p0,add_acc.p1,add_acc.p2); % angle btw points
            
            % points btw p0p1 and p1p2 in the interpolated points
            add_acc.p0_new=[mapinterp.x(k-1), mapinterp.y(k-1), mapinterp.z(k-1)];
            add_acc.p2_new=[mapinterp.x(k), mapinterp.y(k), mapinterp.z(k)];
    

%-----------%------ Adding points


%-----------%------ If delta > deltamax 

            %------ If it is the second point in the mesh
            %------ Check the first delta and add a point before
            if k==2 && add_acc.delta0>manif.grow_info.deltamax  

                %add point p01
                add_acc.add    =[add_acc.add k-1]; %idx of the point we are going to add
                add_acc.loop=true; % We have to check if we need to put more points in the mesh
                Manif.grow_info.runinf.add_deltamax=Manif.grow_info.runinf.add_deltamax+1;
                
            %------ Check the second delta and add a point after 
            elseif add_acc.delta2>manif.grow_info.deltamax %k>2 && add_acc.delta2>manif.grow_info.deltamax
                %add point p12
                add_acc.add    =[add_acc.add k]; %idx of the point we are going to add
                add_acc.loop=true; %idx of the point we are going to add
                Manif.grow_info.runinf.add_deltamax=Manif.grow_info.runinf.add_deltamax+1;
           
                
                
%-----------%------ If alpha > alphamax  or   Delta*alpha > Delta*alpha max

            %------ If alpha fails or BOTH Delta*alpha fail
            %------ Choose where to add a point.
            elseif add_acc.alpha>=manif.grow_info.alphamax || (add_acc.delta0*add_acc.alpha>=manif.grow_info.deltalphamax && add_acc.delta2*add_acc.alpha>=manif.grow_info.deltalphamax) %1 

            %------ Add point where Delta>Deltamin
                
                %-- If only Delta0>Deltamin
                %-- Add a point btw p0 and p1
                if add_acc.delta0>manif.grow_info.deltamin && add_acc.delta2<manif.grow_info.deltamin 
                    % Add a point if we didnt added the point in the previous acc cond checks
                    if numel(add_acc.add)==0 || add_acc.add(end)~=k-1
                        add_acc.add    =[add_acc.add k-1];  %idx of the point we are going to add
                        add_acc.loop=true; %idx of the point we are going to add
                        if (add_acc.delta0*add_acc.alpha>=manif.grow_info.deltalphamax && add_acc.delta2*add_acc.alpha>=manif.grow_info.deltalphamax) 
                            Manif.grow_info.runinf.add_deltalphamax=Manif.grow_info.runinf.add_deltalphamax+1;
                        else
                            Manif.grow_info.runinf.add_alphamax=Manif.grow_info.runinf.add_alphamax+1;
                        end
                    end
                end
                
                %-- If only Delta2>Deltamin
                %-- Add a point btw p1 and p2
                if add_acc.delta2>manif.grow_info.deltamin && add_acc.delta0<manif.grow_info.deltamin
                    %add point p12
                    add_acc.add    =[add_acc.add k]; %idx of the point we are going to add
                    add_acc.loop=true; %idx of the point we are going to add
                    if (add_acc.delta0*add_acc.alpha>=manif.grow_info.deltalphamax && add_acc.delta2*add_acc.alpha>=manif.grow_info.deltalphamax) 
                        Manif.grow_info.runinf.add_deltalphamax=Manif.grow_info.runinf.add_deltalphamax+1;
                    else
                        Manif.grow_info.runinf.add_alphamax=Manif.grow_info.runinf.add_alphamax+1;
                    end
                end
                
                %-- If both Delta0 and Delta 2 > Deltamin
                %-- Choose where to add point
                if add_acc.delta2>manif.grow_info.deltamin && add_acc.delta0>manif.grow_info.deltamin
                    add_acc.alpha0_new = angles(add_acc.p0_new,add_acc.p1,add_acc.p2); % angle btw points
                    add_acc.alpha2_new = angles(add_acc.p0,add_acc.p1,add_acc.p2_new); 

                    if add_acc.alpha0_new < add_acc.alpha2_new 
                        % Add a point if we didnt added the point in the previous acc cond checks
                        if numel(add_acc.add)==0 || add_acc.add(end)~=k-1
                            add_acc.add    =[add_acc.add k-1];  %idx of the point we are going to add
                            add_acc.loop=true; %idx of the point we are going to add
                            if (add_acc.delta0*add_acc.alpha>=manif.grow_info.deltalphamax && add_acc.delta2*add_acc.alpha>=manif.grow_info.deltalphamax) 
                                Manif.grow_info.runinf.add_deltalphamax=Manif.grow_info.runinf.add_deltalphamax+1;
                            else
                                Manif.grow_info.runinf.add_alphamax=Manif.grow_info.runinf.add_alphamax+1;
                            end
                        end
                    else
                        %add point p12
                        add_acc.add    =[add_acc.add k]; %idx of the point we are going to add
                        add_acc.loop=true; %idx of the point we are going to add
                        if (add_acc.delta0*add_acc.alpha>=manif.grow_info.deltalphamax && add_acc.delta2*add_acc.alpha>=manif.grow_info.deltalphamax) 
                            Manif.grow_info.runinf.add_deltalphamax=Manif.grow_info.runinf.add_deltalphamax+1;
                        else
                            Manif.grow_info.runinf.add_alphamax=Manif.grow_info.runinf.add_alphamax+1;
                        end
                    end
                end
                
                
                
%-----------%------ If just one Delta*alpha > Delta*alpha max ( and alpha < alphamax (previous elseif is when alpha > alphamax )

            %------ Delta0*alpha > max, and Delta0 > Deltamin
            %------ Add a point btw p0 and p1
            elseif add_acc.delta0*add_acc.alpha>=manif.grow_info.deltalphamax && add_acc.delta0>manif.grow_info.deltamin
                % Add a point if we didnt added the point in the previous acc cond checks
                if numel(add_acc.add)==0 || add_acc.add(end)~=k-1
                    add_acc.add    =[add_acc.add k-1];  %idx of the point we are going to add
                    add_acc.loop=true; %idx of the point we are going to add
                    Manif.grow_info.runinf.add_deltalphamax=Manif.grow_info.runinf.add_deltalphamax+1;
                end
                
            %------ Delta2*alpha > max, and Delta2 > Deltamin
            %------ Add a point btw p0 and p1
            elseif add_acc.delta2*add_acc.alpha>=manif.grow_info.deltalphamax && add_acc.delta2>manif.grow_info.deltamin
                %add point p12
                add_acc.add    =[add_acc.add k]; %idx of the point we are going to add
                add_acc.loop=true; %idx of the point we are going to add
                Manif.grow_info.runinf.add_deltalphamax=Manif.grow_info.runinf.add_deltalphamax+1;
                
                
%-----------%------              
            end    % (if loop) Adding points     
%-----------%------      
%-------%---------- 
        end       % (for loop) Going through the points that failed 
%-------%----------

        newpoints.x=interp.x(add_acc.add);
        newpoints.y=interp.y(add_acc.add);
        newpoints.z=interp.z(add_acc.add);
        
        
        mapnewpoints.x=mapinterp.x(add_acc.add);
        mapnewpoints.y=mapinterp.y(add_acc.add);
        mapnewpoints.z=mapinterp.z(add_acc.add);

        
       fprintf(' added points: %i) \n',numel(add_acc.add));

        
        % Add the points 
        if ~isempty(add_acc.add)
            
            % get updated idx of failed points
            plus=0:(numel(add_acc.add)-1);
            add_acc.failed = unique(sort([add_acc.add+plus add_acc.add+plus+1 add_acc.add+plus+2]));
            add_acc.failed=add_acc.failed(add_acc.failed>1);
            add_acc.failed=add_acc.failed(add_acc.failed<numel(mappoints.x)+numel(add_acc.add));

            % add points in the mapped manifold and in the old manifold
            mappoints.x=insert(mappoints.x,mapnewpoints.x,add_acc.add);
            mappoints.y=insert(mappoints.y,mapnewpoints.y,add_acc.add);
            mappoints.z=insert(mappoints.z,mapnewpoints.z,add_acc.add); 

            fund.points.x=insert(fund.points.x,newpoints.x,add_acc.add);
            fund.points.y=insert(fund.points.y,newpoints.y,add_acc.add);
            fund.points.z=insert(fund.points.z,newpoints.z,add_acc.add); 


        end
        
%---%--------------       
    end           % (while loop) Checking acc cond 
%---%--------------

%---%----------- Final section: save info

    fund.points.x=mappoints.x;
    fund.points.y=mappoints.y;
    fund.points.z=mappoints.z;
    fund.points.arc = arclength(fund.points);

    fund_dom_arc = [fund_dom_arc fund.points.arc(end)]; % how the arclegth of the fundamental domains increase

    %fprintf('  rate %.2f \n', growth_rate_fund(end));
    fprintf(' Fund domain arclength %.1f \n', fund_dom_arc(end));

    %add new branch to the entire manifold
    if strcmp(manif.orientability,'orientation-preserving') 

        % add indices of fundamental domain
        if stop_arc ~= 1
            iter_fund=numel(Manif.points.idx_fund_dom)/2 + 1;
            Manif.points.idx_fund_dom(iter_fund,:)=numel(Manif.points.x)+[1 numel(fund.points.x)];
            fund.points.idx_fund_dom = Manif.points.idx_fund_dom(iter_fund,:);
        end
        

        Manif.points.x=[Manif.points.x fund.points.x];
        Manif.points.y=[Manif.points.y fund.points.y];
        Manif.points.z=[Manif.points.z fund.points.z];
        Manif.points.arc = arclength(Manif.points);

        total_arc = Manif.points.arc(end);
        fprintf(' Total arclength %.0f \n', total_arc);
        
        % Check once again if we reached the desired arclength
        if total_arc > opts.user_arclength 
            stop_arc = 1; %this is the last iteration
        end

    elseif strcmp(manif.orientability,'orientation-reversing')
        if mod(iter,2)==0

            % add indices of fundamental domain
            if stop_arc_pos ~= 1
                iter_fund=numel(Manif.pointspos.idx_fund_dom)/2 + 1;
                Manif.pointspos.idx_fund_dom(iter_fund,:)=numel(Manif.pointspos.x)+[1 numel(fund.points.x)];
                fund.points.idx_fund_dom = Manif.pointspos.idx_fund_dom(iter_fund,:);
            end
            
            
            Manif.pointspos.x=[Manif.pointspos.x fund.points.x];
            Manif.pointspos.y=[Manif.pointspos.y fund.points.y];
            Manif.pointspos.z=[Manif.pointspos.z fund.points.z];
            Manif.pointspos.arc = arclength(Manif.pointspos);

            total_arc_pos = Manif.pointspos.arc(end);
            fprintf(' Total arclength positive branch %.1f \n', total_arc_pos);
     
            % Check once again if we reached the desired arclength
            if total_arc_pos > opts.user_arclength 
                stop_arc_pos = 1;
            end

        else

            % add indices of fundamental domain
            if stop_arc_neg ~= 1
                iter_fund=numel(Manif.pointsneg.idx_fund_dom)/2 + 1;
                Manif.pointsneg.idx_fund_dom(iter_fund,:)=numel(Manif.pointsneg.x)+[1 numel(fund.points.x)];
                fund.points.idx_fund_dom = Manif.pointsneg.idx_fund_dom(iter_fund,:);
            end
            
            Manif.pointsneg.x=[Manif.pointsneg.x fund.points.x];
            Manif.pointsneg.y=[Manif.pointsneg.y fund.points.y];
            Manif.pointsneg.z=[Manif.pointsneg.z fund.points.z];
            Manif.pointsneg.arc = arclength(Manif.pointsneg);


            total_arc_neg = Manif.pointsneg.arc(end);
            fprintf(' Total arclength negative branch %.1f \n', total_arc_neg);

            % Check once again if we reached the desired arclength
            if total_arc_neg > opts.user_arclength 
                stop_arc_neg = 1;
            end

        end
    end

%---%----------- END Final section: save info
end
%--------------- END adding points

Manif.grow_info.runinf.time=toc;

fprintf('\n elapsed time is %.3f seconds\n\n',Manif.grow_info.runinf.time)
if strcmp(manif.orientability,'orientation-preserving')
    fprintf(' %i final points, arclength %.1f \n',numel(Manif.points.x),Manif.points.arc(end)) %76800 longer
else
    fprintf(' %i final points (pos), arclength (pos) %.1f \n',numel(Manif.pointspos.x),Manif.pointspos.arc(end)) %76800 longer
    fprintf(' %i final points (neg), arclength (neg) %.1f \n',numel(Manif.pointsneg.x),Manif.pointsneg.arc(end)) %76800 longer
end
fprintf('   * %i points removed \n',Manif.grow_info.runinf.rem_deltamin+Manif.grow_info.runinf.rem_nan+Manif.grow_info.runinf.rem_inf) 
fprintf('   * %i points added from deltamax \n',Manif.grow_info.runinf.add_deltamax) 
fprintf('   * %i points added from alpha \n',Manif.grow_info.runinf.add_alphamax) 
fprintf('   * %i points added from delta*alpha \n',Manif.grow_info.runinf.add_deltalphamax) 
%% FUNCTIONS

function interp = makima3D(points,t,tt)
% get interpolation points
% t: parametrization of the points
% tt: parametrization of the interpolated points
interp=struct();

interp.x = interp1(t,points.x,tt,'makima','extrap');
interp.y = interp1(t,points.y,tt,'makima','extrap');
interp.z = interp1(t,points.z,tt,'makima','extrap');
end 

%----------------

function arclen = arclength(points)
%arclength between each point of a vector (px,py,pz)
arclen=((points.x(1:end-1)-points.x(2:end)).^2 + (points.y(1:end-1)-points.y(2:end)).^2 + (points.z(1:end-1)-points.z(2:end)).^2).^(1/2);
arclen=[0 cumsum(arclen)];
end % function arclength

%----------------

function alpha=angles(p0,p1,p2)
%angle between p0p1 and p1p2
n1 = (p1 - p2) / norm(p1 - p2);  % Normalized vectors
n2 = (p0 - p1) / norm(p0 - p1);
alpha = atan2(norm(cross(n1, n2)), dot(n1, n2)); %gives value from 0 to pi
end

%----------------

function Anew=insert(A,B,ind)
% Anew: new vector with the new values
% A: Old vector
% B: vector with new values
%ind: index to insert values after this row


    % Preallocate output
    Anew = zeros(1,numel(A)+numel(B));

    % Find indices for old data
    addRows = ismember(1:numel(A), ind);
    oldDataInd = (1:numel(A)) + cumsum([0, addRows(1:end-1)]);

    % Add in old data
    Anew(oldDataInd) = A;

    % Find indices for new data
    newDataInd = (1:length(ind)) + ind;

    % Add in new data
    Anew(newDataInd) = B;
end

end



