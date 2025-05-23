function orbit=eps_pseudo_orbit(manif,k)
    
        orbit.idx=k; % this is the last index of the eps pseudo orbit
        i=1;

        if strcmp(manif.orientability,'orientation-preserving')

        % find preimages. When the idx==0 then we have reached the starting fundamental domain
        while orbit.idx(end)~=0 
            i=i+1;
            orbit.idx(i)=manif.points.idx_preimages(orbit.idx(end));
        end
        orbit.idx=flip(orbit.idx(1:end-1));

        % Save eps-pseudo-orbit
        orbit.x=manif.points.x(orbit.idx);
        orbit.y=manif.points.y(orbit.idx);
        orbit.z=manif.points.z(orbit.idx);


        elseif strcmp(manif.orientability,'orientation-reversing')


            % find preimages. When the idx==0 then we have reached the starting fundamental domain
            while orbit.idx(end)~=0
                i=i+1;
                if orbit.idx(end)>0 %comes from the positive branch
                    orbit.idx(i)=manif.pointspos.idx_preimages(orbit.idx(end));
                else %comes from the negative branch
                    orbit.idx(i)=manif.pointsneg.idx_preimages(abs(orbit.idx(end)));
                end
            end
            orbit.idx=flip(orbit.idx(1:end-1));
    
            orbit.x=[];orbit.y=[];orbit.z=[];
            % Save eps-pseudo-orbit
            for i=1:numel(orbit.idx)
                if orbit.idx(i)>0 %comes from positive branch
                    orbit.x=[orbit.x manif.pointspos.x(orbit.idx(i))];
                    orbit.y=[orbit.y manif.pointspos.y(orbit.idx(i))];
                    orbit.z=[orbit.z manif.pointspos.z(orbit.idx(i))];
                else %comes from negative branch
                    orbit.x=[orbit.x manif.pointsneg.x(abs(orbit.idx(i)))];
                    orbit.y=[orbit.y manif.pointsneg.y(abs(orbit.idx(i)))];
                    orbit.z=[orbit.z manif.pointsneg.z(abs(orbit.idx(i)))];
                end        
            end
        end

end