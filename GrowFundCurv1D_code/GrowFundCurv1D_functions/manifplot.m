function h=manifplot(manif)

angle=str2num(manif.inter.angle(1:end-2))*pi;

%-- Colormap pmin all by arclength
if strcmp(manif.orientability,'orientation-preserving')
    rad=(manif.points.x.^2+manif.points.y.^2).^(1/2);
elseif  strcmp(manif.orientability,'orientation-reversing')
    radpos=(manif.pointspos.x.^2+manif.pointspos.y.^2).^(1/2);
    radneg=(manif.pointsneg.x.^2+manif.pointsneg.y.^2).^(1/2);
end

if strcmp(manif.stability,'Smanifold')
    RGB1=[131, 195, 251]/255; % light blue
    RGB2=[5, 52, 122]/255; % dark blue
elseif strcmp(manif.stability,'Umanifold')
    RGB1=[240, 120, 98]/255; % light red
    RGB2=[168, 25, 17]/255;
end
R= linspace(RGB1(1),RGB2(1),100);  %// Red from 212/255 to 0
G = linspace(RGB1(2),RGB2(2),100);   %// Green from 212/255 to 0
B = linspace(RGB1(3),RGB2(3),100);  %// Blue from 1 to 170/255
RGB = [R(:), G(:), B(:)];


h=figure;
hold on
colormap(RGB); 
if strcmp(manif.orientability,'orientation-preserving')
    color_line3(manif.points.x,manif.points.y,manif.points.z,rad,'EdgeAlpha',1,'LineWidth',1.5);
elseif  strcmp(manif.orientability,'orientation-reversing')
    color_line3(manif.pointspos.x,manif.pointspos.y,manif.pointspos.z,radpos,'EdgeAlpha',1,'LineWidth',1.5);
    color_line3(manif.pointsneg.x,manif.pointsneg.y,manif.pointsneg.z,radneg,'EdgeAlpha',1,'LineWidth',1.5);
end


%-- Plane
plane.x=[cos(angle),0];
plane.y=[sin(angle),0];
plane.z= [-1,1];
plane.color=[230, 178, 17]/255;
surf(repmat(plane.x,2,1), repmat(plane.y,2,1), repmat(plane.z,2,1)','FaceAlpha',0.4, 'EdgeColor',plane.color,'FaceColor',plane.color,'FaceLighting','gouraud','LineWidth',1.7)

%-- Unit circle
[xunit,yunit] = circle(0,0,1,1000);
plot3(xunit,yunit,ones(size(xunit)),'k','LineWidth',1.5)
plot3(xunit,yunit,-ones(size(xunit)),'k','LineWidth',1.5)
% 

% % %-- intersection points
if strcmp(manif.orientability,'orientation-preserving')
    plot3(manif.inter.points.x,manif.inter.points.y,manif.inter.points.z,'.','color',RGB2*0.7,'MarkerSize',11);
elseif  strcmp(manif.orientability,'orientation-reversing')
    plot3(manif.inter.pointspos.x,manif.inter.pointspos.y,manif.inter.pointspos.z,'.','color',RGB2*0.7,'MarkerSize',11);
    plot3(manif.inter.pointsneg.x,manif.inter.pointsneg.y,manif.inter.pointsneg.z,'.','color',RGB2*0.7,'MarkerSize',11);
end
%-- fixed points
plot3(manif.system_info.fixp.pplu.x, manif.system_info.fixp.pplu.y, manif.system_info.fixp.pplu.z,'marker','o','MarkerFaceColor',[230, 178, 17]/255,'MarkerEdgeColor',[87, 67, 6]/255,'LineWidth',1.4,'MarkerSize',6.5)
plot3(manif.system_info.fixp.pmin.x,manif.system_info.fixp.pmin.y,manif.system_info.fixp.pmin.z,'marker','o','MarkerFaceColor',[230, 178, 17]/255,'MarkerEdgeColor',[87, 67, 6]/255,'LineWidth',1.4,'MarkerSize',6.5)

xlabel('x')
ylabel('y')
zlabel('z')
xlim([-1.01 1.01])
ylim([-1.01 1.01])
zlim([-1.01 1.01])


title(strrep(manif.name,'_','\_'))

daspect([1 1 1])
view([100,30])
if strcmp(manif.orientability,'orientation-preserving')
    clim([min(rad), max(rad)]);
elseif  strcmp(manif.orientability,'orientation-reversing')
    clim([min([radpos,radneg]), max([radpos,radneg])]);  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = color_line3(x, y, z, c, varargin)
% color_line3 plots a 3-D "line" with c-data as color
%
%       h = color_line(x, y, z, c)
%       by default: 'LineStyle','-' and 'Marker','none'
%
%          or
%       h = color_line(x, y, z, c, mark) 
%          or
%       h = color_line(x, y, z, c, 'Property','value'...) 
%             with valid 'Property','value' pairs for a surface object
%
%  in:  x      x-data
%       y      y-data
%       z      z-data
%       c      4th dimension for colouring
%       mark   for scatter plots with no connecting line
%
% out:  h   handle of the surface object
h = surf(...
  'XData',[x(:) x(:)],...
  'YData',[y(:) y(:)],...
  'ZData',[z(:) z(:)],...
  'CData',[c(:) c(:)],...
  'FaceColor','none',...
  'EdgeColor','flat',...
  'Marker','none');
  
if nargin ==5
    switch varargin{1}
        case {'+' 'o' '*' '.' 'x' 'square' 'diamond' 'v' '^' '>' '<' 'pentagram' 'p' 'hexagram' 'h'}
            set(h,'LineStyle','none','Marker',varargin{1})
        otherwise
            error(['Invalid marker: ' varargin{1}])
    end
elseif nargin > 5
    set(h,varargin{:})
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xunit,yunit] = circle(x,y,r,n)
th = linspace(0,2*pi,n);
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
end

end