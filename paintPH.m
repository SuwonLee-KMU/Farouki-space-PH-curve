% Generated on: 190810
% Last modification: 190812
% Author: Suwon Lee from Seoul National University

classdef paintPH < handle
  properties
    spacePHobj
    figureNumber = 1;
  end

  methods
    function obj = paintPH(spacePHobj,figureNumber)
      obj.spacePHobj    = spacePHobj;
      obj.figureNumber  = figureNumber;
    end
    function set.spacePHobj(obj,value)
      obj.spacePHobj = value;
    end
    function fig = visualize(obj)
      fig = figure(obj.figureNumber);
      xi  = linspace(0,1,100);
      r   = obj.spacePHobj.curve(xi);
      CP  = obj.spacePHobj.controlPoints;
      p_i = obj.spacePHobj.initialPosition;
      p_f = obj.spacePHobj.finalPosition;
      t_i = obj.spacePHobj.initialUnitTangent;
      t_f = obj.spacePHobj.finalUnitTangent;
      plot3(r(:,1),r(:,2),r(:,3),'linewidth',2,'color','b');
      hold on;
      plot3(CP(:,1),CP(:,2),CP(:,3),'--o','color','k');
      quiver3(p_i(1),p_i(2),p_i(3),t_i(1),t_i(2),t_i(3),'color','r','autoscalefactor',0.2,'linewidth',1,'maxheadsize',0.5);
      quiver3(p_f(1),p_f(2),p_f(3),t_f(1),t_f(2),t_f(3),'color','r','autoscalefactor',0.2,'linewidth',1,'maxheadsize',0.5);
      scatter3(p_i(1),p_i(2),p_i(3),'markeredgecolor','k','marker','s','sizedata',200);
      scatter3(p_f(1),p_f(2),p_f(3),'markeredgecolor','k','marker','d','sizedata',200);
      % grid on; box on;
      axis equal;
      xlabel('x'); ylabel('y'); zlabel('z');
    end
  end
end
