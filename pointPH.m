% Generated on: 190814
% Last modification: 190816
% Author: Suwon Lee from Seoul National University
%
% pointPH class evaluate the curve at certain point on the curve with the
% curve parameter xi (pointPH.spacePHparameter = xi;).
%
% User needs to hand over a spacePH class instance when allocating a new
% pointPH class instance. (ex) A = pointPH(spacePHobj);
% The default value of the curve parameter is set to zero.
%
% User can set the spacePHparameter as a scalar or vector array.
% When properties spacePHobject or spacePHparameter are set to new
% values, the curve is evaluated automatically and set the following
% properties.
% - positionVector  : point on the curve in the space
% - quaternionA     : quaternion representation of the curve
% - parametricSpeed : parametric speed of the curve, sigma
% - curvature       : curvature of the curve, kappa
% - torsion         : torsion of the curve, tau
% - FrenetBasis     : basis vectors of the Frenet frame, [t,p,b]
% 
% Warning: Since pointPH.spacePHobject is a handle class, the user should
% call pointPH.updateTransients() function manually to update the
% corresponding curve after modifying pointPH.spacePHobj, or the curve is
% not updated automatically.
classdef pointPH < handle
  properties (SetObservable)
    spacePHobject
    spacePHparameter = 0;
  end

  properties (Transient, SetAccess = protected)
    positionVector
    quaternionA
    parametricSpeed
    curvature
    torsion
    tangent
    FrenetBasis
  end

  methods % Methods users can use.
    function set.spacePHobject(obj,value)
      if isa(value,'spacePH')
        obj.spacePHobject = value;
      else
        error('invalid spacePHobj');
      end
    end
    function set.spacePHparameter(obj,value)
      obj.spacePHparameter = value(:);
    end
    function obj = updateTransients(obj)
      obj.positionVector  = curve(obj);
      obj.quaternionA     = getQuaternionAtPoint(obj);
      [obj.parametricSpeed, obj.tangent] = getParametricSpeed(obj);
      obj.curvature       = getCurvature(obj);
      obj.torsion         = getTorsion(obj);
      obj.FrenetBasis     = getFrenetBasis(obj);
    end
    function quaternionA = getQuaternionAtPoint(obj)
      xi = obj.spacePHparameter;
      A0 = obj.spacePHobject.quaternionA0;
      A1 = obj.spacePHobject.quaternionA1;
      A2 = obj.spacePHobject.quaternionA2;
      quaternionA = zeros(numel(xi),4);
      for i = 1:numel(xi)
        quaternionA(i,:) = A0*(1-xi(i))^2 + A1*2*(1-xi(i))*xi(i) + A2*xi(i)^2;
      end
    end
    function [parametricSpeed,tangentVector] = getParametricSpeed(obj)
      xi  = obj.spacePHparameter;
      nxi = numel(xi);
      A  = obj.getQuaternionAtPoint;
      tangentVector = zeros(nxi,3);
      % for i = 1:nxi
        % cA = spacePH.quatConj(A(i,:));
        % tangentVector(i,:) = spacePH.AiA(A(i,:),cA);
      % end
      tangentVector = obj.curveDerivative;
      parametricSpeed = vecnorm(tangentVector,2,2);
    end
    function curvature = getCurvature(obj)
      dr        = obj.curveDerivative;
      ddr       = obj.curvePthDerivative(2);
      num       = vecnorm(cross(dr,ddr,2),2,2);
      den       = vecnorm(dr,2,2).^3;
      curvature = num./den;
    end
    function torsion = getTorsion(obj)
      dr = obj.curveDerivative;
      ddr = obj.curvePthDerivative(2);
      d3r = obj.curvePthDerivative(3);
      num = dot(cross(dr,ddr),d3r,2);
      den = norm(cross(dr,ddr))^2;
      torsion = num./den;
    end
    function [t,p,b] = getFrenetBasis(obj)
      dr = obj.curveDerivative;
      ddr = obj.curvePthDerivative(2);
      t = dr./vecnorm(dr,2,2);
      p = cross(cross(dr,ddr,2),t,2)./vecnorm(cross(dr,ddr,2));
      b = cross(dr,ddr,2)./vecnorm(cross(dr,ddr,2),2,2);
    end
    function E_RMF = computeE_RMF(obj,npoints)
      A     = obj.spacePHobject;
      B     = pointPH(A);
      B.spacePHparameter = linspace(0,1,npoints);
      dxi   = [0;diff(B.spacePHparameter)];
      E_RMF = sum(B.curvature.^2.*B.parametricSpeed.*dxi);
    end
  end

  methods (Hidden) % Points on curve and its derivatives
    function r = curve(obj)
      xi  = obj.spacePHparameter; 
      nxi = numel(xi);
      r   = zeros(nxi,3);
      for j = 1:nxi
        for i = 1:6
          term = obj.spacePHobject.controlPoints(i,:)*pointPH.Bernstein(5,i-1,xi(j));
          r(j,:) = r(j,:) + term;
        end
      end
    end
    function dr = curveDerivative(obj)
      xi  = obj.spacePHparameter; 
      nxi = numel(xi);
      dr   = zeros(nxi,3);
      for j = 1:nxi
        for i = 1:6
          term = obj.spacePHobject.controlPoints(i,:)*pointPH.BernsteinDerivative(5,i-1,xi(j));
          dr(j,:) = dr(j,:) + term;
        end
      end
    end
    function dpr = curvePthDerivative(obj,p)
      xi  = obj.spacePHparameter; 
      nxi = numel(xi);
      dpr = zeros(nxi,3);
      for j = 1:nxi
        for i = 1:6
          term = obj.spacePHobject.controlPoints(i,:)*pointPH.BernsteinPthDerivative(5,i-1,p,xi(j));
          dpr(j,:) = dpr(j,:) + term;
        end
      end
    end
  end

  methods (Hidden)
    function obj = pointPH(spacePHobj)
      obj.spacePHobject = spacePHobj;
      obj.attachListner();
      obj.updateTransients;
    end
    function attachListner(obj)
      addlistener(obj,'spacePHobject','PostSet',@pointPH.propChange);
      addlistener(obj,'spacePHparameter','PostSet',@pointPH.propChange);
    end
  end

  methods (Hidden, Static)  % Bernstein polynomial and its derivatives
    function Bin = Bernstein(n,i,x)
      Bin = nchoosek(n,i)*x^i*(1-x)^(n-i);
    end
    function DBin = BernsteinDerivative(n,i,x)
      if i <= 0
        term1 = 0;
      else
        term1 = pointPH.Bernstein(n-1,i-1,x);
      end
      if n-1 < i
        term2 = 0;
      else
        term2 = pointPH.Bernstein(n-1,i,x);
      end
      DBin = n*(term1-term2);
    end
    function DpBin = BernsteinPthDerivative(n,i,p,x)
      DpBin_ = 0;
      for k = max(0,i+p-n):min(i,p)
        term = (-1)^(k+p)*nchoosek(p,k)*pointPH.Bernstein(n-p,i-k,x);
        DpBin_ = DpBin_ + term;
      end
      DpBin = factorial(n)/factorial(n-p)*DpBin_;
    end
  end

  methods (Static, Hidden)  % For event listner callback
    function propChange(metaProp,eventData)
       h = eventData.AffectedObject;
       h.updateTransients();
    end
  end
end