% Generated on: 190810
% Last modification: 190810
% Author: Suwon Lee from Seoul National University

classdef spacePH < handle
  properties
    initialPosition
    initialUnitTangent
    finalPosition
    finalUnitTangent
    desiredArcLength
    psi0 = pi;
    psi2 = 0;
  end
  properties (Transient, SetAccess=protected)
    spatialRotation
    scaleFactor
    quaternionA0
    quaternionA1
    quaternionA2
    controlPoints
  end

  methods
    function obj = set.initialPosition(obj,value)
      if numel(value)==3
        obj.initialPosition = value(:)';
      else
        error('initialPosition should be a vector with dimension 3');
      end
    end
    function obj = set.finalPosition(obj,value)
      if numel(value)==3
        obj.finalPosition = value(:)';
      else
        error('finalPosition should be a vector with dimension 3');
      end      
    end
    function obj = set.initialUnitTangent(obj,value)
      if numel(value)==3
        if abs(norm(value)-1)>eps
          warning('The vector t_i is normalized into a unit vector');
        end
        obj.initialUnitTangent = value(:)'/norm(value);
      else
        error('Invalid initialUnitTangent');
      end
    end
    function obj = set.finalUnitTangent(obj,value)
      if numel(value)==3
        if abs(norm(value)-1)>eps
          warning('The vector t_f is normalized into a unit vector');
        end
        obj.finalUnitTangent = value(:)'/norm(value);
      else
        error('Invalid finalUnitTangent');
      end
    end
    function obj = updateTransients(obj)
      [R,f] = spacePH.getRotAndScale(obj.initialPosition,obj.finalPosition);
      [tiC,tfC,dpC,SC] = spacePH.cvt2canon(obj.initialUnitTangent,obj.finalUnitTangent,obj.finalPosition-obj.initialPosition,obj.desiredArcLength,R,f);
      [tiP,tiA] = spacePH.getPolarAzimuth(tiC);
      [tfP,tfA] = spacePH.getPolarAzimuth(tfC);
      [c0,c1,c2] = spacePH.getCoeffs27(tiP,tfP,tiA,tfA,SC,obj.psi0,obj.psi2);
      w = spacePH.getws(c0,c1,c2);
      [A0,A1,A2] = spacePH.getQuaternions(tiC,tfC,dpC,w,obj.psi0,obj.psi2);
      controlPoints_canon = spacePH.getBezierControlPoints([0,0,0],A0,A1,A2);

      controlPoints_ = zeros(size(controlPoints_canon));
      for i = 1:6
        controlPoints_(i,:) = (1/f*R'*controlPoints_canon(i,:)')' + obj.initialPosition;
      end

      obj.quaternionA0 = A0;
      obj.quaternionA1 = A1;
      obj.quaternionA2 = A2;
      obj.controlPoints = controlPoints_;
    end
    function r = curve(obj,xi)
      r = spacePH.phCurve(obj.controlPoints,xi);
    end
  end

  methods (Hidden)
    function obj = spacePH(p_i,p_f,t_i,t_f,S)
      obj.initialPosition    = p_i;
      obj.finalPosition      = p_f;
      obj.initialUnitTangent = t_i;
      obj.finalUnitTangent   = t_f;
      minS = norm(obj.initialPosition-obj.finalPosition);
      obj.desiredArcLength = max(minS,S);
      [R,f] = spacePH.getRotAndScale(obj.initialPosition,obj.finalPosition);
      obj.spatialRotation = R;
      obj.scaleFactor     = f;
    end
  end
  
  methods (Static, Hidden)
    function [R,f] = getRotAndScale(pi,pf)
      f = 1/norm(pi-pf);
      dp = pf-pi;
      theta_z = atan2(dp(2),dp(3));
      theta_y = atan2(dp(3),norm(dp(1:2)));
      R3_theta_z = [cos(theta_z),-sin(theta_z),0;...
      sin(theta_z),cos(theta_z),0;...
      0,0,1];
      R2_theta_y = [cos(theta_y),0,-sin(theta_y);...
      0,1,0;
      sin(theta_y),0,cos(theta_y)];
      R = R3_theta_z'*R2_theta_y';
    end
    function [tiC,tfC,dpC,SC] = cvt2canon(ti,tf,dp,S,R,f)
      tiC = (R*ti(:))';
      tfC = (R*tf(:))';
      dpC = f*(R*dp(:))';
      SC  = f*S;
    end
    function [polar,azimuth] = getPolarAzimuth(tangentVector)
      polar   = acos(tangentVector(1));
      azimuth = atan2(tangentVector(3),tangentVector(2));
    end
    function [c0,c1,c2] = getCoeffs27(initialPolar,finalPolar,initialAzimuth,finalAzimuth,SC,psi0,psi2)
      ci = cos(.5*initialPolar);
      si = sin(.5*initialPolar);
      cf = cos(.5*finalPolar);
      sf = sin(.5*finalPolar);
      dphi = finalAzimuth - initialAzimuth;
      dpsi = psi2 - psi0;
      c0 = 36*(SC^2-1);
      c1 = 6*((SC-1)*ci*cf*cos(dphi+dpsi)+(SC+1)*si*sf*cos(dpsi)-3*SC)+9*(ci^2-si^2+cf^2-sf^2);
      c2 = 2*(ci^2*sf^2+si^2*cf^2)-4*ci*si*cf*sf*cos(dphi);
    end
    function [w] = getws(c0,c1,c2)
      w2 = roots([c2,c1,c0]);
      wsquared = min(w2);
      w = sqrt(wsquared);
    end
    function quatOut = AiA(A1,A2)
      quaternionAiA = spacePH.quatMult(spacePH.quatMult(A1,[0,1,0,0]),A2);
      quatOut = quaternionAiA(2:4);
    end
    function cA = quatConj(A)
      cA = zeros(1,4);
      cA(1) = A(1);
      cA(2) = -A(2);
      cA(3) = -A(3);
      cA(4) = -A(4);
    end
    function A = quatMult(A1,A2)
      A = zeros(1,4);
      a1 = A1(1); b1 = A1(2); c1 = A1(3); d1 = A1(4);
      a2 = A2(1); b2 = A2(2); c2 = A2(3); d2 = A2(4);
      A(1) = a1*a2-b1*b2-c1*c2-d1*d2;
      A(2) = a1*b2+b1*a2+c1*d2-d1*c2;
      A(3) = a1*c2-b1*d2+c1*a2+d1*b2;
      A(4) = a1*d2+b1*c2-c1*b2+d1*a2;
    end
    function [A0,A1,A2] = getQuaternions(tiC,tfC,dpC,w,psi0,psi2)
      [tip,tia] = spacePH.getPolarAzimuth(tiC);
      [tfp,tfa] = spacePH.getPolarAzimuth(tfC);
      ci = cos(.5*tip);
      si = sin(.5*tip);
      cf = cos(.5*tfp);
      sf = sin(.5*tfp);
      A01 = w*cos(psi0)*ci*cos(tia);
      A02 = w*ci*(cos(psi0)*sin(tia)+sin(psi0)*cos(tia));
      A03 = w*sin(psi0)*si;
      A04 = w*cos(psi0)*si;
      A0  = [A01,A02,A03,A04];
      A21 = w*cos(psi2)*cf*cos(tfa);
      A22 = w*cf*(cos(psi2)*sin(tfa)+sin(psi2)*cos(tfa));
      A23 = w*sin(psi2)*sf;
      A24 = w*cos(psi2)*sf;
      A2  = [A21,A22,A23,A24];

      cA0 = spacePH.quatConj(A0);
      cA2 = spacePH.quatConj(A2);
      
      d = 120*dpC - 15*w^2*(tiC+tfC) + 5*(spacePH.AiA(A0,cA2)+spacePH.AiA(A2,cA0));

      psi1 = 0;
      term = sqrt(norm(d))/4/norm([norm(d),0,0]+d)*spacePH.quatMult([0,norm(d),0,0]+[0,d],[cos(psi1),sin(psi1),0,0]);
      A1 = -3/4*(A0+A2) + term;
    end
    function [controlPoints_canon] = getBezierControlPoints(p0,A0,A1,A2)
      cA0 = spacePH.quatConj(A0);
      cA1 = spacePH.quatConj(A1);
      cA2 = spacePH.quatConj(A2);
      p1 = p0 + 1/5*spacePH.AiA(A0,cA0);
      p2 = p1 + 1/10*(spacePH.AiA(A0,cA1)+spacePH.AiA(A1,cA0));
      p3 = p2 + 1/30*(spacePH.AiA(A0,cA2)+4*spacePH.AiA(A1,cA1)+spacePH.AiA(A2,cA0));
      p4 = p3 + 1/10*(spacePH.AiA(A1,cA2)+spacePH.AiA(A2,cA1));
      p5 = p4 + 1/5*spacePH.AiA(A2,cA2);
      controlPoints_canon = [p0(:)';p1(:)';p2(:)';p3(:)';p4(:)';p5(:)'];
    end
    function r = phCurve(controlPoints,xi)
      nxi = numel(xi);
      r  = zeros(nxi,3);
      for j = 1:nxi
        for i = 1:6
          term = controlPoints(i,:)*nchoosek(5,i-1)*(1-xi(j))^(5-i+1)*xi(j)^(i-1);
          r(j,:) = r(j,:)+term;
        end
      end
    end
  end
end