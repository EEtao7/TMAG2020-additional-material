function point = point(radius,angle)
% caculate the point by radius and angle

deg = pi/180;
point = [radius*cos(angle*deg) radius*sin(angle*deg)];

end