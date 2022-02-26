function [alpha1, alpha2, beta1, beta2] = PlaneOfSymmetry3(p1, p2)

x1 = p1(2); y1 = p1(1); z1 = p1(3); 
x2 = p2(2); y2 = p2(1); z2 = p1(3); 

alpha1 = (y2-y1)/(z1*y2 - z2*y1); 
beta1 = (z1-z2)/(z1*y2 - z2*y1); 
alpha2 = (x1*y2 - x2*y1)/(z1*y2 - z2*y1); 
beta2 = (z1*x2 - z2*x1)/(z1*y2 - z2*y1); 

end