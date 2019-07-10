function [ Ke ] = PlaneTri2Stiffness( ncoor, Em )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
x1=ncoor(1);
y1=ncoor(2);
x2=ncoor(3);
y2=ncoor(4);
x3=ncoor(5);
y3=ncoor(6);

A=1/2*(  (x2*y3-x3*y2)+(x3*y1-x1*y3)+ (x1*y2-x2*y1)   );

x13=x1-x3;
x21=x2-x1;
x32=x3-x2;
y31=y3-y1;
y23=y2-y3;
y12=y1-y2;


coeff=1/(4*A);

termL=[y23  0  x32;
       0  x32  y23;
       y31  0  x13;
       0  x13  y31;
       y12  0  x21;
       0  x21  y12;];
   
termR=[y23  0  y31  0  y12  0;
       0  x32  0  x13  0  x21;
       x32  y23  x13  y31  x21  y12];
   
Ke=coeff*termL*Em*termR;



end

