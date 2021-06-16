// the projectile is assumed to be 63 degrees
angle = (90-70)*Pi/180;
angle2 = (90-85)*Pi/180;

E = 5e-3;
e = 0.2e-3;
eps = 0.001e-3;

Point(1) = {0, 0, 0, E};
Point(2) = {0, 25e-3-eps, 0, E/2};
Point(3) = {50e-3, 25e-3, 0, e};
Point(4) = {0, 25e-3+eps, 0, E/2};
Point(5) = {0, 100e-3, 0, E};
Point(6) = {50e-3+75e-3*Tan(angle), 100e-3, 0, e};
Point(7) = {100e-3, 100e-3, 0, E};
Point(8) = {100e-3, 0, 0, E};
Point(9) ={50e-3+75e-3*Tan(angle2),0,0,e};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 1};

Line(10) = {3, 6};
Line(11) = {3, 9};

Line Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9};

Plane Surface(1) = {1};
Line{10} In Surface {1};
Line{11} In Surface{1};

Physical Surface("all") = {1};
Physical Line("load") = {1};
Physical Line("bottom") = {8,9};
