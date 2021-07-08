// the projectile is assumed to be 63 degrees
angle = (90-70)*Pi/180;
angle2 = (90-85)*Pi/180;

E = 8;
e = 0.8;
eps = 0.001;

Point(1) = {0, 0, 0, E};
Point(2) = {0, 25-eps, 0, E/10};
Point(3) = {50, 25, 0, e};
Point(4) = {0, 25+eps, 0, E/10};
Point(5) = {0, 100, 0, E};
Point(6) = {100, 26, 0, e};
Point(7) = {100, 100, 0, E};
Point(8) = {100, 0, 0, E};
Point(9) ={100,5,0,E/10};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 7};
Line(6) = {7, 6};
Line(8) = {6, 9};
Line(9) = {1, 8};
Line(12) = {9,8};

Line(10) = {3, 6};
Line(11) = {3, 9};

Curve Loop(1) = {4, 5, 6, 8, 12, -9, 1, 2, 3};
//+
Plane Surface(1) = {1};

//+
Line{10} In Surface {1};
Line{11} In Surface{1};

Physical Surface("all") = {1};
//+
Physical Curve("load") = {1};
//+
Physical Curve("bottom") = {9};
//+
Physical Curve("other") = {4, 5, 6};
//+
