// the projectile is assumed to be 63 degrees
angle = (90-70)*Pi/180;
angle2 = (90-85)*Pi/180;

E = 4;
e = 0.2;
eps = 0.001;

Point(1) = {0, 0, 0, E};
Point(2) = {0, 25-eps, 0, E};
Point(3) = {50, 25, 0, e};
Point(4) = {0, 25+eps, 0, E};
Point(5) = {0, 100, 0, E};
Point(7) = {100, 100, 0, E};
Point(8) = {100, 0, 0, E};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 7};
Line(7) = {7, 8};
Line(8) = {1, 8};

Curve Loop(1) = {4, 5, 7, -8, 1, 2, 3};

Plane Surface(1) = {1};

Physical Surface("all") = {1};
Physical Line("load") = {1};
Physical Line("bottom") = {8};
//+
Physical Curve("other") = {4, 5, 7};
//+
