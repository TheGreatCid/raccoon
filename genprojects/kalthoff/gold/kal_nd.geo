// the projectile is assumed to be 63 degrees
angle = (90-70)*Pi/180;
angle2 = (90-85)*Pi/180;

E = 10;
e = 0.071;
eps = 0.001;

Point(1) = {0, 0, 0, E};
Point(2) = {0, 25-eps, 0, E};
Point(3) = {50, 25, 0, e};
Point(4) = {0, 25+eps, 0, E};
Point(5) = {0, 100, 0, E};
Point(7) = {100, 100, 0, E};
Point(8) = {100, 0, 0, E};

Point(9) = {50, 28, 0, e};
Point(10) = {50, 22 , 0, e};
Point(11) = {50, 25.5, 0 ,e};
Point(12) = {50, 24.5, 0 ,e};
Point(13) = {49.5, 25+.01, 0, e};
Point(14) = {49.5, 25-.01,0,e};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 7};
Line(7) = {7, 8};
Line(8) = {1, 8};

Circle(9) = {10,3,9};
Line(10) = {9,11};
Line(11) = {10,12};
Circle(12) = {11, 3, 13};
Circle(13) = {14, 3, 12};

Curve Loop(1) = {4, 5, 7, -8, 1, 2, 3};

Plane Surface(1) = {1};

Point{9} In Surface{1};
Point{10} In Surface{1};
Point{11} In Surface{1};
Curve{9} In Surface{1};
Line{10} In Surface{1};
Line{11} In Surface{1};
Curve{12} In Surface{1};
Curve{13} In Surface{1};


Physical Surface("all") = {1};
Physical Line("load") = {1};
Physical Line("bottom") = {8};
//+
Physical Curve("other") = {4, 5, 7};
//+
//+
Physical Point("hairpoint") = {3};
