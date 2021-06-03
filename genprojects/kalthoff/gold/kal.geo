//+
e = 0.07;
E = 10;
eps=0.0001;
Point(1) = {0, 0, 0, E};
//+
Point(2) = {100, 0, 0, E}; //Middle Line
//+
Point(3) = {100, 100, 0, e};
//+
Point(4) = {0, 100, 0, E};
//+
Point(5) = {0, 25+eps, 0, E};
//+
Point(6) = {50, 25, 0, e};
//+
Point(7) = {0, 25-eps,0,E};
Point(8) = {58, 0, 0, e};

Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 1};
//+
Line(8) = {6,2};
//+
Line(9) = {6,3};
//+
Line(10) = {6,8};
Curve Loop(1) = {4, 5, 6, 7, 1, 2, 3};
//+
Plane Surface(1) = {1};
Point{8} In Surface{1};
Line{8} In Surface{1};
//Line{9} In Surface{1};
//Line{10} In Surface{1};
//+
Physical Curve("load") = {7};
//+
Physical Curve("bottom") = {1};
//+
Physical Surface("domain") = {1};
