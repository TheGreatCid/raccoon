//+
Mesh.MshFileVersion = 2.2;
e = 150;

Point(1) = {0,0,0};
Point(2) = {0.5,0,0};
Point(3) = {0.5,1.5,0};
Point(4) = {0,1.5,0};

Point(5) = {0,0.75,0};
Point(6) = {0.25, 0.75,0};
Point(7) = {0,0.75,0};

Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,4};
Line(9) = {5,6};
Line(10) = {7,6};
Line(11) = {4,5};
Line(12) = {7,1};

Curve Loop(1) = {5,6,7,11,9,10,12};
Plane Surface(1) = {1};


Physical Curve("left") = {8};
Physical Curve("top") = {7};
Physical Curve("right") = {6};
Physical Curve("bottom") = {5};
Physical Surface("body") = {1};
Physical Surface("Domain") = {1};
