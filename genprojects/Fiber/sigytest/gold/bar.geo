//+
Mesh.MshFileVersion = 2.2;
e = 150;

Point(1) = {0,0,0};
Point(2) = {0.5,0,0};
Point(3) = {0.5,1.5,0};
Point(4) = {0,1.5,0};

Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,4};
Line(8) = {4,1};

Curve Loop(1) = {5,6,7,8};
Plane Surface(1) = {1};

Transfinite Curve{5,6,7,8} = e Using Progression 1;
Transfinite Surface{1};

Recombine Surface{1};

Physical Curve("left") = {8};
Physical Curve("top") = {7};
Physical Curve("right") = {6};
Physical Curve("bottom") = {5};
Physical Surface("body") = {1};
Physical Surface("Domain") = {1};
