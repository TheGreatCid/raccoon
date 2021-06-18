Mesh.MshFileVersion = 2.2;
Mesh.Smoothing = 100;
//Mesh.Algorithm = 8; // Delaunay for quads
e = 0.006;
E = 0.05;
Point(1) = {-1.5,-1.5,0,E};
Point(2) = {-1.5,1.5,0,E};
Point(3) = {1.5,1.5,0,E};
Point(4) = {1.5,-1.5,0,E};
Point(5) = {0,0,0,e};
Point(6) = {0,0.5,0,e};
Point(7) = {0,-0.5,0,e};
Point(8) = {0.5,0,0,e};
Point(9) = {-0.5,0,0,e};
Circle(1) = {9,5,6}; //Sectioned inner cirlce
Circle(2) = {6,5,8};
Circle(3) = {8,5,7};
Circle(4) = {7,5,9};

Line(10) = {1,4}; //on surface 1  //Outer Lines
Line(11) = {4,3}; //''2
Line(12) = {3,2}; //''3
Line(13) = {2,1}; //''4



//+
Curve Loop(1) = {13, 10, 11, 12};
//+
Curve Loop(2) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1, 2};

//Mesh Refinement

Point(10) = {0,1.47,0,e};
Point(11) = {0,-1.45,0,e};
Point(12) = {-0.6,0,0,e};
Point(13) = {0.6,0,0,e};
Point(14) = {0,-1,0,e};
Point(15) = {0,0.8,0,E/2};
Point(16) = {0,0.75,0,E/2};

//Line(14) = {6,10};
Line(15) = {7,11};
Line(16) = {10,12};
//Line(17) = {12,14};
Line(18) = {10,13};
//Line(19) = {13,14};
Line(20) = {15,16};

Point{10,11,12,13,14,15,16} In Surface{1};
Line{15,16,18,20} In Surface{1};
//+
Physical Surface("domain") = {1};
//+
Physical Curve("Top") = {12};
//+
Physical Curve("Bottom") = {10};
//+
Physical Curve("Hole") = {1, 4, 3, 2};
