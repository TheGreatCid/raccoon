// Gmsh project created on Thu May 27 11:33:28 2021
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {.75, 0, 0, 1.0};
//+
Point(3) = {.75, 2, 0, 1.0};
//+
Point(4) = {0, 2, 0, 1.0};
//+
Point(5) = {0, 1.00000000001, 0, 1.0};
//+
Point(6) = {.25, 1, 0, 1.0};
//+
Point(7) = {0, 0.999999999, 0, 1.0};
//+
Line(1) = {4, 5};
//+
Line(2) = {5, 6};
//+
Line(3) = {6, 7};
//+
Line(4) = {7, 1};
//+
Line(5) = {1, 2};
//+
Line(6) = {2, 3};
//+
Line(7) = {3, 4};
//+
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7};
//+
Plane Surface(1) = {1};