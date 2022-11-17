// Gmsh project created on Thu Oct 13 12:35:15 2022
SetFactory("OpenCASCADE");

S = 0.001;
H = S*1.0;
L = S*12.0;
L2 = S*10.0;
b = S*0.1;

h = S*0.1;
l = S*3.0;

bi = 10;
bj = 80;
bj2 = 200;
ui = 140;
vi = 100;
vy = 20;
E = L/20.0;

Point(1) = {0, -H+b, 0, 1.0};
Point(2) = {0, H-b, 0, 1.0};
Point(3) = {L, H-b, 0, 1.0};
Point(4) = {L, -H+b, 0, 1.0};
Point(5) = {0, H, 0, 1.0};
Point(6) = {L, H, 0, 1.0};
Point(7) = {0, -H, 0, 1.0};
Point(8) = {L, -H, 0, 1.0};
Point(9) = {(L-l)/2, -h, 0, 1.0};
Point(10) = {(L-l)/2, h, 0, 1.0};
Point(11) = {(L+l)/2, h, 0, 1.0};
Point(12) = {(L+l)/2, -h, 0, 1.0};
Point(13) = {(L-l)/2-b, -h-b, 0, 1.0};
Point(14) = {(L-l)/2-b, h+b, 0, 1.0};
Point(15) = {(L+l)/2+b, h+b, 0, 1.0};
Point(16) = {(L+l)/2+b, -h-b, 0, 1.0};
Point(17) = {-L2, -H+b, 0, 1.0};
Point(18) = {-L2, H-b, 0, 1.0};
Point(19) = {-L2, H, 0, 1.0};
Point(20) = {-L2, -H, 0, 1.0};
Point(21) = {L+L2, -H+b, 0, 1.0};
Point(22) = {L+L2, H-b, 0, 1.0};
Point(23) = {L+L2, H, 0, 1.0};
Point(24) = {L+L2, -H, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {2, 5};
Line(6) = {5, 6};
Line(7) = {6, 3};
Line(8) = {4, 8};
Line(9) = {8, 7};
Line(10) = {7, 1};
Line(11) = {9, 10};
Line(12) = {10, 11};
Line(13) = {11, 12};
Line(14) = {12, 9};
Line(15) = {13, 16};
Line(16) = {16, 15};
Line(17) = {15, 14};
Line(18) = {14, 13};
Line(19) = {9, 13};
Line(20) = {12, 16};
Line(21) = {11, 15};
Line(22) = {10, 14};
Line(23) = {2, 18};
Line(24) = {18, 17};
Line(25) = {17, 1};
Line(26) = {17, 20};
Line(27) = {20, 7};
Line(28) = {5, 19};
Line(29) = {19, 18};
Line(30) = {4, 21};
Line(31) = {21, 22};
Line(32) = {22, 3};
Line(33) = {8, 24};
Line(34) = {24, 21};
Line(35) = {22, 23};
Line(36) = {23, 6};

Curve Loop(1) = {4, 1, 2, 3};
Curve Loop(2) = {17, 18, 15, 16};
Curve Loop(3) = {2, -7, -6, -5};
Curve Loop(4) = {9, 10, -4, 8};
Curve Loop(5) = {15, -20, 14, 19};
Curve Loop(6) = {20, 16, -21, 13};
Curve Loop(7) = {12, 21, 17, -22};
Curve Loop(8) = {19, -18, -22, -11};
Curve Loop(9) = {24, 25, 1, 23};
Curve Loop(10) = {27, 10, -25, 26};
Curve Loop(11) = {23, -29, -28, -5};
Curve Loop(12) = {30, 31, 32, 3};
Curve Loop(13) = {33, 34, -30, 8};
Curve Loop(14) = {32, -7, -36, -35};

Plane Surface(1) = {1, 2};
Plane Surface(2) = {3};
Plane Surface(3) = {4};
Plane Surface(4) = {5};
Plane Surface(5) = {6};
Plane Surface(6) = {7};
Plane Surface(7) = {8};
Plane Surface(8) = {9};
Plane Surface(9) = {10};
Plane Surface(10) = {11};
Plane Surface(11) = {12};
Plane Surface(12) = {13};
Plane Surface(13) = {14};

Transfinite Curve {18, 11, 13, 16, 19, 22, 21, 20,29, 35, 26, 34, 10, 8, 5, 7, 29, 26, 32, 35, 34} = bi Using Progression 1;
Transfinite Curve {15, 14, 12, 17} = bj Using Progression 1;
Transfinite Curve {9, 4, 2, 6} = bj2 Using Progression 1;
Transfinite Curve {1, 3} = vy Using Progression 1;
Transfinite Curve {27, 25, 23, 28, 33, 30, 32, 36} = vi Using Progression 1;
Transfinite Curve {31, 24} = vy Using Progression 1;
Transfinite Surface {7};
Transfinite Surface {4};
Transfinite Surface {5};
Transfinite Surface {6};
Transfinite Surface {3};
Transfinite Surface {2};
Transfinite Surface {12};
Transfinite Surface {13};
Transfinite Surface {9};
Transfinite Surface {10};
Transfinite Surface {8};
Transfinite Surface {11};
Recombine Surface {2, 3, 6, 4, 5, 7, 9, 12, 13, 10, 8, 11};

Extrude {0, 0, E} {
    Surface{1};
    Surface{2};
    Surface{3};
    Surface{4};
    Surface{5};
    Surface{6};
    Surface{7};
    Surface{8};
    Surface{9};
    Surface{10};
    Surface{11};
    Surface{12};
    Surface{13};
    Layers {1};
    Recombine;
}


//+
Physical Surface("front", 97) = {46, 22, 56, 52, 26, 49, 30, 62, 59};
//+
Physical Surface("back", 98) = {11, 1, 8, 10, 2, 13, 12, 3, 9};
//+
Physical Surface("wall_top", 99) = {60, 24, 51};
//+
Physical Surface("wall_base", 100) = {47, 27, 57};
//+
Physical Surface("inlet", 101) = {43};
//+
Physical Surface("outlet", 102) = {54};
//+
Physical Surface("plate_top", 103) = {38};
//+
Physical Surface("plate_base", 104) = {32};
//+
Physical Surface("plate_trailing", 105) = {36};
//+
Physical Surface("plate_leading", 106) = {41};
//+
Physical Volume("fluid", 107) = {8, 10, 9, 11, 13, 12, 3, 2, 6, 4, 1, 5, 7};
