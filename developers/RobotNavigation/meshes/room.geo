// Room with two doors and one shelf obstacle
// Room: [0, 10] x [0, 6]
// Door 1: bottom wall, x in [0.5, 2.0]
// Door 2: top wall, x in [7.0, 9.0]
// Shelf: rectangle [4, 7] x [1.5, 2.0]

lc = 0.4;  // mesh element size

// Room corners
Point(1)  = {0,  0, 0, lc};
Point(2)  = {0.5, 0, 0, lc};  // door1 start
Point(3)  = {2.0, 0, 0, lc};  // door1 end
Point(4)  = {10, 0, 0, lc};
Point(5)  = {10, 6, 0, lc};
Point(6)  = {9.0, 6, 0, lc};  // door2 end
Point(7)  = {7.0, 6, 0, lc};  // door2 start
Point(8)  = {0,  6, 0, lc};

// Shelf obstacle corners
Point(9)  = {4.0, 1.5, 0, lc};
Point(10) = {7.0, 1.5, 0, lc};
Point(11) = {7.0, 2.0, 0, lc};
Point(12) = {4.0, 2.0, 0, lc};

// Room boundary lines (bottom wall with door1 gap)
Line(1) = {1, 2};    // wall
Line(2) = {2, 3};    // door1
Line(3) = {3, 4};    // wall
// Right wall
Line(4) = {4, 5};    // wall
// Top wall with door2 gap
Line(5) = {5, 6};    // wall
Line(6) = {6, 7};    // door2
Line(7) = {7, 8};    // wall
// Left wall
Line(8) = {8, 1};    // wall

// Shelf boundary
Line(9)  = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 9};

// Room boundary (outer)
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
// Shelf hole
Curve Loop(2) = {9, 10, 11, 12};

// Room surface with shelf hole
Plane Surface(1) = {1, 2};

// Physical groups
Physical Curve("doors") = {2, 6};
Physical Curve("walls") = {1, 3, 4, 5, 7, 8, 9, 10, 11, 12};
Physical Surface("room") = {1};
