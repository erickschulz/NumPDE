// Warehouse room with three shelves and two doors
// Room: [0, 10] x [0, 6]
// Door 1: bottom wall, x in [0.5, 2.0]
// Door 2: top wall, x in [7.0, 9.0]
// Three shelves evenly spaced, each 1.0 wide x 3.5 tall

lc = 0.3;  // mesh element size

// Room corners
Point(1)  = {0,   0, 0, lc};
Point(2)  = {0.5, 0, 0, lc};  // door1 start
Point(3)  = {2.0, 0, 0, lc};  // door1 end
Point(4)  = {10,  0, 0, lc};
Point(5)  = {10,  6, 0, lc};
Point(6)  = {9.0, 6, 0, lc};  // door2 end
Point(7)  = {7.0, 6, 0, lc};  // door2 start
Point(8)  = {0,   6, 0, lc};

// Shelf 1: [2.0, 3.0] x [1.25, 4.75]
Point(101) = {2.0, 1.25, 0, lc};
Point(102) = {3.0, 1.25, 0, lc};
Point(103) = {3.0, 4.75, 0, lc};
Point(104) = {2.0, 4.75, 0, lc};

// Shelf 2: [4.5, 5.5] x [1.25, 4.75]
Point(201) = {4.5, 1.25, 0, lc};
Point(202) = {5.5, 1.25, 0, lc};
Point(203) = {5.5, 4.75, 0, lc};
Point(204) = {4.5, 4.75, 0, lc};

// Shelf 3: [7.0, 8.0] x [1.25, 4.75]
Point(301) = {7.0, 1.25, 0, lc};
Point(302) = {8.0, 1.25, 0, lc};
Point(303) = {8.0, 4.75, 0, lc};
Point(304) = {7.0, 4.75, 0, lc};

// Room boundary lines
Line(1) = {1, 2};    // wall
Line(2) = {2, 3};    // door1
Line(3) = {3, 4};    // wall
Line(4) = {4, 5};    // wall
Line(5) = {5, 6};    // wall
Line(6) = {6, 7};    // door2
Line(7) = {7, 8};    // wall
Line(8) = {8, 1};    // wall

// Shelf 1
Line(101) = {101, 102};
Line(102) = {102, 103};
Line(103) = {103, 104};
Line(104) = {104, 101};

// Shelf 2
Line(201) = {201, 202};
Line(202) = {202, 203};
Line(203) = {203, 204};
Line(204) = {204, 201};

// Shelf 3
Line(301) = {301, 302};
Line(302) = {302, 303};
Line(303) = {303, 304};
Line(304) = {304, 301};

// Curve loops
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Curve Loop(2) = {101, 102, 103, 104};
Curve Loop(3) = {201, 202, 203, 204};
Curve Loop(4) = {301, 302, 303, 304};

// Room surface with shelf holes
Plane Surface(1) = {1, 2, 3, 4};

// Physical groups
Physical Curve("doors") = {2, 6};
Physical Curve("walls") = {1, 3, 4, 5, 7, 8,
                           101, 102, 103, 104,
                           201, 202, 203, 204,
                           301, 302, 303, 304};
Physical Surface("room") = {1};
