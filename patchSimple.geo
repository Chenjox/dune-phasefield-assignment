// Gmsh project created on Sun Jan 21 16:21:53 2024

width = 1;
height = 1;

// Zuerst die Punkte
Point(1) = {-width,-height,0};
Point(2) = {width,-height,0};
Point(3) = {width,height,0};
Point(4) = {-width,height,0};

// Dann die äußeren Linien
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

// unteres Element
Curve Loop(100) = {1,2,3,4};


Surface(100) = {100};

Physical Curve(4)   = {4};
Physical Curve(2)  = {2};
Physical Curve(1) = {1};
Physical Curve(3)    = {3};

Physical Surface(1) = {100};


// // Viele Linien, und vorallem Knoten
Transfinite Curve {1}  = 3 Using Progression 1;
Transfinite Curve {2}  = 3 Using Progression 1;
Transfinite Curve {3}  = 3 Using Progression 1;
Transfinite Curve {4}  = 3 Using Progression 1;
Transfinite Surface{100};

Coherence;
// SetOrder 1;
Recombine Surface {100};
Mesh 2;

Mesh.MshFileVersion = 2.0;
RenumberMeshNodes;
Save StrCat("patchA.msh");