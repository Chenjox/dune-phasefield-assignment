// Gmsh project created on Sun Jan 21 16:21:53 2024

width = 1;
height = 1;
patchoffset = (width/2)*0.8;
patchheightset = (height/2)*0.6;

// Zuerst die Punkte
Point(1) = {-width,-height,0};
Point(2) = {width,-height,0};
Point(3) = {width,height,0};
Point(4) = {-width,height,0};
Point(5) = {-width+patchoffset,-height+patchheightset,0};
Point(6) = {width-patchheightset,-height+patchoffset,0};
Point(7) = {width-patchoffset,height-patchheightset,0};
Point(8) = {-width+patchheightset,height-patchoffset,0};


// Dann die äußeren Linien
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
// diagonalen
Line(5) = {1,5};
Line(6) = {2,6};
Line(7) = {3,7};
Line(8) = {4,8};
// innere
Line(9) = {5,6};
Line(10) = {6,7};
Line(11) = {7,8};
Line(12) = {8,5};

// unteres Element
Curve Loop(100) = {1,6,-9,-5};
Curve Loop(101) = {2,7,-10,-6};
Curve Loop(102) = {3,8,-11,-7};
Curve Loop(103) = {4,5,-12,-8};

// Element der Mitte
Curve Loop(104) = {9,10,11,12};

Surface(100) = {100};
Surface(101) = {101};
Surface(102) = {102};
Surface(103) = {103};
Surface(104) = {104};

Physical Curve(4)   = {4};
Physical Curve(2)  = {2};
Physical Curve(1) = {1};
Physical Curve(3)    = {3};

Physical Surface("Domain",1) = {100,101,102,103,104};


// // Viele Linien, und vorallem Knoten
Transfinite Line {1,6,-9,-5}  = 2 Using Progression 1;
Transfinite Line {2,7,-10,-6}  = 2 Using Progression 1;
Transfinite Line {3,8,-11,-7}  = 2 Using Progression 1;
Transfinite Line {4,5,-12,-8}  = 2 Using Progression 1;
Transfinite Line {9,10,11,12}  = 2 Using Progression 1;


Coherence;
// SetOrder 1;
Recombine Surface {100,101,102,103,104};
Mesh 2;

Mesh.MshFileVersion = 2.0;
RenumberMeshNodes;
Save StrCat("patchC.msh");