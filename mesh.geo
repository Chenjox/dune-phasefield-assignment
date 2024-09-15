
width = 2;
a = 0.3;
sepa = 0.0001;
numNodes = 16;
numNodesOuter = 24;
//num_nodes= If (!Exists(num_nodes)) 
// 3
//EndIf

// Zuerst die Punkte
Point(1) = {0,0,0};
Point(2) = {-a,-sepa,0};
Point(3) = {-width/2,-sepa,0};
Point(4) = {-width/2,-width/2,0};
Point(5) = {0,-width/2,0};
Point(6) = {width/2,-width/2,0};
Point(7) = {width/2,0,0};
Point(8) = {width/2,width/2,0};
Point(9) = {0.0,width/2,0};
Point(10) = {-width/2,width/2,0};
Point(11) = {-width/2,sepa,0};
// innerer Kreis
Point(12) = {-a,sepa,0};
Point(13) = {-a,a,0};
Point(14) = {0,a,0};
Point(15) = {a,a,0};
Point(16) = {a,0,0};
Point(17) = {a,-a,0};
Point(18) = {0,-a,0};
Point(19) = {-a,-a,0};

// Points to Lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};
Line(10) = {10,11};
Line(11) = {11,12};
Line(12) = {12,13};
Line(13) = {13,14};
Line(14) = {14,15};
Line(15) = {15,16};
Line(16) = {16,17};
Line(17) = {17,18};
Line(18) = {18,19};
// Now the non-hamiltonian path
Line(19) = {19,2};
// Inner connection
Line(20) = {1,12};
Line(21) = {1,14};
Line(22) = {1,16};
Line(23) = {1,18};
// outer connections
Line(24) = {13,10};
Line(25) = {14,9};
Line(26) = {15,8};
Line(27) = {16,7};
Line(28) = {17,6};
Line(29) = {18,5};
Line(30) = {19,4};

// Now, create the regions
// x- y- 
Curve Loop(100) = {1,-19,-18,-23};
Curve Loop(101) = {2,3,-30,19};
Curve Loop(102) = {30,4,-29,18};
// x+ y-
Curve Loop(103) = {23,-17,-16,-22};
Curve Loop(104) = {17,29,5,-28};
Curve Loop(105) = {16,28,6,-27};
// x+ y+
Curve Loop(106) = {22,-15,-14,-21};
Curve Loop(107) = {15,27,7,-26};
Curve Loop(108) = {14,26,8,-25};
// x- y+
Curve Loop(109) = {20,12,13,-21};
Curve Loop(110) = {-11,-10,-24,-12};
Curve Loop(111) = {13,25,9,-24};

// Now Surfaces
Surface(100) = {100};
Surface(101) = {101};
Surface(102) = {102};
Surface(103) = {103};
Surface(104) = {104};
Surface(105) = {105};
Surface(106) = {106};
Surface(107) = {107};
Surface(108) = {108};
Surface(109) = {109};
Surface(110) = {110};
Surface(111) = {111};

// Place the nodes
fineLines[] = {1,3,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22,23}; 
For cur In {0:#fineLines[]-1}
  Transfinite Curve{fineLines[cur]} = numNodes;
EndFor

finerLines[] = {11,24,25,26,27,28,29,30,2};
For cur In {0:#finerLines[]-1}
  Transfinite Curve{finerLines[cur]} = numNodesOuter;
EndFor

For sur In {100:111}
  Transfinite Surface{sur};
EndFor

// Export only boundaries and elements
bounds[] = {1,2,3,4,5,6,7,8,9,10,11,20};
For boun In {0:#bounds[]-1}
Physical Curve(boun) = {bounds[boun]};
EndFor

For sur In {100:111}
  Physical Surface(sur)={sur};
EndFor

Coherence;
// SetOrder 1;

For sur In {100:111}
  Recombine Surface {sur};
EndFor
Mesh 2;


Mesh.MshFileVersion = 2.0;
RenumberMeshNodes;
Save StrCat("mode1slit.msh");