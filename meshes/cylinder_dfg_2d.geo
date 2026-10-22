// DFG 2D-2, Schaefer & Turek (1996), Table 1. Use -setnumber coarse 1 for tests and
// -setnumber refine r to divide every production size by r: level 1 is r = 1 (this .msh),
// level 2 is r = 2.2 (_level2.msh), level 3 is r = 4 (_level3.msh).
If (!Exists(coarse))
  coarse = 0;
EndIf
If (!Exists(refine))
  refine = 1;
EndIf
If (coarse)
  near = 0.05; bulk = 0.22; wake = 0.16;
Else
  near = 0.004 / refine; bulk = 0.06 / refine; wake = 0.025 / refine;
EndIf
Point(1) = {0, 0, 0, bulk}; Point(2) = {2.2, 0, 0, bulk};
Point(3) = {2.2, 0.41, 0, bulk}; Point(4) = {0, 0.41, 0, bulk};
Point(5) = {0.2, 0.2, 0, near};
Point(6) = {0.25, 0.2, 0, near}; Point(7) = {0.2, 0.25, 0, near};
Point(8) = {0.15, 0.2, 0, near}; Point(9) = {0.2, 0.15, 0, near};
Line(1) = {1, 2}; Line(2) = {2, 3}; Line(3) = {3, 4}; Line(4) = {4, 1};
Circle(5) = {6, 5, 7}; Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9}; Circle(8) = {9, 5, 6};
Curve Loop(1) = {1, 2, 3, 4}; Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(1) = {1, 2};
Physical Curve("inlet", 11) = {4}; Physical Curve("outlet", 12) = {2};
Physical Curve("walls", 13) = {1, 3}; Physical Curve("cylinder", 14) = {5, 6, 7, 8};
Physical Surface("fluid", 20) = {1};
Field[1] = Distance; Field[1].CurvesList = {5, 6, 7, 8}; Field[1].Sampling = 100;
Field[2] = Threshold; Field[2].InField = 1;
Field[2].SizeMin = near; Field[2].SizeMax = bulk;
Field[2].DistMin = 0.015; Field[2].DistMax = 0.3;
Field[3] = Box; Field[3].VIn = wake; Field[3].VOut = bulk;
Field[3].XMin = 0.25; Field[3].XMax = 1.5;
Field[3].YMin = 0.1; Field[3].YMax = 0.3;
Field[4] = Min; Field[4].FieldsList = {2, 3}; Background Field = 4;
Mesh.Algorithm = 6;
Mesh.MshFileVersion = 2.2;
