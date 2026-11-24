// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

// Unit square with named physical groups.
// Generate the mesh with: gmsh unitsquare-v4-physicalnames.geo -2 -format msh41

Point(1) = {0,0,0,0.5}; Point(2) = {1,0,0,0.5}; Point(3) = {1,1,0,0.5}; Point(4) = {0,1,0,0.5};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Curve Loop(1) = {1,2,3,4}; Plane Surface(1) = {1};
Physical Curve("boundary") = {1,2,3,4};
Physical Surface("my domain") = {1};
