// global length scale
lc = 5e-2;

// domain size
domain_width = 5;
domain_height = 5;

// pipe buried depth
pipe_depth = 2.5;

// pipe diamter and radius
pipe_d = 10e-2;
pipe_r = pipe_d / 2;

// domain corner point
Point(1) = {0, 0, 0, lc};
Point(2) = {domain_width, 0, 0, lc};
Point(3) = {domain_width, domain_height, 0, lc};
Point(4) = {0, domain_height, 0, lc};

// arc center
Point(5) = {0, pipe_depth, 0, lc};
// arc upper point
Point(6) = {0, pipe_depth + pipe_r, 0, lc};
// arc lower point
Point(7) = {0, pipe_depth - pipe_r, 0, lc};

// domain bottom line
Line(1) = {1, 2};
// domain right line
Line(2) = {2, 3};
// domain upper line
Line(3) = {3, 4};
// domain left upper line
Line(4) = {4, 6};
// doamin left lower line
Line(5) = {7, 1};

// half circle
Circle(6) = {7, 5, 6};

// create a surface from the lines
Curve Loop(1) = {4, -6, 5, 1, 2, 3};
Plane Surface(1) = {1};

// make matrix elements conforming to lowdim domain
Line {6} In Surface {1};

// add phsical group so that gmsh generates the target mesh
Physical Surface("flow_domain", 7) = {1};
Physical Curve("lowdim", 8) = {6};
