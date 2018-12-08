//USER INPUT
#define Mode 'S'       // Choice S == Serial Mode and Choice P == Parallel mode
#define nxc 65         // Number of cells in North and South Boundary
#define nyc 65         // Number of cells in East  and West  Boundary
#define Re  100.0      // Number of cells in East  and West  Boundary
#define pRefCell 1     // Pressure Reference cell // Assume center of cell
#define MAXitr 20000     // Maximum number of Main Loop iterations// Time marching iterations
#define MAXnormux (1e-8)     // Desired normalized L2 norm for U-Velocity
#define MAXnormuy (1e-8)     // Desired normalized L2 norm for V-Velocity
#define MAXnormp (1e-8)      // Desired normalized L2 norm for Pressure


#define nvar 3         // u,v,p // dont change it // this in only 2D solver
