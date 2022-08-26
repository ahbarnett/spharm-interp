# spharm-interp

MATLAB+Fortran routines for interpolation from scattered points on the sphere via spherical harmonics.

Authors: Alex Barnett (MATLAB); Leslie Greengard, Zydrunas Gimbutas, Marina Spivak (Fortran). 2015-2016. Repackaged 8/25/22.

The main task is to use CG to solve the normal equations for the least-squares problem of matching a spherical harmonic expansion of given degree to given data at given scattered points on the sphere. This expansion may then be evaluated on a regular (theta,phi) tensor-product grid on the sphere, or at arbitrary new points. There are also some helper and visualization routines. The tests are not formal unit tests, so the user has to parse the output. The drivers are in MATLAB, with a MEX library calling Fortran90 routines.


### Installation:

Edit `makefile` for your system, then `make`.
From MATLAB try `lsqsolvespharm` to test (takes a few seconds).


### Main routines available from MATLAB:

`lsqsolvespharm` : iterative LSQ solve of sph harm coeffs to match data at arbitrary scattered points on sphere  
`spharmeval` : evaluate spherical harmonic expansion at arbitrary sphere points (not performance code)  
`spharmgrideval` : evaluate sph harm expansion on grid (pure MATLAB version)  
`spharmgridevalf` : evaluate sph harm expansion on grid (MEX interface version)  
`spharmproj` : project grid data on the sphere onto spherical harmonic coeffs (MEX interface)  
`spharmprojfunc` : same as spharmproj but acts on function handle  
`testall.m` : run all tests  

Info about some other files:

`spharm.mw` : mWrap-annotated MATLAB used to generate gateway.c MEX interface  
`gateway.c` : MEX interface
`*.f` : Fortran sources
