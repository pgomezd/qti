  function str = jf_hct_arg(pn, cg, ig)
%|function str = jf_hct_arg(pn, cg, ig)
%|
%| name/value pairs needed for hct command line
%|
%| Copyright 2009-11-21, Jeff Fessler, University of Michigan

str = [
	sprintf(' nx %d ny %d nz %d', ig.nx, ig.ny, ig.nz) ...
	sprintf(' dx %g dy %g dz %g', ig.dx, ig.dy, ig.dz) ...
	sprintf(' offset_x %g offset_y %g offset_z %g', ...
		ig.offset_x, ig.offset_y, ig.offset_z) ...
	sprintf(' ns %d nt %d na %d', cg.ns, cg.nt, cg.na) ...
	sprintf(' ds %g dt %g', cg.ds, cg.dt) ...
	sprintf(' offset_s %g offset_t %g', cg.offset_s, cg.offset_t) ...
	sprintf(' dfs %g dso %g dsd %g', cg.dfs, cg.dso, cg.dsd) ...
	sprintf(' orbit %g orbit_start %g', cg.orbit, cg.orbit_start) ...
	sprintf(' pitch %g source_z0 %g', cg.pitch, cg.source_z0) ...
];
