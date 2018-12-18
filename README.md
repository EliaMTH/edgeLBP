File list: 

- docs_and_readmes:
	-edgeLBP_cite_me.bib;
	-output_data_readme.txt
	-sphere_edge_inters.pdf
- edgeLBP.m;
- edgeLBPp.m;
- README.txt; 



Description of the package:

This is an implementation of the edgeLBP as presented in the paper "Description and retrieval of geometric patterns on surface meshes using an edge-based LBP approach". You can find the complete reference in the docs_and_readmes folder. 

The only function you need to run the code is edgeLBP.m. No additional toolboxes are required. The code is not fully optimized. I'm working on a better structured version of the code, eventually I will update this archive with that version. 
This code was tested on Matlab 2014a and it should work even for previous versions. You can also use Octave, but I never tested it. 
If you have the Parallel Computing Toolbox, you can use the edgeLBPp.m function, which is faster than the classic edgeLBP.m.

Both functions are fully commented. If you have any problem and/or questions, please contact me at elia.moscoso@ge.imati.cnr.it.
If you use this code, or its variations, in your work, PLEASE cite us, using the bibtex reference in the doc_and_readmes folder.


Copyright (c) Moscoso~Thompson Elia, 2018


FAQ:
- Q: I'm using your code on mymesh.off, and I get an error. Why?
A: Usually, this is due to the mesh geometry. The mesh must be a SINGLE SHELL, with or without boundaries doesn't matter and must be a MANIFOLD surface. Another possibility is that you are trying to use a radius (R) too big for the mesh.

- Q: Which radius should I use? P? n_rad? 
A: A quick read to the paper mentioned above can help you in this. About the radius, try to set it such that a sphere of that radius cointains the feature that is repeated in the pattern. This is not a general rule, but practically speaking it is the most efficent one. About the other parameter, the most efficient settings we found across multiple datasets are P=15, nrad=5. 

- Q: I'm looking your code, the structure at line xx is odd/doesn't make sense. Are you hiding something?
A: This is because the code I'm giving you is a simplified/cutted version of the advanced code I'm currently working on. Right now, I have no time for developing a new and simplier version from scrathes. I'm "rushing" this release because I had many requests for the edgeLBP code and I do not want to make them wait that long.

- Q: I'm using your code on this dataset you have used on your paper but my results are different!
A: Of course, other factors must be considered, like how the descriptor is computed, which implementation of the distance measure between descriptor is used, and so on. If you need help in this, contact me. I'll answer when I can. 
