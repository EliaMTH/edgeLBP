STRUCTURE of the "data" output.
------------------------------------------------------------------------
data.prel 	-> a list of most of the variables used multiple time in the 
				edgeLBP function: 
				P: 		number of samples taken on each e-ring;
				R: 		size of the maximum radius;
				v: 		vertex list;
				f: 		face list;
				nrad: 	number of rings considered per vertex;
				vv:		vertex-vertex relation;
				fe:		face-edge relation;
				e:		edge list;
				ve:		vertex-edge relation;
				bv:		list of the vertex that are on the boundary of the
						mesh;
				iv:		list of the vertex NOT on the boundary;
				vld:	list of admissible vertices;
				nv:		number of vertices;
				ni:		size of iv;
				nb:		size of nb;
				nf:		number of faces;
				ne:		number of edges;
				ef:		edge-face relation;
				ee:		edge-edge relation;
				h:		pattern descriptor, sampled on the vertices;
				
data.tr		-> transitions between the 0 and 1 in the array of each ring,
				the structure is the same as that of data.a1; 
				
data.a1		-> the same as the a1 output of the function edgeLBPp.
------------------------------------------------------------------------
