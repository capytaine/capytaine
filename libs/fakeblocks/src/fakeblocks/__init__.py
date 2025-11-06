    #
    # def sliced_by_plane(self, plane):
    #     """Return the same body, but replace the mesh by a set of two meshes
    #     corresponding to each sides of the plane."""
    #     return FloatingBody(mesh=self.mesh.sliced_by_plane(plane),
    #                         lid_mesh=self.lid_mesh.sliced_by_plane(plane)
    #                                     if self.lid_mesh is not None else None,
    #                         dofs=self.dofs, name=self.name)
    #
    # def minced(self, nb_slices=(8, 8, 4)):
    #     """Experimental method decomposing the mesh as a hierarchical structure.
    #
    #     Parameters
    #     ----------
    #     nb_slices: Tuple[int, int, int]
    #         The number of slices in each of the x, y and z directions.
    #         Only powers of 2 are supported at the moment.
    #
    #     Returns
    #     -------
    #     FloatingBody
    #     """
    #     minced_body = self.copy()
    #
    #     # Extreme points of the mesh in each directions.
    #     x_min, x_max, y_min, y_max, z_min, z_max = self.mesh.axis_aligned_bbox
    #     sizes = [(x_min, x_max), (y_min, y_max), (z_min, z_max)]
    #
    #     directions = [np.array(d) for d in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]]
    #
    #     def _slice_positions_at_depth(i):
    #         """Helper function.
    #
    #         Returns a list of floats as follows:
    #         i=1 -> [1/2]
    #         i=2 -> [1/4, 3/4]
    #         i=3 -> [1/8, 3/8, 5/8, 7/8]
    #                ...
    #         """
    #         denominator = 2**i
    #         return [numerator/denominator for numerator in range(1, denominator, 2)]
    #
    #     # GENERATE ALL THE PLANES THAT WILL BE USED TO MINCE THE MESH
    #     planes = []
    #     for direction, nb_slices_in_dir, (min_coord, max_coord) in zip(directions, nb_slices, sizes):
    #         planes_in_dir = []
    #
    #         depth_of_treelike_structure = int(np.log2(nb_slices_in_dir))
    #         for i_depth in range(1, depth_of_treelike_structure+1):
    #             planes_in_dir_at_depth = []
    #             for relative_position in _slice_positions_at_depth(i_depth):
    #                 slice_position = (min_coord + relative_position*(max_coord-min_coord))*direction
    #                 plane = Plane(normal=direction, point=slice_position)
    #                 planes_in_dir_at_depth.append(plane)
    #             planes_in_dir.append(planes_in_dir_at_depth)
    #         planes.append(planes_in_dir)
    #
    #     # SLICE THE MESH
    #     intermingled_x_y_z = chain.from_iterable(zip_longest(*planes))
    #     for planes in intermingled_x_y_z:
    #         if planes is not None:
    #             for plane in planes:
    #                 minced_body = minced_body.sliced_by_plane(plane)
    #     return minced_body
    #
    # def cluster_bodies(*bodies, name=None):
    #     """
    #     Builds a hierarchical clustering from a group of bodies
    #
    #     Parameters
    #     ----------
    #     bodies: list
    #         a list of bodies
    #     name: str, optional
    #         a name for the new body
    #
    #     Returns
    #     -------
    #     FloatingBody
    #         Array built from the provided bodies
    #     """
    #     from scipy.cluster.hierarchy import linkage
    #     nb_buoys = len(bodies)
    #
    #     if any(body.center_of_buoyancy is None for body in bodies):
    #         raise ValueError("The center of buoyancy of each body needs to be known for clustering")
    #     buoys_positions = np.stack([body.center_of_buoyancy for body in bodies])[:,:2]
    #
    #     ln_matrix = linkage(buoys_positions, method='centroid', metric='euclidean')
    #
    #     node_list = list(bodies)  # list of nodes of the tree: the first nodes are single bodies
    #
    #     # Join the bodies, with an ordering consistent with the dendrogram.
    #     # Done by reading the linkage matrix: its i-th row contains the labels
    #     # of the two nodes that are merged to form the (n + i)-th node
    #     for ii in range(len(ln_matrix)):
    #         node_tag = ii + nb_buoys # the first nb_buoys tags are already taken
    #         merge_left = int(ln_matrix[ii,0])
    #         merge_right = int(ln_matrix[ii,1])
    #         # The new node is the parent of merge_left and merge_right
    #         new_node_ls = [node_list[merge_left], node_list[merge_right]]
    #         new_node = FloatingBody.join_bodies(*new_node_ls, name='node_{:d}'.format(node_tag))
    #         node_list.append(new_node)
    #
    #     # The last node is the parent of all others
    #     all_buoys = new_node
    #
    #     if name is not None:
    #         all_buoys.name = name
    #
    #     return all_buoys
