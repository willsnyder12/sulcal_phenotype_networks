##Sulcal phenotype measurement functions that require brainvisa's python distribution
##Note that trimesh must be installed by entering into BrainVISA's bash distribution
##(e.g., `bv bash`) and using `pip install trimesh`

#load in relevant modules
from soma import aims, aimsalgo
import numpy as np
from scipy.spatial import distance
import scipy.sparse.csgraph as csg
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
from skimage.segmentation import find_boundaries
from skimage.morphology import binary_dilation
from scipy import ndimage
import trimesh
import skimage
from scipy.sparse.csgraph import connected_components
from scipy.spatial import distance, ConvexHull, convex_hull_plot_2d
from smallestenclosingcircle import *
from scipy.signal import lombscargle
from scipy import interpolate
from scipy import optimize
from FractalDimension import fractal_dimension

#relabel a sulcal graph according to how you specify the merging of brainvisa sulcal labels
#N.B., this will overwrite original sulcal labels, so be cautious if you want to preserve those
#One good approach is to copy all data to a temporary directory when computing sulcal phenotypes
#so that you can save original brainvisa labeling if you want to merge sulcal labels differently
#in the future
def relabel_graph(graph, relabel_coding, hemi,graph_path):
    for v in graph.vertices():
        if 'label' in v:
            for sulcus in relabel_coding:
                if v['label'].replace("._" + hemi, "") in sulcus.split(","):
                    v['label'] = sulcus.split(",")[0] + "._" + hemi
    aims.write(graph, graph_path,
                   options={"save_only_modified": False})

#write out a 3d array in which sulcal exteriors (hull junctions) are labeled by 
#their index in the list of default BrainVISA sulcal labels
def write_hull_junctions(graph,skel,sulcus_list, hemi):
    skel_array = np.asarray(skel)
    hj_image = np.squeeze(np.zeros(skel_array.shape))
    hj_image.fill(-1)
    #loop through sulcal branches, denoted as vertices in BrainVISA sulcal graphs
    for v in graph.vertices():
        if 'label' in v:
            if v['label'] in sulcus_list:
                label_id = sulcus_list.index(v['label'])
                for e in v.edges():
                    if e.getSyntax() == 'hull_junction':
                        if 'aims_junction' in e:
                            hj = e['aims_junction']
                            for i in xrange(hj[0].size()):
                                    coord = hj[0].keys()[i]
                                    hj_image[coord[0], coord[1], coord[2]] = label_id
    return hj_image

#A recreation of the BrainVISA function used to get depth profiles,
#which we will use these instead to get median / depth variability
#This function computes the geodesic distance at each point along the sulcal medial wall
#to the nearest hull junction / sulcal exterior point. Therefore, hull junction points
#will have zero depth, and accessing sulcal fundal points will get geodesic depth from 
#which you can get a depth profile
def compute_depthmap(graph, skel):
    fat = aims.FoldGraphAttributes(skel, graph)
    fat.prepareDepthMap()
    depth = fat.getDepth()
    depth_factor = fat.getDepthfactor()
    converter = aims.Converter_Volume_S16_Volume_FLOAT()
    depth2 = aims.Volume_FLOAT(depth.dimX(), depth.dimY(), depth.dimZ(), depth.dimT())
    converter.convert(depth, depth2)
    depth = depth2
    depth[:] /= depth_factor
    return depth

#This function loops through sulcal fundal points to collect the geodesic depth profile for each sulcus
def compute_depth_profiles(graph,skel,sulcus_list):
    geodesic_depth_profile = {new_list: np.array([]) for new_list in range(len(sulcus_list))}
    depth = compute_depthmap(graph,skel)
    depth_array = depth.arraydata()
    #loop through sulcal branches
    for v in graph.vertices():
        if 'label' in v:
            #only use labeled branches
            if v['label'] in sulcus_list:
                #use sulcal fundal or "bottom" points
                if 'aims_Tmtktri' in v and 'aims_bottom' in v:
                    sulcus_id = sulcus_list.index(v['label'])
                    bottom_map = v['aims_bottom']
                    depth_vals_geo = []
                    #loop through each bottom point / sulcal fundal point and get depth map value
                    for i in xrange(bottom_map[0].size()):
                        coord = bottom_map[0].keys()[i]
                        p2 = (0, coord[2], coord[1], coord[0])
                        sulcus_depth_geo = depth_array[p2]
                        try:
                            depth_vals_geo.append(sulcus_depth_geo)
                        except:
                            pass
                    geodesic_depth_profile[sulcus_id] = np.append(geodesic_depth_profile[sulcus_id], depth_vals_geo)
    return geodesic_depth_profile

##A helper function for compute_largest_cc_length. This function takes in an array defining
#which branches are connected to which other sulcal branches (by keeping track of their corresponding indices),
#converts this into an adjacency matrix, uses graph theory to extract connected components, and then writes out
#the branch lengths in terms of contiguous sulcal branches
def collapse_cc_lengths(sulcus_indices, hj_lengths):
    #hj_lengths must be key value pair, with sulcus ID as key and length as value
    sulcus_indices = np.asarray(sulcus_indices)
    max_idx = int(sulcus_indices.max())
    adj = np.zeros((max_idx, max_idx))
    for row in sulcus_indices:
        adj[int(row[0]-1), int(row[1]-1)] = 1
        adj[int(row[1]-1), int(row[0]-1)] = 1
    connectivity = np.asarray(csg.connected_components(adj)[1])
    u, c = np.unique(connectivity, return_counts=True)
    repeated_vals = np.asarray(u[c > 1])
    collapsed_branch_lengths = []
    for x in repeated_vals:
        tmp_component = np.where(connectivity == x)
        branch_length = 0
        for y in tmp_component[0]:
            z = y + 1
            branch_length = branch_length + hj_lengths[z]
        collapsed_branch_lengths.append(branch_length)
    return(collapsed_branch_lengths)

#This function collects all of the branches with a given sulcal label,
#identifies which of these branches are connected using an adjacency matrix of their
#connectivity (connection = physically touching, distance == 0), then reformats all of the
#lengths of branches in terms of their connected components so that the longest connected
#sulcal branch can be extracted
def compute_largest_cc_length(graph,sulcus_list):
    cc_id_list = {new_list: np.empty((0,2), float) for new_list in range(len(sulcus_list))}
    hull_junction_length_list = {new_list: {} for new_list in range(len(sulcus_list))}
    branch_lengths = {new_list: np.array([]) for new_list in range(len(sulcus_list))}
    #loop through all possible branches
    for v in graph.vertices():
        if 'label' in v:
            if v['label'] in sulcus_list:
                sulcus_id = sulcus_list.index(v['label'])
                pp_flag = 0
                #loop through all edges for that branch (edge = relation with other sulci)
                for e in v.edges():
                    if e.getSyntax() == 'cortical':
                        #if two branches have zero distance between them, e.g. are connected
                        if e.vertices()[0]['label'] == e.vertices()[1]['label'] and e['dist'] == 0:
                            try:
                                del v1_length
                            except:
                                pass
                            try:
                                del v2_length
                            except:
                                pass
                            #get lengths of each branch
                            v1 = e.vertices()[0]
                            v2 = e.vertices()[1]       
                            for e2 in v1.edges():
                                if e2.getSyntax() == 'hull_junction':
                                    v1_length = e2['length']
                            for e3 in v2.edges():
                                if e3.getSyntax() == 'hull_junction':
                                    v2_length = e3['length']
                            #collect lengths and indices for collapse of contiguous branches in collapse_cc_lengths
                            if 'v1_length' in locals() and 'v2_length' in locals():
                                pp_flag = 1
                                cc_id_list[sulcus_id] = np.vstack((cc_id_list[sulcus_id], [v1['index'], v2['index']]))
                                hull_junction_length_list[sulcus_id].update({v1['index']: v1_length, v2['index']: v2_length})
                #if a sulcal branch is isolated, i.e. no plis de passage (pp) or junction edges with other sulci, its length is 
                #already in the format needed for final computation of max branch length
                if pp_flag == 0:
                    for e in v.edges():
                        if e.getSyntax() == 'hull_junction':
                            hj_length = e['length']
                            branch_lengths[sulcus_id] = np.append(branch_lengths[sulcus_id], hj_length)
    #collapse all branch lengths to get contiguous sulci lengths, from which the max is later computed in combine_hemi_features
    for x in hull_junction_length_list:
        if len(hull_junction_length_list[x]) > 0:
            new_lengths = collapse_cc_lengths(cc_id_list[x], hull_junction_length_list[x])
            branch_lengths[x] = np.append(branch_lengths[x], new_lengths)
    return branch_lengths

#Branch span is computed by projecting hull junction points (sulcus exterior) from 3d space into a 2d plane 
#defined by the point in the center of mass of the longest sulcal branch and the direction tangent to the smoothed
#hull of the cortex at that point,
#and then subsequently finding the circumscribed area about these projected points and the convex hull area of these points
#to get a ratio of how "circular" the extent of the points are
def compute_branch_span(hj_image,skel,sulcus_list):
    branch_span_values = {new_list: np.array([]) for new_list in range(len(sulcus_list))}
    skel_array = np.asarray(skel)
    threshed = skel_array != 11
    binary_hull = threshed.astype(int)
    verts, faces = skimage.measure.marching_cubes_classic(np.squeeze(binary_hull),0)
    tri_mesh = trimesh.Trimesh(vertices = verts, faces = faces)
    #use trimesh to smooth the hull further (hypothetical convex hull / smoothed or unfolded representation of the cortex)
    #to provide good estimates of the normal at each point in the hull
    smoothed_mesh = trimesh.smoothing.filter_laplacian(tri_mesh, lamb=1, iterations=100, implicit_time_integration=False, volume_constraint=True, laplacian_operator=None)
    for i in range(0, len(sulcus_list)):
        if np.sum(hj_image == i) > 0:
            #collect all sulcal branches for a given sulcal label as the connected components of the sulcus' hull junction
            sulcus_coords = np.where(hj_image == i)
            sulcus_coords_t = np.transpose(sulcus_coords)
            overall_centroid = np.mean(sulcus_coords_t, axis = 0)
            #structure defines connectivity
            coord_graph = distance.cdist(sulcus_coords_t,sulcus_coords_t)
            coord_adj = (coord_graph < 2).astype('int')
            ccs = connected_components(coord_adj)
            n_cc = ccs[0]
            
            longest_cc = 0
            longest_cc_length = 0
            centroids = np.empty((0,3), float)
            #find the longest sulcal branch
            for j in range(0,n_cc):
                n_voxels = np.sum(ccs[1] == j)
                #branch_coords_j = np.where(labels_out == j + 1)
                #branch_coords_j_t = np.transpose(branch_coords_j)
                cc_idx = np.where(ccs[1] == j)[0]
                branch_coords_j_t = sulcus_coords_t[cc_idx,:]
                centroid_j =  np.mean(branch_coords_j_t, axis = 0)
                centroids = np.vstack((centroids, centroid_j))

                if n_voxels > longest_cc_length:
                    longest_cc_length = n_voxels
                    longest_cc = j
            #remap all points to the center of mass of the longest sulcal branch
            new_coordinates = np.empty((0,3), float)
            for j in range(0,n_cc):
                shift = centroids[longest_cc] - centroids[j]
                #add shift
                cc_idx = np.where(ccs[1] == j)[0]
                branch_coords_j_t = sulcus_coords_t[cc_idx,:]

                shifted_sulc = branch_coords_j_t + shift
                new_coordinates = np.vstack((new_coordinates, shifted_sulc))
            #get the point on smoothed cortex closest to the center of mass of the longest sulcal branch
            #and also find the normal to this point to define the tangent plane
            dists = distance.cdist([overall_centroid], smoothed_mesh.vertices)
            closest_pts = np.argsort(dists[0])[:10]
            smooth_norm_sample = smoothed_mesh.vertex_normals[closest_pts]
            smooth_norm = np.mean(smooth_norm_sample, axis = 0)
            ref_pt = smoothed_mesh.vertices[closest_pts[0]]

            #remap all points into a 2d basis, anchored at the point on the smoothed mesh and along the tangent plane
            basis_1 = np.array([-1 * smooth_norm[1], smooth_norm[0], 0])
            basis_1_norm = basis_1 / np.sqrt(np.sum(basis_1**2))
            basis_2 = np.cross(basis_1,smooth_norm)
            basis_2_norm = basis_2 / np.sqrt(np.sum(basis_2**2))

            new_x = np.dot(new_coordinates - ref_pt, basis_1 )
            new_y = np.dot(new_coordinates - ref_pt, basis_2 )
            
            two_d_point_cloud = np.transpose(np.vstack((new_x,new_y)))
            #get circumscribed area and convex hull area of this 2d point set
            if two_d_point_cloud.shape[0] > 3:
                hull = ConvexHull(two_d_point_cloud)

                enscribed_circle = make_circle(two_d_point_cloud)
                r = enscribed_circle[2]
                theta = np.linspace(0, 2*np.pi, 100)
                x1 = r*np.cos(theta) + enscribed_circle[0]
                x2 = r*np.sin(theta) + enscribed_circle[1]

                #branch span = convex hull area (called "volume" in module's documentation) / circle area
                branching_index = hull.volume / (np.pi * r**2)
                branch_span_values[i] = np.append(branch_span_values[i], branching_index)

    return(branch_span_values)

#compute fractal dimension from the patterning of hull junction points using a 3d fractal dimension function
def compute_fractal_dimension(hj_image,sulcus_list):
    fd_values = {new_list: np.array([]) for new_list in range(len(sulcus_list))}
    for i in range(0, len(sulcus_list)):
        if np.sum(hj_image == i) > 0:
            sulc_image = (hj_image == i).astype('int')
            fd = fractal_dimension(sulc_image)
            fd_values[i] = np.append(fd_values[i], fd)

    return fd_values

#combine all features measured on each hemisphere so everything can be put into one, cleaned up csv file
def combine_hemi_features(sulcus_list,
                        geodesic_depth_left,
                        geodesic_depth_right,
                        left_longest_branch,
                        right_longest_branch,
                        branch_span_left,
                        branch_span_right,
                        fd_left,
                        fd_right):
    #initialize each with NaN
    average_depth = np.empty(len(sulcus_list))
    average_depth[:] = np.NaN
    depth_variability = np.empty(len(sulcus_list))
    depth_variability[:] = np.NaN
    longest_branch = np.empty(len(sulcus_list))
    longest_branch[:] = np.NaN
    branch_span = np.empty(len(sulcus_list))
    branch_span[:] = np.NaN
    fd = np.empty(len(sulcus_list))
    fd[:] = np.NaN
    for i in range(0, len(sulcus_list)):
        if np.any(geodesic_depth_left[i]):
            #average depth is given by median
            average_depth[i] = np.median(geodesic_depth_left[i])
            #depth variability is given by median absolute deviation
            depth_variability[i] = np.median(np.abs(geodesic_depth_left[i] - np.median(geodesic_depth_left[i])))
            if np.any(left_longest_branch[i]):
                longest_branch[i] = left_longest_branch[i].max()
            if np.any(branch_span_left[i]):
                branch_span[i] = branch_span_left[i]
            if np.any(fd_left[i]):
                fd[i] = fd_left[i]
        elif np.any(geodesic_depth_right[i]):
            average_depth[i] = np.median(geodesic_depth_right[i])
            depth_variability[i] = np.median(np.abs(geodesic_depth_right[i] - np.median(geodesic_depth_right[i])))
            if np.any(right_longest_branch[i]):
                longest_branch[i] = right_longest_branch[i].max()
            if np.any(branch_span_right[i]):
                branch_span[i] = branch_span_right[i]
            if np.any(fd_right[i]):
                fd[i] = fd_right[i]
    return average_depth, depth_variability, longest_branch, branch_span, fd
