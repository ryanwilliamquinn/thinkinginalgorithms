"""
Student template code for Project 3
Student will implement five functions:

slow_closest_pair(cluster_list)
fast_closest_pair(cluster_list)
closest_pair_strip(cluster_list, horiz_center, half_width)
hierarchical_clustering(cluster_list, num_clusters)
kmeans_clustering(cluster_list, num_clusters, num_iterations)

where cluster_list is a 2D list of clusters in the plane
"""

import math
import alg_cluster



######################################################
# Code for closest pairs of clusters

def pair_distance(cluster_list, idx1, idx2):
    """
    Helper function that computes Euclidean distance between two clusters in a list

    Input: cluster_list is list of clusters, idx1 and idx2 are integer indices for two clusters

    Output: tuple (dist, idx1, idx2) where dist is distance between
    cluster_list[idx1] and cluster_list[idx2]
    """
    return cluster_list[idx1].distance(cluster_list[idx2]), min(idx1, idx2), max(idx1, idx2)


def slow_closest_pair(cluster_list):
    """
    Compute the distance between the closest pair of clusters in a list (slow)

    Input: cluster_list is the list of clusters

    Output: tuple of the form (dist, idx1, idx2) where the centers of the clusters
    cluster_list[idx1] and cluster_list[idx2] have minimum distance dist.
    """
    closest_points = (float("inf"), -1, -1)
    for idxa, _ in enumerate(cluster_list):
        for idxb, _ in enumerate(cluster_list):
            if (idxa == idxb):
                continue

            cluster_distance = pair_distance(cluster_list, idxa, idxb)
            if cluster_distance[0] < closest_points[0]:
                closest_points = (cluster_distance[0], idxa, idxb)

    return closest_points


def fast_closest_pair(cluster_list):
    """
    Compute the distance between the closest pair of clusters in a list (fast)

    Input: cluster_list is list of clusters SORTED such that horizontal positions of their
    centers are in ascending order

    Output: tuple of the form (dist, idx1, idx2) where the centers of the clusters
    cluster_list[idx1] and cluster_list[idx2] have minimum distance dist.
    """
    cluster_list_length = len(cluster_list)
    if cluster_list_length <= 3:
        return slow_closest_pair(cluster_list)
    else:
        cluster_midpoint = int(cluster_list_length / 2)
        left_half_cluster_list = cluster_list[0:cluster_midpoint]
        right_half_cluster_list = cluster_list[cluster_midpoint:]
        left_closest_pair = fast_closest_pair(left_half_cluster_list)
        right_closest_pair = fast_closest_pair(right_half_cluster_list)
        closest_pair = get_min_closest_pair(left_closest_pair, right_closest_pair, cluster_midpoint)
        midpoint_x_coord = 0.5 * (cluster_list[cluster_midpoint - 1].horiz_center() +
                                  cluster_list[cluster_midpoint].horiz_center())
        return min(closest_pair, closest_pair_strip(cluster_list, midpoint_x_coord, closest_pair[0]))


def get_min_closest_pair(left_closest_pair, right_closest_pair, midpoint):
    """
    Return a tuple with the distance and indices of the closest pair for two different pairs
    """
    if left_closest_pair[0] < right_closest_pair[0]:
        return left_closest_pair
    else:
        return right_closest_pair[0], right_closest_pair[1] + midpoint, right_closest_pair[2] + midpoint


def check_middle_width(xcoord, horiz_center, half_width):
    """
    Check if the xcoord of a point is within half_width of the horiz_center xcoord
    """
    return math.fabs(xcoord - horiz_center) < half_width



def closest_pair_strip(cluster_list, horiz_center, half_width):
    """
    Helper function to compute the closest pair of clusters in a vertical strip

    Input: cluster_list is a list of clusters produced by fast_closest_pair
    horiz_center is the horizontal position of the strip's vertical center line
    half_width is the half the width of the strip (i.e; the maximum horizontal distance
    that a cluster can lie from the center line)

    Output: tuple of the form (dist, idx1, idx2) where the centers of the clusters
    cluster_list[idx1] and cluster_list[idx2] lie in the strip and have minimum distance dist.
    """

    filtered_cluster_list = filter(lambda cluster: check_middle_width(cluster.horiz_center(), horiz_center, half_width), cluster_list)
    filtered_cluster_list_length = len(filtered_cluster_list)
    filtered_cluster_list.sort(key = lambda cluster: cluster.vert_center())
    closest_points = (float("inf"), -1, -1)
    for idxa in range(0, filtered_cluster_list_length - 1):
        for idxb in range(idxa+1, min(idxa+4, filtered_cluster_list_length)):
            cluster_distance = pair_distance(filtered_cluster_list, idxa, idxb)
            if cluster_distance[0] <= closest_points[0]:
                closest_points = (cluster_distance[0], idxa, idxb)
    # find the indices of idxa and idxb in the original list
    if closest_points[0] == float("inf"):
        return closest_points
    closest_indices = (cluster_list.index(filtered_cluster_list[closest_points[1]]), cluster_list.index(filtered_cluster_list[closest_points[2]]))
    return (closest_points[0], min(closest_indices), max(closest_indices))



######################################################################
# Code for hierarchical clustering


def hierarchical_clustering(cluster_list, num_clusters):
    """
    Compute a hierarchical clustering of a set of clusters
    Note: the function may mutate cluster_list

    Input: List of clusters, integer number of clusters
    Output: List of clusters whose length is num_clusters
    """
    while len(cluster_list) > num_clusters:
        # sort the list
        cluster_list.sort(key=lambda cluster: cluster.horiz_center())
        closest_cluster_indices = fast_closest_pair(cluster_list)
        closest_cluster_a = cluster_list[closest_cluster_indices[1]]
        closest_cluster_b = cluster_list[closest_cluster_indices[2]]
        closest_cluster_a.merge_clusters(closest_cluster_b)
        del cluster_list[closest_cluster_indices[2]]
    return cluster_list


######################################################################
# Code for k-means clustering


def kmeans_clustering(cluster_list, num_clusters, num_iterations):
    """
    Compute the k-means clustering of a set of clusters
    Note: the function may not mutate cluster_list

    Input: List of clusters, integers number of clusters and number of iterations
    Output: List of clusters whose length is num_clusters
    """
    cluster_list_copy = map(lambda cluster: cluster.copy(), cluster_list)

    # position initial clusters at the location of clusters with largest populations
    sorted_cluster_list = sorted(cluster_list_copy, key=lambda cluster: cluster.total_population(), reverse=True)
    clusters = sorted_cluster_list[0:num_clusters]
    for _ in range(num_iterations):
        cluster_copy = map(lambda cluster: cluster.copy(), cluster_list)
        new_clusters = [None for _ in range(num_clusters)]
        min_pair_distance = float("inf")
        closest_cluster_idx = None
        for county in cluster_copy:
            for idx, old_cluster in enumerate(clusters):
                distance = county.distance(old_cluster)
                if distance < min_pair_distance:
                    min_pair_distance = distance
                    closest_cluster_idx = idx
            if new_clusters[closest_cluster_idx] is not None:
                new_clusters[closest_cluster_idx].merge_clusters(county)
            else:
                new_clusters[closest_cluster_idx] = county
            min_pair_distance = float("inf")
            closest_cluster_idx = None
        clusters = list(new_clusters)
    return clusters
