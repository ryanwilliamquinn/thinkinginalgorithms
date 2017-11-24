import random
import time
import closestpair
from matplotlib import pyplot
import alg_project3_viz
from alg_cluster import Cluster


DIRECTORY = "http://commondatastorage.googleapis.com/codeskulptor-assets/"
DATA_3108_URL = DIRECTORY + "data_clustering/unifiedCancerData_3108.csv"
DATA_896_URL = DIRECTORY + "data_clustering/unifiedCancerData_896.csv"
DATA_290_URL = DIRECTORY + "data_clustering/unifiedCancerData_290.csv"
DATA_111_URL = DIRECTORY + "data_clustering/unifiedCancerData_111.csv"


def gen_random_clusters(num_clusters):
    """
    creates a list of clusters where each cluster in this list corresponds to
    one randomly generated point in the square with corners (plus/minus1,plus/minus1)

    """
    return [ Cluster([], random.uniform(-1, 1), random.uniform(-1, 1), 1, 0) for _ in range(num_clusters)]

def compute_closest_pair_runtimes():
    """
    compute the running times of the functions slow_closest_pair and fast_closest_pair
    for lists of clusters of size 2 to 200
    """
    slow_closest_pair_runtimes = []
    cluster_sizes = range(2, 200)
    for num_clusters in cluster_sizes:
        start = time.time()
        closestpair.slow_closest_pair(gen_random_clusters(num_clusters))
        slow_closest_pair_runtimes.append(time.time() - start)
    print "slow closest pair runtimes: " + str(slow_closest_pair_runtimes)

    fast_closest_pair_runtimes = []
    for num_clusters in cluster_sizes:
        start = time.time()
        closestpair.fast_closest_pair(gen_random_clusters(num_clusters))
        fast_closest_pair_runtimes.append(time.time() - start)

    print "fast closest pair runtimes: " + str(fast_closest_pair_runtimes)

    pyplot.plot(cluster_sizes, fast_closest_pair_runtimes, '-b', label='fast closest pair runtimes')
    pyplot.plot(cluster_sizes, slow_closest_pair_runtimes, '-r', label='slow closest pair runtimes')
    pyplot.legend(loc='upper right')
    pyplot.ylabel('runtime in seconds')
    pyplot.xlabel('number of clusters in graph')
    pyplot.title('slow closest pair vs fast closest pair running time')
    pyplot.show()


def compute_distortion(cluster_list, data_table):
    """
    takes a list of clusters and uses cluster_error to compute its distortion
    """
    return sum(map(lambda cluster: cluster.cluster_error(data_table), cluster_list))


def test_distortion():
    """
    sum up the errors for all clusters

    """
    data_urls = [DATA_111_URL, DATA_290_URL, DATA_896_URL]
    for idx, url in enumerate(data_urls):
        if idx == 0:
            title = '111'
        elif idx == 1:
            title = '290'
        else:
            title = '896'

        data_table = alg_project3_viz.load_data_table(url)

        singleton_list = []
        for line in data_table:
            singleton_list.append(Cluster(set([line[0]]), line[1], line[2], line[3], line[4]))

        ascending_cluster_range = range(6, 21)

        cluster_lists = map(lambda cluster_count: closestpair.kmeans_clustering(singleton_list, cluster_count, 5),
                            ascending_cluster_range)

        kmeans_cluster_distortion = map(lambda cluster_list: compute_distortion(cluster_list, data_table),
                                        cluster_lists)

        singleton_list = []
        for line in data_table:
            singleton_list.append(Cluster(set([line[0]]), line[1], line[2], line[3], line[4]))

        cluster_range = range(20, 5, -1)
        hierarchical_cluster_distortion = []

        previous_clusters = singleton_list
        for num_clusters in cluster_range:
            new_clusters = closestpair.hierarchical_clustering(previous_clusters, num_clusters)
            hierarchical_cluster_distortion.append(compute_distortion(new_clusters, data_table))
            previous_clusters = new_clusters

        hierarchical_cluster_distortion.reverse()

        pyplot.plot(ascending_cluster_range, hierarchical_cluster_distortion, '-b', label='hierarchical clustering')
        pyplot.plot(ascending_cluster_range, kmeans_cluster_distortion, '-r', label='k-means clustering')
        pyplot.legend(loc='upper right')
        pyplot.ylabel('cluster distortion')
        pyplot.xlabel('number of clusters in graph')
        pyplot.title('distortion from hierarchical and k-means clustering: ' + title + ' counties')
        pyplot.show()


test_distortion()

