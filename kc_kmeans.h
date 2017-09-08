//
// Copyright 2017 Wyatt Miller
//

    #ifndef KC_KMEANS_H
    #define KC_KMEANS_H

    #include <fstream>
    #include <list>
    #include <string>
    #include <vector>

namespace kc {

    struct cluster_wrapper {
        std::vector<double> centroid;
        unsigned int count;
        bool move_flag;   // Signifies if anything moved in or out of a cluster - avoids unneeded centroid recomputation
        cluster_wrapper() : count(0), move_flag(1) {}
    };


    struct point_wrapper {
        std::vector<double> point;
        cluster_wrapper *cluster_pointer;

        point_wrapper() : cluster_pointer(nullptr) {}

    };


    class kmeans_set {
        std::list<point_wrapper> universe;     // The universe of points, a vector of vectors (vector of points)
        std::list<cluster_wrapper> clusters;   // Clusters with pointers to each point

        double epsilon;            // The maximum level of error (default 0)
        unsigned int hard_limit;   // Hard limit on the number of iterations (defaults to max of unsigned int)

        // I/O
        void import_points(std::list<std::vector<double>> &);          // Used by both import functions
        void print_point(std::ostream &, std::vector<double>, char);   // A quick way to observe output

        // Primary Algo Functions
        void move_point_to_nearest_cluster();                                      // "Moves" points to the nearest cluster/centroid via pointer manipulation
        double calculate_distance(std::vector<double> &, std::vector<double> &);   // 2-norm distance calculation (not rooted to save computation)
        double recompute_centroids();                                              // Iterates through all clusters and recomputes centroids
        bool check_move_state();                                                   // Checks if points have moved in/out of a cluster - avoids centroid recomputing when no movement

        // Reporting (commented out in .cpp)
        void report_point_move(std::ostream &, std::vector<double>, cluster_wrapper &, cluster_wrapper &);

        void report_universe_state(std::ostream &, double, double);

    public:
        kmeans_set(std::ifstream &, char);              // Import points from file with specified delimiter
        kmeans_set(std::list<std::vector<double>> &);   // Import points from vector on construction

        void compute_centroids(unsigned int);           // Primary Algorithm
        void print_centroids(std::ostream &, char);     // Output stream, centroid delimiter
        void burst_clusters(std::string, char);         // Centroid output stream, starting string for each cluster file, point delimiter

        void set_epsilon(double eps) { this->epsilon = eps; }

        void set_hard_limit(unsigned int lim) { this->hard_limit = lim; }
    };

}

#endif //KC_KMEANS_H