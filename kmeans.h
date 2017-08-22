//
// Copyright 2017 Wyatt Miller
//

    #ifndef KMEANS_KMEANS_H
    #define KMEANS_KMEANS_H

    #include <fstream>
    #include <list>
    #include <vector>



    struct cluster_wrapper {
        std::vector<double> centroid;
        unsigned int count;
        bool move_flag;   // Signifies if anything moved in or out of a cluster - avoids unneeded centroid recomputation
        cluster_wrapper() : count(0), move_flag(1) {}
    };



    struct point_wrapper {
        std::vector<double> point;
        cluster_wrapper * centroid_pointer;
        point_wrapper() : centroid_pointer(nullptr) {}
    };



    class kmeans_set {
        std::list<point_wrapper> universe;      // The universe of points, a vector of vectors (vector of points)
        std::list<cluster_wrapper> clusters;   // Clusters with pointers to each point

        double epsilon;            // The maximum level of error (default 0)
        unsigned int hard_limit;   // Hard limit on the number of iterations (defaults to max of unsigned int)

        void import_points(std::list<std::vector<double>> &);   // Used by both import functions

        double calculate_distance(std::vector<double> &, std::vector<double> &);   // 2-norm distance calculation
        void move_to_nearest();                                                    // "Moves" points to the nearest cluster/centroid via pointer manipulation
        double recompute_centroids();
        bool check_move_state();                                                   // Checks if points have moved in/out of a cluster - avoids centroid recomputing when no movement
        void print_point(std::ostream &, std::vector<double>, char);               // A quick way to observe output

    public:
        kmeans_set(std::ifstream &, char);                // Import points from file with specified delimiter
        kmeans_set(std::list<std::vector<double>> &);   // Import points from vector on construction

        void compute_centroids(unsigned int);
        void print_centroids(std::ostream &, char);       // output stream, centroid delimiter
        void print_clusters(std::ostream &, char, char);  // output stream, point delimiter, centroid designator

        void set_epsilon(double eps) { this->epsilon = eps; }
        void set_hard_limit(unsigned int lim) { this->hard_limit = lim; }
    };

    #endif //KMEANS_KMEANS_H
