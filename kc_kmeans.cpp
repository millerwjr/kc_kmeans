//
// Copyright 2017 Wyatt Miller
//

#include "kc_kmeans.h"
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>



// PUBLIC CONSTRUCTORS

// Import points from file with specified delimiter into a temporary list - calls import_points()
kc::kmeans_set::kmeans_set(std::ifstream & input_file, char delimiter) :
        epsilon(0),
        hard_limit(std::numeric_limits<unsigned int>::max()) {
    std::string line;
    std::list<std::vector<double>> temp_point_list;
    while (std::getline(input_file, line)) {
        while ((line.length() == 0) && !(input_file.eof())) {
            std::getline(input_file, line);   // Skips blank lines in file
        }
        std::string parameter;
        std::stringstream ss(line);
        std::vector<double> temp_point;
        if ((line.length() != 0)) {
            while (std::getline(ss, parameter, delimiter)) {
                temp_point.push_back(atof(parameter.c_str()));
            }
            temp_point_list.push_back(temp_point);
        }
    }
    this->import_points(temp_point_list);
}



// Import points from list on construction
kc::kmeans_set::kmeans_set(std::list<std::vector<double>> & point_list) :
        epsilon(0),
        hard_limit(std::numeric_limits<unsigned int>::max()) {
    this->import_points(point_list);
}



// PRIVATE FUNCTIONS

// Primary point import function - Assures dimensional integrity by comparing all intake to the first point in the [points] list
void kc::kmeans_set::import_points(std::list<std::vector<double>> & point_list) {
    for (auto point_list_iter = point_list.begin(); point_list_iter != point_list.end(); ++point_list_iter) {
        point_wrapper new_point;
        new_point.point = *point_list_iter;
        if (this->universe.empty()) {
            this->universe.push_back(new_point);
        } else if (this->universe.front().point.size() == point_list_iter->size()) {
            // Assures dimensional integrity
            // Sorting on intake assures well-dispursed centroid selection
            // Well-dispursed centroid selection helps mitigate centroids "fighting" over points
            auto universe_list_iter = this->universe.begin();
            while (universe_list_iter != this->universe.end() && (universe_list_iter->point < *point_list_iter)) {   // PRECEDENCE MATTERS!!
                // Sorts points on insert - this ensures evenly distributed clusters
                ++universe_list_iter;
            }
            this->universe.insert(universe_list_iter, new_point);
        }
    }
}



void kc::kmeans_set::print_point(std::ostream & output, std::vector<double> vec, char delim) {
    for (auto iter = vec.begin(); iter != vec.end(); ++iter) {
        if (iter != vec.begin()) { output << delim; }
        output << *iter;
    }
}



void kc::kmeans_set::move_point_to_nearest_cluster() {
    for (auto universe_iter = this->universe.begin(); universe_iter != this->universe.end(); ++universe_iter) {
        for (auto cluster_iter = this->clusters.begin(); cluster_iter != this->clusters.end(); ++cluster_iter) {
            if (calculate_distance(universe_iter->point, cluster_iter->centroid) < calculate_distance(universe_iter->point, universe_iter->cluster_pointer->centroid)) {
                // report_point_move(std::cout, universe_iter->point, *universe_iter->cluster_pointer, *&(*cluster_iter));   // REPORTING - comment out for silent run
                --universe_iter->cluster_pointer->count;             // Decrement the count in the old cluster
                universe_iter->cluster_pointer->move_flag = true;    // Flag the old cluster for centroid recalculation
                universe_iter->cluster_pointer = &(*cluster_iter);   // Move the point; This is why lists over vectors - iterator stays valid with a list
                ++universe_iter->cluster_pointer->count;             // Increment the count in the new cluster
                universe_iter->cluster_pointer->move_flag = true;    // Flag the new cluster for centroid recalculation
            }
        }
    }
}



double kc::kmeans_set::calculate_distance(std::vector<double> & a_vec, std::vector<double> & b_vec) {
    auto b_iter = b_vec.begin();
    double r_val = 0;
    for (auto a_iter = a_vec.begin(); a_iter != a_vec.end(); ++a_iter, ++b_iter) {
        r_val += (*a_iter - *b_iter) * (*a_iter - *b_iter);
    }
    return r_val;   // No need to root - actual distance not necessary
}



double kc::kmeans_set::recompute_centroids() {

    // The following is used for delta calculation
    std::list<std::vector<double>> old_centroids;
    for (auto cluster_iter = this->clusters.begin(); cluster_iter != this->clusters.end(); ++cluster_iter) {
        old_centroids.push_back(cluster_iter->centroid);   // Copy the current centroids state
    }
    for (auto cluster_iter = this->clusters.begin(); cluster_iter != this->clusters.end(); ++cluster_iter) {
        std::fill(cluster_iter->centroid.begin(), cluster_iter->centroid.end(), 0);   // Zero out the primary centroid list
    }
    auto old_centroid_iter = old_centroids.begin();
    for (auto cluster_iter = this->clusters.begin(); cluster_iter != this->clusters.end(); ++cluster_iter, ++old_centroid_iter) {
        if (!(cluster_iter->move_flag)) {
            cluster_iter->centroid = *old_centroid_iter;  // Recopy back to the current cluster if there was no movement - will avoid recomputation.
        }
    }

    // Recompute centroids
    for (auto universe_iter = this->universe.begin(); universe_iter != this->universe.end(); ++universe_iter) {
        if (universe_iter->cluster_pointer->move_flag) {
            // If the cluster has seen movement...
            auto cluster_pointer_centroid_iter = universe_iter->cluster_pointer->centroid.begin();
            for (auto point_iter = universe_iter->point.begin(); point_iter != universe_iter->point.end(); ++point_iter, ++cluster_pointer_centroid_iter) {
                // Division on the fly useful for accuracy super-large dimensions and prevents a second loop - doubles computation
                *cluster_pointer_centroid_iter += *point_iter / universe_iter->cluster_pointer->count;
            }
        }
    }

    // Compute max delta
    double max_delta = 0;
    auto old_centroids_iter = old_centroids.begin();
    for (auto centroids_iter = this->clusters.begin(); centroids_iter != this->clusters.end(); ++centroids_iter) {
        double test_delta = calculate_distance(*old_centroids_iter, centroids_iter->centroid);
        if (max_delta < test_delta) { max_delta = test_delta; }
    }
    return max_delta;
}



bool kc::kmeans_set::check_move_state() {
    for (auto cluster_iter = this->clusters.begin(); cluster_iter != this->clusters.end(); ++cluster_iter) {
        if (cluster_iter->move_flag) { return true; }
    }
    return false;
}



void kc::kmeans_set::report_point_move(std::ostream & output, std::vector<double> point, cluster_wrapper & fm, cluster_wrapper & to) {
    output << "\nMoving (";
    print_point(output, point, ',');
    output << ")   from cluster: " << &fm << "   to cluster: " << &to;
}



void kc::kmeans_set::report_universe_state(std::ostream & output, double d1, double d2) {
    output << "\n\n[Iteration Complete]";
    output << "\n    delta: " << fabs(d2 - d1);
    unsigned int ct = 0;
    for (auto cluster_iter = this->clusters.begin(); cluster_iter != this->clusters.end(); ++cluster_iter) {
        output << "\n    [Centroid " << ++ct << "]";
        output << " Address: " << &(*cluster_iter) << "   Move State: " << ((cluster_iter->move_flag) ? ("INVALID") : ("  valid"));
        output << "   No. of pts.: " << cluster_iter->count;
        output << " Dimensions: ";
        this->print_point(std::cout, cluster_iter->centroid, ',');
    }
}



// PUBLIC FUNCTIONS

// Primary algorithm
void kc::kmeans_set::compute_centroids(unsigned int k) {
    // Error handling
    if (this->universe.empty() || k == 0) {   // Exits program if no points have been imported or if no centroids are requested
        return;
    } else if (k > this->universe.size()) {   // Forces a maximum amount of centroids in a low-population universe.
        k = this->universe.size();
    }

    // Distributes centroids evenly among points
    {   // Encapsulation for variable destruction
        unsigned int subset = this->universe.size() / k;   // Used for even distribution
        unsigned int rem = this->universe.size() % k;      // Used for even distribution
        unsigned int i = 0;                                // Used for even distribution
        for (auto centroid_pick_iter = this->universe.begin(); centroid_pick_iter != this->universe.end(); ++centroid_pick_iter, --i) {
            if (i == 0) {
                cluster_wrapper new_cluster;
                new_cluster.centroid = centroid_pick_iter->point;   // Forces the dimensionality onto the centroid - this point is never actually used as a centroid
                this->clusters.push_back(new_cluster);
                i = subset + (rem ? 1 : 0);
                if (rem) { --rem; }
            }
            centroid_pick_iter->cluster_pointer = &(this->clusters.back());
            ++centroid_pick_iter->cluster_pointer->count;
        }
    }

    // Primary algo
    unsigned int n = this->hard_limit;
    double delta_max_1 = recompute_centroids();   // The first computation of centroids
    double delta_max_2 = 0;

    while ((fabs(delta_max_2 - delta_max_1) > fabs(this->epsilon)) && n && check_move_state()) {
        delta_max_2 = delta_max_1;

        // Reset move flags
        for (auto cluster_iter = this->clusters.begin(); cluster_iter != this->clusters.end(); ++cluster_iter) {
            cluster_iter->move_flag = false;
        }

        // Reassign all points to nearest centroid
        this->move_point_to_nearest_cluster();

        // Recompute all centroids, returns the highest movement (delta)
        delta_max_1 = recompute_centroids();
        --n;

        this->report_universe_state(std::cout, delta_max_1, delta_max_2);   // REPORTING - comment out for silent run

    }
}



void kc::kmeans_set::print_centroids(std::ostream & output, char delimiter) {
    for (auto cluster_iter = this->clusters.begin(); cluster_iter != this->clusters.end(); ++cluster_iter) {
        if (cluster_iter != this->clusters.begin()) {
            output << "\n";
        }
        print_point(output, cluster_iter->centroid, delimiter);
    }
}



void kc::kmeans_set::burst_clusters(std::string str, char delimiter) {
    unsigned int ct = 0;
    for (auto cluster_iter = this->clusters.begin(); cluster_iter != this->clusters.end(); ++cluster_iter) {
        std::stringstream ss;   // Avoids using std::to_string() in case compiling with older minGW
        ss << ct++;
        std::string fn = str + ss.str() + ".dat";
        std::ofstream file_out;
        file_out.open(fn);
        bool first_flag = 1;
        for (auto universe_iter = this->universe.begin(); universe_iter != this->universe.end(); ++universe_iter) {
            if (universe_iter->cluster_pointer == &(*cluster_iter)) {
                if (first_flag) {
                    print_point(file_out, universe_iter->point, delimiter);
                    first_flag = 0;
                } else {
                    file_out << "\n";
                    print_point(file_out, universe_iter->point, delimiter);
                }
            }
        }
    }
}





