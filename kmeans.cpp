//
// Copyright 2017 Wyatt Miller
//

    #include "kmeans.h"
    #include <cstdlib>
    #include <limits>
    #include <sstream>
    #include <cmath>
    #include <iostream>


    // Import points from file with specified delimiter into a temporary list - calls import_points()
    kmeans_set::kmeans_set(std::ifstream & input_file, char delimiter) :
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

    // REPORTING - STARTS
    /*
        std::cout << "\nFILE INTAKE\ntemp_list ct: " << temp_point_list.size();   // Reporting
        unsigned int ct = 0;
        for(auto temp_list_iter = temp_point_list.begin(); temp_list_iter != temp_point_list.end(); ++temp_list_iter) {
            std::cout << "\n" << ++ct << ": ";
            this->print_point(std::cout, *temp_list_iter, ',');
        }
        */
    // REPORTING - ENDS

        this->import_points(temp_point_list);
    }



    // Import points from list on construction
    kmeans_set::kmeans_set(std::list<std::vector<double>> & point_list) :
            epsilon(0),
            hard_limit(std::numeric_limits<unsigned int>::max()) {
        this->import_points(point_list);
    }



    // Primary point import function - Assures dimensional integrity by comparing all intake to the first point in the [points] list
    void kmeans_set::import_points(std::list<std::vector<double>> & point_list) {
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
                while(universe_list_iter != this->universe.end() && (universe_list_iter->point < *point_list_iter)) {   // PRECEDENCE MATTERS!!
                    // Sorts points on insert - this ensures evenly distributed centroid picking
                        ++universe_list_iter;
                }
                this->universe.insert(universe_list_iter,new_point);
            }
        }

        // REPORTING
        /*
        std::cout << "\n\nIMPORTING\nuniverse ct: " << this->universe.size();   // Reporting
        unsigned int ct = 0;
        for(auto universe_iter = this->universe.begin(); universe_iter != this->universe.end(); ++universe_iter) {
            std::cout << "\n" << ++ct << ": ";
            this->print_point(std::cout, universe_iter->point, ',');
            std::cout << "   centroid ptr: " << universe_iter->centroid_pointer;
        }
         */
        // REPORTING - ENDS

    }



    // Primary algorithm
    void kmeans_set::compute_centroids(unsigned int k) {

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
                    new_cluster.centroid = centroid_pick_iter->point;
                    this->clusters.push_back(new_cluster);
                    i = subset + (rem ? 1 : 0);
                    if (rem) { --rem; }
                }
                centroid_pick_iter->centroid_pointer = &(this->clusters.back());
                ++centroid_pick_iter->centroid_pointer->count;
            }
        }

        // REPORTING
        /*
        {
            std::cout << "\n\nASSIGNING CENTROIDS\nuniverse ct: " << this->universe.size();   // Reporting
            unsigned int ct = 0;
            for (auto universe_iter = this->universe.begin(); universe_iter != this->universe.end(); ++universe_iter) {
                std::cout << "\n" << ++ct << ": ";
                this->print_point(std::cout, universe_iter->point, ',');
                std::cout << "   centroid ptr: " << universe_iter->centroid_pointer;
            }
            std::cout << "\n\ncluster ct: " << this->clusters.size();   // Reporting
            ct = 0;
            for (auto cluster_iter = this->clusters.begin(); cluster_iter != this->clusters.end(); ++cluster_iter) {
                std::cout << "\n" << ++ct << ": ";
                this->print_point(std::cout, cluster_iter->centroid, ',');
                std::cout << "   cluster ct: " << cluster_iter->count << "   Address: " << &(*cluster_iter);
            }
        }
         */
        // REPORTING

        // Primary algo
        unsigned int n = this->hard_limit;
        double delta_max_1 = std::numeric_limits<double>::max();
        double delta_max_2 = 0;

        unsigned int iter_ct = 0; // REPORTING

        while ((fabs(delta_max_2 - delta_max_1) > fabs(this->epsilon)) && n && check_move_state()) {
            delta_max_2 = delta_max_1;

            // Reset move flags
            for(auto cluster_iter = this->clusters.begin(); cluster_iter != this->clusters.end(); ++cluster_iter) {
                cluster_iter->move_flag = false;
            }

            // Reassign all points to nearest centroid
            this->move_to_nearest();

            // Recompute all centroids, returns the highest movement (delta)
            delta_max_1 = recompute_centroids();
            --n;

            // REPORTING
            {
                std::cout << "\n\nITERATION " << ++iter_ct << "   cluster ct: " << this->clusters.size();   // Reporting
                unsigned int ct = 0;
                for (auto cluster_iter = this->clusters.begin(); cluster_iter != this->clusters.end(); ++cluster_iter) {
                    std::cout << "\n"  << "Move State: " << ((cluster_iter->move_flag)?("INVALID    "):("  valid    ")) << ++ct << ": ";
                    this->print_point(std::cout, cluster_iter->centroid, ',');
                    std::cout << "   cluster ct: " << cluster_iter->count << "   Address: " << &(*cluster_iter);
                }
            }
            // REPORTING

        }
    }





    void kmeans_set::move_to_nearest() {
        for (auto universe_iter = this->universe.begin(); universe_iter != this->universe.end(); ++universe_iter) {
            for (auto cluster_iter = this->clusters.begin(); cluster_iter != this->clusters.end(); ++cluster_iter) {

                // REPORTING
                /*
                std::cout << "\nComparing: (";
                this->print_point(std::cout,universe_iter->point,',');
                std::cout << ")-(";
                this->print_point(std::cout,cluster_iter->centroid,',');
                std::cout << ") dist: ";
                std::cout << calculate_distance(universe_iter->point, cluster_iter->centroid);
                std::cout << "\n       To: (";
                this->print_point(std::cout,universe_iter->point,',');
                std::cout << ")-(";
                this->print_point(std::cout,universe_iter->centroid_pointer->centroid,',');
                std::cout << ") dist: ";
                std::cout << calculate_distance(universe_iter->point, universe_iter->centroid_pointer->centroid);
                 */
                // REPORTING

                if (calculate_distance(universe_iter->point, cluster_iter->centroid) < calculate_distance(universe_iter->point, universe_iter->centroid_pointer->centroid)) {

                    // REPORTING
                    /*
                    std::cout << "\n-----------------------------------------\nSending: (";
                    this->print_point(std::cout, universe_iter->point, ',');
                    std::cout << ") From: " << universe_iter->centroid_pointer;
                     */
                    // REPORTING

                    --universe_iter->centroid_pointer->count;
                    universe_iter->centroid_pointer->move_flag = true;
                    universe_iter->centroid_pointer = &(*cluster_iter);  // This is why I'm using lists vs. vectors - iterator would be invalid
                    ++universe_iter->centroid_pointer->count;
                    universe_iter->centroid_pointer->move_flag = true;

                    // REPORTING
                    /*
                    std::cout << "    To: " << universe_iter->centroid_pointer << "\n-----------------------------------------";
                    // REPORTING
                     */
                }
            }
        }
    }



    double kmeans_set::calculate_distance(std::vector<double> & a_vec, std::vector<double> & b_vec) {
        auto b_iter = b_vec.begin();
        double r_val = 0;
        for (auto a_iter = a_vec.begin(); a_iter != a_vec.end(); ++a_iter, ++b_iter) {
            r_val += (*a_iter - *b_iter) * (*a_iter - *b_iter);
        }
        // return std::sqrt(r_val);
        return r_val;   // No need to root - actual distance not necessary
    }



    double kmeans_set::recompute_centroids() {

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
            if (universe_iter->centroid_pointer->move_flag) {
                auto centroid_pointer_centroid_iter = universe_iter->centroid_pointer->centroid.begin();
                for (auto point_iter = universe_iter->point.begin(); point_iter != universe_iter->point.end(); ++point_iter, ++centroid_pointer_centroid_iter) {
                    // Division on the fly useful for accuracy super-large dimensions - doubles computation
                    *centroid_pointer_centroid_iter += *point_iter / universe_iter->centroid_pointer->count;
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



    void kmeans_set::print_centroids(std::ostream & output, char delimiter) {
        for(auto cluster_iter = this->clusters.begin(); cluster_iter != this->clusters.end(); ++cluster_iter) {
            for (auto centroid_iter = cluster_iter->centroid.begin(); centroid_iter != cluster_iter->centroid.end(); ++centroid_iter) {
                if(cluster_iter != this->clusters.begin() && centroid_iter == cluster_iter->centroid.begin()) {
                    output << "\n";
                } else if(centroid_iter != cluster_iter->centroid.begin()) {
                    output << delimiter;
                }
                output << *centroid_iter;
            }
        }
    }



    void kmeans_set::print_point(std::ostream & output, std::vector<double> vec, char delim) {
        for(auto iter = vec.begin(); iter != vec.end(); ++iter) {
            if(iter != vec.begin()) { output << delim; }
            output << *iter;
        }
    }



    bool kmeans_set::check_move_state() {
        for (auto cluster_iter = this->clusters.begin(); cluster_iter != this->clusters.end(); ++cluster_iter) {
            if (cluster_iter->move_flag) { return true; }
        }
        return false;
    }