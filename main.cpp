//
// Copyright 2017 Wyatt Miller
//



#include <iostream>
#include "kmeans.h"

int main() {

    std::ifstream input;
    input.open("test_data_1.dat");
    std::ostream &output_console = std::cout;
    std::ofstream file_out;
    file_out.open("centroids.dat");


    kc::kmeans_set km(input, ',');
    km.compute_centroids(5);

    km.print_centroids(file_out, ',');
    km.burst_clusters("cl_", ',');

}