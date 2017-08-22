//
// Copyright 2017 Wyatt Miller
//



#include <iostream>
#include "kmeans.h"

int main() {

    std::ifstream input;
    input.open("test_data_7.dat");
    std::ostream &output_console = std::cout;


    kmeans_set km(input, ',');
    km.compute_centroids(5);

    std::cout << "\n\nOutput in main():\n";
    km.print_centroids(output_console, ',');


}