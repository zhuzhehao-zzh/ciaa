# ciaa_implement

C++ implementation of ciaa and baseline

## Structure

* `CIAA.cpp`: implementation of ciaa algorithm
  * `class point`: class representing a point, includes its properties and the cluster it belongs to
  * `class cluster`: a cluster, contains points
  * `class CIAA`: main implementation of the algorithm
  * `BuildGraph, quicksort, Calmatch, Calmismatch, Selfdiffusion,...`: some helper functions, also used to calculate the mismatch score and match score
  * able to show the running time
* `base_medoids.cpp`: implementation of k-medoids
* `base_means.cpp`: implementation of k-means

## Data

* `dataset_Iris.txt`: Iris dataset
* `dataset_linear`: a simple dataset
* `dataset_sample`: artificially generated data

## Environment

* C++11
* gcc: 11.3.0
* CPU: AMD Ryzen 7 5800H with Radeon Graphics
* RAM: 32GB
