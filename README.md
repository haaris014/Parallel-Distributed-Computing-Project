# Parallel-Distributed-Computing-Project
MapReduce, Scatter, and Gather Concepts for Randomly Generated Array Processing using MPI
This brief report explores the application of MapReduce, Scatter, and Gather concepts in the context of processing a randomly generated array using the Message Passing Interface (MPI) on a Beowulf cluster. By leveraging these concepts, we can efficiently distribute and process large-scale data sets across a cluster of interconnected computers.

MapReduce:
MapReduce is a programming model and framework designed for processing large-scale data sets in parallel across a distributed system. In the context of processing a randomly generated array using MPI, MapReduce can be applied as follows:
Map Phase: Each node in the Beowulf cluster receives a subset of the randomly generated array. Using MPI, each node performs a mapping operation on its subset, applying a specified function or operation to transform the data.

Reduce Phase: The results from the Map phase are combined and reduced to produce a single output or aggregated result. This step involves a communication and synchronization process across the cluster to consolidate the data and compute the final result.

Scatter:
Scatter is a data distribution technique that allows for the partitioning and distribution of data across multiple processes or nodes in a parallel computing environment. In the context of processing a randomly generated array using MPI, Scatter can be applied as follows:
Data Partitioning: The randomly generated array is divided into equal-sized chunks, and each chunk is assigned to a specific process in the Beowulf cluster.

Data Distribution: Using MPI Scatter operations, the master process scatters the chunks of the array to the individual processes, ensuring that each process receives a unique subset of the data.

Gather:
Gather is a data collection technique that enables the gathering of data from multiple processes or nodes in a parallel computing environment. In the context of processing a randomly generated array using MPI, Gather can be applied as follows:
Data Collection: After performing the desired computations on their assigned subsets of the array, each process uses the MPI Gather operation to send their results back to the master process.

Data Aggregation: The master process collects the results from all the processes using the MPI Gather operation, combining the individual results into a single data structure or array.
