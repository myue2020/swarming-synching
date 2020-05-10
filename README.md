# Swarming and Synchronization with Big Compute O'Keefe Model Implementations

Hybrid implementation of two algorithms based on O'Keefe's Swarmalator model.

- Description of problem and the need for HPC and/or Big Data
- Description of solution and comparison with existing work on the problem
- Description of your model and/or data in detail: where did it come from, how did you acquire it, what does it mean, etc.
- Technical description of the parallel application, programming models, platform and infrastructure
- Links to repository with source code, evaluation data sets and test cases
- Technical description of the software design, code baseline, dependencies, how to use the code, and system and environment needed to reproduce your tests
- Performance evaluation (speed-up, throughput, weak and strong scaling) and discussion about overheads and optimizations done
- Description of advanced features like models/platforms not explained in class, advanced functions of modules, techniques to mitigate overheads, challenging parallelization or implementation aspects...
- Final discussion about goals achieved, improvements suggested, lessons learnt, future work, interesting insights…
- Citations

Your web page should include screenshots of your software that demonstrate how it functions. You should include a link to your source code.

<table>
<tr>
<td><img src="Images/refs/bees.jpg"/></td>
<td><img src="Images/refs/frogs.jpg"/></td>
</tr>
</table>


### Background
The swarm behavior is the collective motion of a large number of self-propelled entities. Many swarming systems in Nature show the remarkable ability to spontaneously fall into synchrony, giving rise to a large number of bio-inspired models. Many researchers have studied the close relation between swarming and synchronization, which interestingly represent two effects that stand as "spatiotemporal opposites". Among them, O'Keefe et al proposed a ‘bottom-up’ models without reference to a background medium, which they called the ‘swarmalators’, to capture their twin identities as swarming oscillators. In the swarmalator model, the paradigmatic model of biological aggregation and synchronization have the following forms:
<img src="Images/refs/formulae.png"/ width = "300" />

Here, the terms *K* and *J* capture the tendency to synchronize/desynchronzie and the spatial attraction/repulsion between entites of similar phases respectively. A combination of different *K* and *J* values yield a series of different swarmalator states shown as follows:

<img src="Images/refs/states.ong"/ width = "600" />

The swarmalator model simulation intuitively falls into the category of a special kind of N-body problem, which is traditionally compute intensive. The dynamic scale that needs to be resolved for studying a real-world swarming system in a self-consistent manner is enormous and spans many orders of magnitudes, thus necessitating the use of high performance computing and carefully tailored codes that can utilize various HPC programming models.

### Description
Naively, the Swarming-Synching model can be simulated by pairwise calculation of aggregation and synchronization forces of individual points and adding up all such contributions on all the entities in the system. Such an approach has a quadratic time complexity and scales up very quickly. Most of the time, when we deal with realistic problems, the entities in a swarm system can be hundreds of millions. Under these scenarios, the quadratic scaling of the naive algorithm is not feasible and approximate models using tree-based data structures, for example, the Barnes-Hut algorithm is often utilized to facilitate the computation. As far as big compute techniques are concerned, there is large potential for these algorithms to benefit from parallelization methods such as OpenMP and MPI and as the two algorithms use different data structure, there are differences in where the parallelization can take place and eventually contribute to an overall speedup. Our project aims at investigating the performance of parallelized implementation of both the naive algorithm and the Barnes-Hut algorithm that address the same swarming-synching model, to look into the potential of parallelization in both models and to compare the consequent speedups.

### Github Repository
https://github.com/myue2020/swarming-synching

### Infrastructure


### Dependents
### Implementation
#### Naive Algorithm
<img src="Images/refs/pairwise.png" width="400"/>

#### Barnes-Hut Algorithm
1. Quadtree
2. 
3. Integration
<img src="Images/refs/barnes8.png" width="800"/>

### Example
#### Plotting Into Graph
This is example using matplotlib to plot graphs

### Runtime Analysis

### References
1. O’Keeffe and Bettstetter. *A review of swarmalators and their potential in bio-inspired computing}*  https://arxiv.org/pdf/1903.11561.pdf. 2019. 
2. O’Keeffe, Hong, and Strogatz. *Oscillators that sync and swarm*,  https://www.nature.com/articles/s41467-017-01190-3. 2017.
3. Gan and Xu. *Efficient Implementation of the Barnes-Hut Octree Algorithm for Monte Carlo Simulations of Charged Systems*  https://arxiv.org/pdf/1305.1825.pdf. 2013.
