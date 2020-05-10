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



### Background
The swarm behavior is the collective motion of a large number of self-propelled entities. Many swarming systems in Nature show the remarkable ability to spontaneously fall into synchrony, giving rise to many bio-inspired models. Many researchers have studied the close relation between swarming and synchronization, which interestingly represent two effects that stand as "spatiotemporal opposites".  
N-body problem and is traditionally compute intensive

<table>
<tr>
<td><img src="Images/refs/bees.jpg"/></td>
<td><img src="Images/refs/frogs.jpg"/></td>
</tr>
</table>

### Description
Naively, the Swarming-Synching model can be simulated by pairwise calculation of aggregation and synchronization forces of individual points and adding up all such contributions on all the entities in the system. Such an approach has a quadratic time complexity and scales up very quickly. Most of the time, when we deal with realistic 

### Infrastructure
### Dependents
### Implementation
#### Naive Algorithm
<img src="Images/refs/pairwise.png"/>

#### Barnes-Hut Algorithm
<img src="Images/refs/barnes8.png"/>

### Example
#### Plotting Into Graph
This is example using matplotlib to plot graphs

### Runtime Analysis

### References
1. O’Keeffe and Bettstetter. *A review of swarmalators and their potential in bio-inspired computing}*  https://arxiv.org/pdf/1903.11561.pdf. 2019. 
2. O’Keeffe, Hong, and Strogatz. *Oscillators that sync and swarm*,  https://www.nature.com/articles/s41467-017-01190-3. 2017.
3. Gan and Xu. *Efficient Implementation of the Barnes-Hut Octree Algorithm for Monte Carlo Simulations of Charged Systems*  https://arxiv.org/pdf/1305.1825.pdf. 2013.
