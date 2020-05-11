import numpy as np
import matplotlib.pyplot as plt


threads = [1,2,3,4,5,6,7,8]

naive_200 = [1,
1.902418316,
2.770583215,
3.593385291,
4.363267402,
5.320462219,
6.086717785,
6.314405632]

bh_200 = [1,
1.707503158,
2.29308753,
2.800334601,
3.333483046,
3.800381114,
4.202925353,
4.485626752]

mpi_cores = [1,2,4,8]

naivempi_200 = [1,
1.719008124,
2.336504086,
1.84885406]



naivempi_600 = [1,
1.91667029,
3.57399645,
5.672715268]

naive_600 = [1,
1.913102955,
2.786122694,
3.713057601,
4.608205028,
5.506942101,
6.390209825,
7.155778284]

bh_600 = [0.9999999998,
1.714829618,
2.358118033,
2.985696396,
3.517673651,
4.077462303,
4.530600021,
4.918268292]

plt.plot(threads, naive_600, label='Naive, n=600')
plt.plot(threads, naive_200, label='Naive, n=200')
plt.plot(mpi_cores, naivempi_600, label='Naive MPI, n=600')
plt.plot(threads, bh_600, label='Barnes-Hut, n=600')
plt.plot(threads, bh_200, label='Barnes-Hut, n=200')
plt.plot(mpi_cores, naivempi_200, label='Naive MPI, n=200')

plt.legend()
plt.title('Speedup Ratios')
plt.xlabel('Parallel threads/cores')
plt.ylabel('Speedup')
plt.savefig('speedup_vs_threads.png')
plt.show()























# n_b =[ 100,
# 	200,
#     	300	,
#         400	,
#         500	,
#         600,
#         	1000
# ,            	1800
# ,                3000
# ,                	4000]
#
# t_b = [1.6809258,	4.468020667,	8.056646,	11.545364,	15.64071333,	19.80075633	,38.669263,	84.460101,	155.81142,	226.409964]
#
# n_n =[ 100,
# 	200,
#     	300	,
#         400	,
#         500	,
#         600,
#         	1000
# ,            	1800
# ,                3000]
# t_n =[0.480207,	1.996413,	4.286402,	7.705542,	12.110508,	17.622435,	48.121326,	154.50399,	433.893076]
#
# plt.plot(n_b, t_b, label='Barnes-Hut')
# plt.plot(n_n, t_n, label='Naive')
# plt.title('Comparing Naive and Barnes-Hut Complexities')
# plt.xlabel('simulation size')
# plt.ylabel('runtime')
# plt.legend()
# plt.savefig('barnes_hut_naive_comparison.png')
# plt.show()


# thetas = [0,
# 0.25,
# 0.5,
# 0.75,
# 1,
# 1.5,
# 2]
# times = [95.451988,
# 31.964402,
# 15.956491,
# 10.026705,
# 6.86496,
# 4.008479,
# 2.552203]
#
# times_8 = [14.603416,
# 5.542921,
# 3.129218,
# 2.347713,
# 1.917794,
# 1.535374,
# 1.562756]
#
# plt.plot(thetas, times, label='1 thread')
# plt.plot(thetas, times_8, label='8 threads')
# plt.legend()
# plt.title('Barnes-Hut precision tradeoff')
# plt.xlabel('theta threshold')
# plt.ylabel('runtime')
# plt.savefig('theta_threshold.png')
# plt.show()
