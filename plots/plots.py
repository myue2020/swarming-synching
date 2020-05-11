import numpy as np
import matplotlib.pyplot as plt

#
# threads = [1,2,3,4,5,6,7,8]
#
# naive_200 = [1,
# 1.902418316,
# 2.770583215,
# 3.593385291,
# 4.363267402,
# 5.320462219,
# 6.086717785,
# 6.314405632]
#
# bh_200 = [1,
# 1.707503158,
# 2.29308753,
# 2.800334601,
# 3.333483046,
# 3.800381114,
# 4.202925353,
# 4.485626752]
#
# mpi_cores = [1,2,4,8]
#
# naivempi_200 = [1,
# 0.7047452799,
# 0.9579013633,
# 0.7579784841]
#
#
#
# naivempi_600 = [1,
# 0.7794331933,
# 1.453401496,
# 2.306866548]
#
# naive_600 = [1,
# 1.913102955,
# 2.786122694,
# 3.713057601,
# 4.608205028,
# 5.506942101,
# 6.390209825,
# 7.155778284]
#
# bh_600 = [0.9999999998,
# 1.714829618,
# 2.358118033,
# 2.985696396,
# 3.517673651,
# 4.077462303,
# 4.530600021,
# 4.918268292]
#
# plt.plot(threads, naive_600, label='Naive, n=600')
# plt.plot(threads, naive_200, label='Naive, n=200')
# plt.plot(threads, bh_600, label='Barnes-Hut, n=600')
# plt.plot(threads, bh_200, label='Barnes-Hut, n=200')
# plt.plot(mpi_cores, naivempi_600, label='Naive MPI, n=600')
#
# plt.plot(mpi_cores, naivempi_200, label='Naive MPI, n=200')
#
# plt.legend()
# plt.title('Speedup Ratios vs Threads/Tasks')
# plt.xlabel('Parallel threads/tasks')
# plt.ylabel('Speedup')
# plt.savefig('speedup_vs_threads.png')
# plt.show()
#
#
#
#
#
#
#
n = [100	,200	,300,	400	,500,	600	,1000,	1800,	3000]

hybrid = [1.160155684,	2.95571782,	3.964797454,	4.459140537	,5.034549445,	5.332119097	,5.809260155,	6.000282181,	6.164402128]

bh_2 = [1.653943368,	1.707503158,	1.77340194,	1.710811489,	1.722160122	,1.714829618,	1.764083443	,1.819176349	,1.806516695]
bh_8 = [4.191688569,	4.485626752	,4.876284492,	4.814725046,	4.906051619,	4.918268293	,5.162892289,	5.190037806	,5.503589186]

naive_2 = [1.874402792,	1.902418316,	1.867352024,	1.914925152	,1.90491911,	1.913102955,	1.87307007,	1.893697435	,1.897263513]
naive_8 = [4.919397634,	6.314405632,	6.69967005,	7.057145579	,7.096011869,	7.155778284,	7.222033439,	7.283994136,	7.367682087]

naive_mpi_2 = [0.5008829472	,0.7047452799,	0.7210326318,	0.7466389801,	0.7627926968,	0.7794331933,	0.7792739631,	0.7739602591,	0.7803132518]

naive_mpi_4 = [0.4458779638,	0.9579013633,	1.20204626,	1.314059206,	1.376933097	,1.453401496,	1.507227445,	1.532864069	,1.561133551]


plt.plot(n, naive_8, label='Naive, 8 threads')
plt.plot(n, hybrid, label='Naive Hybrid, 2 tasks 8 threads')
plt.plot(n, bh_8, label='Barnes-Hut, 8 threads')
plt.plot(n, naive_2, label='Naive, 2 threads')
plt.plot(n, bh_2, label='Barnes-Hut, 2 threads')
plt.plot(n, naive_mpi_4, label='Naive MPI, 4 tasks')
plt.plot(n, naive_mpi_2, label='Naive MPI, 2 tasks')
plt.legend()
plt.xlabel('n')
plt.ylabel('Speedup')
plt.title('Speedup Ratios vs Simulation Size')
plt.savefig('speedup_vs_n.png')
plt.show()




















#
#
n_b =[ 100,
	200,
    	300	,
        400	,
        500	,
        600,
        	1000
,            	1800
,                3000
,                	4000]

t_b = [1.6809258,	4.468020667,	8.056646,	11.545364,	15.64071333,	19.80075633	,38.669263,	84.460101,	155.81142,	226.409964]

n_n =[ 100,
	200,
    	300	,
        400	,
        500	,
        600,
        	1000
,            	1800
,                3000]
t_n =[0.480207,	1.996413,	4.286402,	7.705542,	12.110508,	17.622435,	48.121326,	154.50399,	433.893076]

naive_8 = [0.097615,	0.316168,	0.639793,	1.091878,	1.706664,	2.462686,	6.663127,	21.211438,	58.891395]

naive_mpi_4 = [1.076992,	2.084153,	3.565921,	5.863923,	8.795277,	12.12496,	31.92705,	100.794319,	277.93463]

bh_8 = [0.401014,	0.996075,	1.65221,	2.397928,	3.188045,	4.025961	,7.489845	,16.273504,	28.310874	,40.628309]

hybrid = [0.413916,	0.675441,	1.081115,	1.728033,	2.40548,	3.304959,	8.283555,	25.749454,	70.386887]


plt.plot(n_n, t_n, label='Naive, 1 thread')
plt.plot(n_n, naive_mpi_4, label='Naive MPI, 4 tasks')
plt.plot(n_b, t_b, label='Barnes-Hut, 1 thread')
plt.plot(n_n, hybrid, label='Naive Hybrid, 2 tasks 8 threads')
plt.plot(n_n, naive_8, label='Naive, 8 threads')
plt.plot(n_b, bh_8, label='Barnes-Hut, 8 threads')

plt.title('Comparing Naive and Barnes-Hut Complexities')
plt.xlabel('Simulation size')
plt.ylabel('Runtime')
plt.legend()
plt.savefig('barnes_hut_naive_comparison.png')
plt.show()











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
