import numpy as np
import matplotlib.pyplot as plt

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


thetas = [0,
0.25,
0.5,
0.75,
1,
1.5,
2]
times = [95.451988,
31.964402,
15.956491,
10.026705,
6.86496,
4.008479,
2.552203]

times_8 = [14.603416,
5.542921,
3.129218,
2.347713,
1.917794,
1.535374,
1.562756]

plt.plot(thetas, times, label='1 thread')
plt.plot(thetas, times_8, label='8 threads')
plt.legend()
plt.title('Barnes-Hut precision tradeoff')
plt.xlabel('theta threshold')
plt.ylabel('runtime')
plt.savefig('theta_threshold.png')
plt.show()
