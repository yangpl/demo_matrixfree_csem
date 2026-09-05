import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['axes.xmargin'] = 0

it1, fk1 = np.loadtxt('iterate_fbicgstab_freq1_airbc1.txt', skiprows=1, unpack=True)
it2, fk2 = np.loadtxt('iterate_fbicgstab_freq2_airbc1.txt', skiprows=1, unpack=True)
it3, fk3 = np.loadtxt('iterate_fbicgstab_freq3_airbc1.txt', skiprows=1, unpack=True)


it1_, fk1_ = np.loadtxt('iterate_fbicgstab_freq1_airbc0.txt', skiprows=1, unpack=True)
it2_, fk2_ = np.loadtxt('iterate_fbicgstab_freq2_airbc0.txt', skiprows=1, unpack=True)
it3_, fk3_ = np.loadtxt('iterate_fbicgstab_freq3_airbc0.txt', skiprows=1, unpack=True)

plt.plot(it1, fk1, 'r', label='Air-Earth BC: 0.25 Hz')
plt.plot(it2, fk2, 'g', label='Air-Earth BC: 0.75 Hz')
plt.plot(it3, fk3, 'b', label='Air-Earth BC: 1.25 Hz')

plt.plot(it1_, fk1_, 'r--', label='Dirichlet BC: 0.25 Hz')
plt.plot(it2_, fk2_, 'g--', label='Dirichlet BC: 0.75 Hz')
plt.plot(it3_, fk3_, 'b--', label='Dirichlet BC: 1.25 Hz')

plt.xlabel('# iteration k')
plt.ylabel('$|r|/|r_0|$')
plt.yscale('log')
plt.legend()
#plt.grid(True)
#plt.ylim([1e-8,1e3])

#plt.xticks(np.arange(0, 20, 2))
plt.savefig('convergence.png', bbox_inches='tight')
plt.show()


