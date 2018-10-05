from pytriqs.gf import *
from pytriqs.gf.tools import tail_fit
import numpy as np

# Read in F(i\omega_n)
f_iw_table = np.loadtxt("selfenergy_17")
freqs = f_iw_table[:,0]
n_iw = len(freqs)
beta = np.pi / freqs[0]

n_tau = 2049

block_names = ("up","dn")

# F in Matsubara frequency representation
F_iw = BlockGf(name_block_generator = [(bn, GfImFreq(indices=[0], beta = beta, n_points = n_iw)) for bn in block_names],
               make_copies = True
)

F_iw["up"].data[n_iw:,0,0] = f_iw_table[:,1] + 1j*f_iw_table[:,2]
F_iw["up"].data[n_iw-1::-1,0,0] = np.conj(F_iw["up"].data[n_iw:,0,0])
F_iw["dn"].data[n_iw:,0,0] = f_iw_table[:,3] + 1j*f_iw_table[:,4]
F_iw["dn"].data[n_iw-1::-1,0,0] = np.conj(F_iw["dn"].data[n_iw:,0,0])

# Fit tails
for bn, f_iw in F_iw:
  f_iw.fit_tail(known_moments = TailGf(1,1,0),
                max_moment = 3,
                n_min = 400,
                n_max = 1024,
                replace_by_fit = False)
  print "Norm =", f_iw.tail[1][0,0].real


# F in imaginary time representation
F_tau = BlockGf(name_block_generator = [(bn, GfImTime(indices=[0], beta = beta, n_points = n_tau)) for bn in block_names],
                make_copies = True
)

# Transform!
for bn, f_iw in F_iw:
  F_tau[bn] << InverseFourier(f_iw)

# Save results
f_tau_table = np.transpose(np.vstack((list(F_tau["up"].mesh),
                                      F_tau["up"].data[:,0,0].real,
                                      F_tau["up"].data[:,0,0].imag,
                                      F_tau["dn"].data[:,0,0].real,
                                      F_tau["dn"].data[:,0,0].imag)))
np.savetxt("selfenergy_17.tau", f_tau_table)
