cases = {
    'WD1.0LT5.0':{'shedding_freq':0.0,'strouhal_num':0.0},
    'WD1.0LT6.0':{'shedding_freq':0.0,'strouhal_num':0.0},
#   'WD1.0LT7.0':{'shedding_freq':0.0,'strouhal_num':0.0},
    'WD1.0LT8.0':{'shedding_freq':0.0,'strouhal_num':0.0},
#   'WD2.0LT3.0':{'shedding_freq':0.0,'strouhal_num':0.0},
#   'WD2.0LT4.0':{'shedding_freq':0.0,'strouhal_num':0.0},
    'WD2.0LT5.0':{'shedding_freq':0.0,'strouhal_num':0.0},
#   'WD2.0LT6.0':{'shedding_freq':0.0,'strouhal_num':0.0},
#   'WD2.0LT7.0':{'shedding_freq':0.0,'strouhal_num':0.0},
    'WD2.0LT8.0':{'shedding_freq':0.0,'strouhal_num':0.0},
#   'WD3.0LT3.0':{'shedding_freq':0.0,'strouhal_num':0.0},
#   'WD3.0LT4.0':{'shedding_freq':0.0,'strouhal_num':0.0},
    'WD3.0LT5.0':{'shedding_freq':0.0,'strouhal_num':0.0},
    'WD3.0LT6.0':{'shedding_freq':0.0,'strouhal_num':0.0},
#   'WD3.0LT7.0':{'shedding_freq':0.0,'strouhal_num':0.0},
    'WD3.0LT8.0':{'shedding_freq':0.0,'strouhal_num':0.0}
}

u0 = 1.4776
T = 0.01

import numpy as np
import matplotlib.pyplot as plt
import scipy



for case in cases:

    filename = 'Data/{case}/export_data.csv'.format(case=case)
    data = np.loadtxt(filename, delimiter='\t', skiprows=13)

    time = data[:,0]
    Cd = data[:,1]
    Cl = data[:,4]
    N = len(time)//2
    time = time[N:]
    Cd = Cd[N:]
    Cl = Cl[N:]
    N = len(time)

    sampling_freq = time[1]-time[0]
    duration = time[-1]-time[0]

    Cl_f = abs(scipy.fft.fft(Cl))[1:N//2]
    freq = scipy.fft.fftfreq(N,sampling_freq)[1:N//2]


    peak_idx = np.argmax(Cl_f)
    shedding_freq = freq[peak_idx]

    strouhal_num = shedding_freq*T/u0

    cases[case]['shedding_freq'] = shedding_freq
    cases[case]['strouhal_num'] = strouhal_num
    print(case, strouhal_num)

for i,case in enumerate(cases):
    plt.plot([i],[cases[case]['strouhal_num']],'x',label=case)
plt.legend()
plt.show()