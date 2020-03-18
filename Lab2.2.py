import random
from math import sin, cos, pi
import matplotlib.pyplot as plot
import timeit

N = 64
n = 12
freq = 1100


def generate_signal():
    xS_arr = [0] * N
    for i in range(1, n + 1):
        A = random.random()
        phase = random.random()
        for j in range(N):
            xS_arr[j] += A * sin(freq/i * (j + 1) + phase)
    return xS_arr


# function of DFT for comparing with FFT
def get_discrete_fourier_transform(signal: list):
    length = len(signal)
    return [sum((signal_X[j] * complex(cos(-2*pi*i*j/length), sin(-2*pi*i*j/length)) for j in range(length))) for i in range(N)]


def get_fast_fourier_transform(signal, length_of_s, start_count=0, step_count=1):
    # try:
    #     assert length_of_s == 1, "Signal is from 1 element!"
    # except AssertionError:
    #     return [signal[start_count]]
    if length_of_s == 1:
        return [signal[start_count]]
    hns, sns = length_of_s//2, step_count*2
    rss = get_fast_fourier_transform(signal, hns, start_count, sns) + get_fast_fourier_transform(signal, hns, start_count + step_count, sns)
    for p in range(hns):
        cons_e = complex(cos(-2*pi*p/length_of_s), sin(-2*pi*p/length_of_s))
        rss[p], rss[p + hns] = rss[p] + cons_e * rss[p + hns], rss[p] - cons_e * rss[p + hns]
    return rss


if __name__ == "__main__":
    signal_X = generate_signal()
    time_start = timeit.default_timer()
    dft = get_discrete_fourier_transform(signal_X)
    print("Calculating DFT time: {}".format(timeit.default_timer() - time_start))
    time_start = timeit.default_timer()
    fft = get_fast_fourier_transform(signal_X, len(signal_X))
    print("Calculating FFT time: {}".format(timeit.default_timer() - time_start))
    figure, ((xy1, xy2), (xy3, xy4)) = plot.subplots(2, 2, figsize=(9, 9))
    xy1.plot(range(N), signal_X, "b")
    xy1.title.set_text("Signal")
    xy2.plot(range(N), dft, "r")
    xy2.title.set_text("DFT")
    xy3.plot(range(N), signal_X, "b")
    xy3.title.set_text("Signal")
    xy4.plot(range(N), fft, "r")
    xy4.title.set_text("FFT")
    plot.show()