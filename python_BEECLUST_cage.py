# this file should contain adapted "work-in-progress" SAGE code.
# thomas uses some sort of self written pulse function that seems to work only in the SAGE Code. 
# the internet has some very simple solutions.
# libraries in general don't provide built-in signal functions due to one-line-implementations.

# following lines are ways to create a pulse:
# ver 1 
from scipy import signal
import matplotlib.pyplot as plt
t = np.linspace(0, 1, 500, endpoint=False)
plt.plot(t, signal.square(2 * np.pi * 5 * t))
plt.ylim(-2, 2)

# ver 2 
N = 100 # sample count
P = 10  # period
D = 5   # width of pulse
sig = np.arange(N) % P < D
plot(sig)

# ver 3 (fully configurable)
def rect(T):
    """create a centered rectangular pulse of width $T"""
    return lambda t: (-T/2 <= t) & (t < T/2)

def pulse_train(t, at, shape):
    """create a train of pulses over $t at times $at and shape $shape"""
    return np.sum(shape(t - at[:,np.newaxis]), axis=0)

sig = pulse_train(
    t=np.arange(100),              # time domain
    at=np.array([0, 10, 40, 80]),  # times of pulses
    shape=rect(10)                 # shape of pulse
)
