from obspy import read
import matplotlib.pyplot as plt

st = read('BK.FARB.00.50.0000.ZDD', debug_headers=True)
print(st[0].data)
plt.plot(st[0].data)
plt.show()
