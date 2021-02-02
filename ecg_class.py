import numpy as np
import matplotlib.pyplot as plt
from ECGTot import ECGDualGuide as ECGen
from NoiseAdd import NoiseAdd
from ECGAnomalias import AnomFunc

class Ecg():
	""" ECG Class """
	def __init__(self, bpm, n, amplitud):
		""" Initializated method 
		Args:
			bpm: frecuencia cardiaca.
			n: número de armónicos.
		Attr:
			frec: frecuencia por seg.
			T: periodo.
			TPW, TW, QRS, PQ, ST: duración de ondas.
			DsfP & T: desfase de ondas.
			A, B, C, D: Limites de integración.
		"""
		self.bpm = bpm
		self.n = n
		self.amplitud = amplitud

	def derivadas(self):
		""" Método para obtener derivadas """
		t = np.linspace(0.0, 3.0 * (60 / self.bpm), int(3.0 * (60 / self.bpm) * 100.0))
		return ECGen(self.amplitud, self.bpm, self.n, t), t

	def anomalias(self, anomalia):
		""" Método para obtener """
		anom = {'Fibrilación auricular': 2,
			'Perdida de latido': 8,
			'Bradicarida extrema': 11,
			'Taquicardia extrema': 12,
			'Flúter auricular': 13,
			'Taquicárdia paroxistica': 14,
			'Ritmo nodal': 15,
			'Taquicardia supraventricular': 16,
			'Taquicardia ventricular': 17
		}
		t = np.linspace(0.0, 3.0 * (60 / self.bpm), int(3.0 * (60 / self.bpm) * 100.0))
		return AnomFunc(self.n, (60 / self.bpm), t, anom['Fibrilación auricular'])

	@staticmethod
	def add_noise(signal, fre, t):
		""" Método para añadir ruido """
		return NoiseAdd(signal, fre, t)

########################################### Creación de objeto ####################################################

mio = Ecg(90, 60, 180)
array, t = mio.derivadas()
withnoise = Ecg.add_noise(array, [1, 60], t)

name = ['I','II','III', 'avR', 'avF', 'avL', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6']

fig = plt.figure(figsize=(10, 10))
fig.tight_layout()
fig.subplots_adjust(left=0.125, bottom=0.1, right=None, top=0.9, wspace=None, hspace=None)
for i in range(0, 12):
	ax = plt.subplot(3, 4, i + 1)
	ax.plot(t, withnoise[i])
	ax.set_xlabel('x')
	ax.set_ylabel('y')
	ax.set_title('Derivada:' + name[i])
plt.show()
