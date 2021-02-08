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
		return AnomFunc(self.n, (60 / self.bpm), t, anom[anomalia]) , t

	@staticmethod
	def add_noise(signal, fre, t):
		""" Método para añadir ruido """
		return NoiseAdd(signal, fre, t)

	@staticmethod
	def fourier(signal, t):
		""" Método para representar espectro de fourier """
		ECGFourier = np.fft.fft(signal[0,:])
		Freq = np.fft.fftfreq(t.shape[-1]) * 100
		EjX=max(Freq)
		EjY=max(ECGFourier)
		fig = plt.figure(figsize=(10, 10))
		plt.xlim(0, EjX)
		plt.ylim(0, EjY)
		plt.stem(Freq, ECGFourier.real)
		plt.show()

	@staticmethod
	def show_ecg(signal, t):
		""" Método para mostrar derivadas """
		name = ['I','II','III', 'avR', 'avF', 'avL', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6']

		fig = plt.figure(figsize=(10, 10))
		for i in range(0, 12):
			fig.tight_layout(pad=2.0)
			ax = plt.subplot(3, 4, i + 1)
			ax.plot(t, signal[i])
			ax.set_xlabel('x')
			ax.set_ylabel('y')
			ax.set_title('Derivada:' + name[i])
		plt.show()


########################################### Creación de objeto ####################################################

########################################## ECG NORMAL #############################################################

mio = Ecg(70, 60, 5) ### armónicos se deja fijo
array, t = mio.derivadas()
Ecg.show_ecg(array, t)

######################################### ECG ANOMALÍA ############################################################

array2, t2 = mio.anomalias('Fibrilación auricular')
Ecg.show_ecg(array2, t2)

######################################### ECG CON RUIDO ###########################################################
withnoise = Ecg.add_noise(array, [1, 60], t)
Ecg.show_ecg(withnoise, t)

######################################## Espectro de Fourier ######################################
Ecg.fourier(array, t)
