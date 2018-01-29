
# coding: utf-8

# In[1]:

#importing libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import square


# ## Creating a class for modulation

# In[2]:

class modulation(object):
    #receives the signal and the frequecy of the carrier and modulates(FM and PM) 
    def __init__(self,frequency_carrier=200., t_max=1., kind='fm'): 
        self.frequency_carrier = frequency_carrier
        self.t_max = t_max
        self.kind = kind
        #Freq should be im Mhz, but I'm dont have this kind of prcessing power, it must be scaled, but...
        c = 10. #constant that makes things "continuous"
        self.t = t = np.linspace(0, self.t_max,  self.t_max * frequency_carrier * c)
        self.mod = np.zeros(self.t.shape[0])
        self.sig = np.zeros(self.t.shape[0])
        if (self.kind not in["fm", "pm"]): 
            raise  NameError('%s is not implemented and problably will never be, deal with it!'%(self.kind))
    
    def _mod_fm(self,kfm, amplitude):
        print "#########################FM Modulation############################"
        return amplitude*np.cos(self.frequency_carrier*2*np.pi*self.t 
                                + 2*np.pi*kfm*np.cumsum(self.sig)) #the jump of the cat
    
    def _mod_pm(self, kpm, amplitude):
        print "#########################PM Modulation############################"
        return amplitude*np.cos(self.frequency_carrier*2*np.pi*self.t + kpm*self.sig)
    
    def modulate(self, signal, k=1., amplitude=1., showing_options='both', periods_to_show=10):
        self.sig = signal
        if (self.kind=='fm'): self.mod = self._mod_fm(k, amplitude)
        else: self.mod = self._mod_pm(k, amplitude)
        
        if (showing_options in['time', 'both']):
            plt.subplot(211)
            plt.title("Signal/Modulated signal (%s)"%(self.kind.upper()))
            plt.plot(self.t, self.sig)
            plt.xlim(0, periods_to_show/float(self.frequency_carrier))
            
            plt.subplot(212)
            plt.plot(self.t, self.mod)
            plt.xlim(0, periods_to_show/float(self.frequency_carrier))
            plt.xlabel("Time(s)")
            plt.ylabel("Amplitude(V)")
            plt.show()
            
        if (showing_options in ["frequency", "freq", "both"]):
            fs = self.t.shape[0]/self.t_max #another way to retrieve the frequency 
            f = np.linspace(-fs/2.,fs/2.,self.mod.shape[0])
            M = np.fft.fftshift(np.abs(np.fft.fft(self.mod)))
            plt.title("Signal in frequency (%s)"%(self.kind.upper()))
            plt.plot(f, M)
            plt.xlabel("Frequency(Hz)")
            plt.ylabel("Absolute value")
            plt.grid()
            plt.show()

    def demodulate(self, frequency, periods_to_show=10):
        print "#########################%s Demodulation############################"%(self.kind.upper())
        s_diff = np.diff(np.hstack((self.mod[1], self.mod)))
        s_diode = np.zeros(self.sig.shape[0])
        for i in xrange(self.sig.shape[0]):
            if (self.sig[i]>=0): s_diode[i] = self.sig[i]
        #ideal low pass, very well compressed in two line
        f_lp = int(round(1.5*frequency*self.t_max)) # (f/fs)*samples[], fs = samples[]/t and 10% more
        s_filtered = np.fft.ifft(np.multiply(np.fft.fft(s_diode),
                                             np.hstack((np.ones(f_lp), 
                                                        np.zeros(s_diode.shape[0]-2*f_lp),
                                                        np.ones(f_lp))))).real
        
        s_nodc = s_filtered - s_filtered.mean() # removing DC in the easiest way
        plt.subplots_adjust(hspace=.7, wspace=.7)#adjusting spacing
        plt.subplot(211)
        plt.title("Original signal")
        plt.plot(self.t, self.sig)
        plt.xlim(0, periods_to_show/float(self.frequency_carrier))
        
        plt.subplot(212)
        plt.title("Demodulated signal")
        plt.plot(self.t, 2*s_nodc)
        plt.xlim(0, periods_to_show/float(self.frequency_carrier))
        plt.xlabel("Time(s)")
        plt.ylabel("Amplitude(V)")
        plt.show()
        
        print "#########################Showing steps############################"
        plt.subplots_adjust(hspace=1., wspace=1.)#adjusting spacing
        plt.subplot(411)
        plt.title("Differentiating")
        plt.plot(self.t, s_diff)
        plt.xlim(0, periods_to_show/float(self.frequency_carrier))
       
        plt.subplot(412)
        plt.title("'Diode'")
        plt.plot(self.t, s_diode)
        plt.xlim(0, periods_to_show/float(self.frequency_carrier))
       
        plt.subplot(413)
        plt.title("Ideal Low pass (+50%)")
        plt.plot(self.t, s_filtered)
        plt.xlim(0, periods_to_show/float(self.frequency_carrier))
      
        plt.subplot(414)
        plt.title("Removing DC level")
        plt.plot(self.t, s_nodc)
        plt.xlim(0, periods_to_show/float(self.frequency_carrier))
        plt.show()


# ## Testing the modulations

# ### FM

# In[3]:

# creating things
# "Let there be light..."  Just kidding
freq1 = 20.
fm = modulation(kind='fm')
signal = np.cos(freq1*2*np.pi*fm.t)
fm.modulate(signal, periods_to_show=20, k=.03, amplitude=2.)


# ### PM

# In[4]:

freq2 =40.
pm = modulation(kind='pm')
signal2 = np.sin(freq2*2*np.pi*pm.t)
pm.modulate(signal2, periods_to_show=10, k=0.5*np.pi, amplitude=7.)


# ## Testing the demodulations
# 

# ### FM

# In[5]:

fm.demodulate(frequency=freq1, periods_to_show=20)


# ### PM

# In[6]:

pm.demodulate(frequency=freq2, periods_to_show=10)


# # Square wave

# ### FM Modulation 

# In[13]:

fq = 10.
sqfm = modulation(kind='fm')
sq = square(2 * np.pi * fq * sqfm.t)
sqfm.modulate(signal=sq, k=.0001, periods_to_show=40 )


# ### FM Demodulation

# In[14]:

sqfm.demodulate(frequency=fq+40, periods_to_show=40)


# ### PM Modulation

# In[9]:

fq2 = 10.
sqpm = modulation(kind='pm')
sq2 = square(2 * np.pi * fq2 * sqpm.t)
sqpm.modulate(signal=sq2, k=5*np.pi, periods_to_show=40 )


# ### PM Demodulation

# In[11]:

sqpm.demodulate(frequency=fq2+10, periods_to_show=40)


# In[ ]:



