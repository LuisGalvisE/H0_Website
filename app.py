# -- Use Agg backend to be thread safe
import matplotlib as mpl
mpl.use("agg")

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import requests, os
from gwpy.timeseries import TimeSeries
from gwosc.locate import get_urls
from gwosc import datasets
from gwosc.api import fetch_event_json
from copy import deepcopy
import io
from scipy import signal
from scipy.io import wavfile
from freqdomain2 import showfreqdomain

# -- Need to lock plots to be more thread-safe
from matplotlib.backends.backend_agg import RendererAgg
lock = RendererAgg.lock

# -- Helper functions in this git repo
from helper import *

apptitle = 'Signal Processing Tutorial'

st.set_page_config(page_title=apptitle, page_icon=":headphones:",
                               initial_sidebar_state='collapsed')

# Title the app
st.title(apptitle)

fs = 32000
noisedt = 8
noise = deepcopy(makewhitenoise(fs, noisedt))

#-- Try to color the noise
noisefreq = noise.fft()
color = 1.0 / (noisefreq.frequencies)**2
indx = np.where(noisefreq.frequencies.value < 30)
color[indx] = 0  #-- Apply low frequency cut-off at 30 Hz

#-- Red noise in frequency domain
weightedfreq = noisefreq * color.value

# -- Try returning to time domain
colorednoise = weightedfreq.ifft()

###
# -- Inject the signal
###
secret = TimeSeries.read('LOZ_Secret.wav')

# -- Normalize and convert to float
secret -= secret.value[0]  #-- Remove constant offset
secret = np.float64(secret)
secret = secret/np.max(np.abs(secret)) * 1*1e-8   #-- Set amplitude
secret.t0 = 4

volume = st.sidebar.radio("Secret sound volume", ["Default", "Louder"])

if volume == 'Louder':
    maze = colorednoise.inject(10*secret)
else:
# -- Might be useful to make easier to hear option
    maze = colorednoise.inject(secret)


# -------
# Begin Display Here
# -------
st.markdown("## Introduction")

st.markdown("""
In this demo, we will try to find a **secret sound** hidden in noisy
data.  To do this, we will practice with a few signal processing concepts:

 * Plotting in the time domain and frequency domain
 * Highpass and bandpass filtering
 * Whitening
""")

sectionnames = [
                'Introduction to the frequency domain',
                'White Noise',
                'Red Noise',
                'Find the Secret Sound',
                'Whitening',
                'Gravitational Wave Data',
]

def headerlabel(number):
    return "{0}: {1}".format(number, sectionnames[number-1])
    
page = st.radio('Select Section:', [1,2,3,4,5,6], format_func=headerlabel)

st.markdown("## {}".format(headerlabel(page)))

if page==1:
    
    showfreqdomain()
    
if page==2:

    # White Noise
    
    st.markdown("""
    Next, let's take a look at some **white noise**.  Any 
    signal can be represented based on its frequency content.  When we 
    say that noise is *white*, we mean the signal has about the same 
    amplitude at all frequencies.  
    
    Below, we'll represent the **same signal three different ways**:
    
    * A time-domain signal
    * A frequency-domain signal
    * An audio file
    """)

    st.markdown("### Time domain")

    st.markdown("""
    In the **time domain**, we see a signal as a function of time.  The 
    x-axis represents time, and the y-axis represents the value of
    the signal at each time.  For an audio signal, the signal value 
    corresponds to the amount of pressure felt on your eardrum at any 
    moment.  For a 
    gravitatonal-wave signal, the signal value represents the strain - 
    or fractional change in length - of the observatory's arms.
    """)

    with lock:
        tplot = noise.plot(ylabel='Pressure')
        st.pyplot(tplot)
    
    st.markdown("### Frequency domain")

    st.markdown("""
    In the **frequency domain**, the x-axis represents a frequency 
    value, and the y-axis shows the 
    **amplitude**,
    or the closely related amplitude spectral density,
    of the signal at each
    frequency.  Since white noise has about the same amplitude at each 
    frequency, this plot is mostly flat as you move from left to right.
    """)

    with lock:
        figwn = noise.asd(fftlength=1).plot(ylim=[1e-10, 1], ylabel='Amplitude Spectral Density')
        st.pyplot(figwn)

    st.markdown("### Audio player")
    st.markdown("""
    :point_right: **Use the audio player to listen the signal.  You should hear
    a hiss of white noise**.
    """)
    
    st.audio(make_audio_file(noise), format='audio/wav')

    st.markdown("")
    st.markdown("""
    When ready, go to the next section using the controls at the 
    top.
    """)
    
if page == 3:

    # st.markdown("## 3: Red Noise")
    
    st.markdown("""
    Next, we'll look at some **red noise**.  Red noise 
    has more power at low frequencies than high frequencies.
    
    Imagining random noise at different frequencies can be a hard thing
    to understand.  A silly way to picture this is as a sports stadium
    full of animals cheering. Some animals (like birds and kittens)
    cheer with higher pitches, and other animals (like bullfrogs and 
    lions) will cheer with lower pitches.  If the stadium has animals of 
    all kinds in equal numbers, you might get white noise cheering.  If 
    the stadium is full of low pitch creatures (say, lots of bullfrogs), 
    you might get red noise cheering.  Can you imagine the difference?
    
    A similar idea can be seen in noise in the LIGO and Virgo instruments.
    Low frequency noise sources contribute noise at low frequencies.  These 
    are big, slowly vibrating things, especially motion from the constant 
    shaking of the ground, called seismic motion.  At higher frequencies, 
    there are lots of noise souces from vibrating instrument parts, like 
    shaking mirrors and tables.  
    """)

    ###
    # -- Show red noise with signal
    ###

    st.markdown("In the time-domain, you can see the red noise looks random.")

    with lock:
        figrnt = maze.plot(ylabel='Pressure')
        st.pyplot(figrnt)

    st.markdown("In the frequency-domain, the red noise has lots of power at low frequencies.")

    with lock:
        figrn = maze.asd(fftlength=1).plot(ylabel='Amplitude Spectral Density', ylim=[1e-11, 1e-4], xlim=[30, fs/2])
        st.pyplot(figrn)
        
    st.audio(make_audio_file(maze), format='audio/wav')
    st.markdown("""
    Can you hear the bullfrogs cheering?

    :point_right: **How does this compare with the white noise sound?**
    """)

if page == 4:

    # ----
    # Try to recover the signal
    # ----
    # st.markdown("## 4: Find the Secret Sound")
    
    st.markdown("""
    The red noise above isn't just noise - there's a secret sound 
    inside.  Did you hear it?  Probably not!  All of that low-frequency
    noise is making the secret sound very hard to hear.  But ... if the
    secret sound is at higher frequencies, maybe we could still hear it.
    
    What we need is a way to get rid of some of the low frequency noise, 
    while keeping the high frequency part of the signal.  In signal processing,
    this is known as a **high pass filter** - a filter that removes
    low frequency sounds, and keeps (or allows to *pass*) the high frequency 
    sounds.  The term **cutoff frequency** marks the boundary: frequencies
    below the cuttoff frequency are removed, and frequencies above the 
    cutoff frequency are passed.

    See if you can use a high pass filter to find the secret sound.  

    :point_right: **Adjust the 
    cutoff frequency using the slider below, and see if you can remove 
    some noise to find the secret sound.**

    """)

    lowfreq = st.slider("High pass filter cutoff frequency (Hz)", 0, 3000, 0, step=100)
    if lowfreq == 0: lowfreq=1

    highpass = maze.highpass(lowfreq)
    #st.pyplot(highpass.plot())

    with lock:
        fighp = highpass.asd(fftlength=1).plot(ylabel='Amplitude Spectral Density',
                                           ylim=[1e-12, 1e-5],
                                           xlim=[30, fs/2]
                                           )
        ax = fighp.gca()
        ax.axvspan(1, lowfreq, color='red', alpha=0.3, label='Removed by filter')
        st.pyplot(fighp)

    st.audio(make_audio_file(highpass), format='audio/wav')

    st.markdown("Can you hear the sound now?  What value of the cutoff frequency makes it easiest to hear?")

    st.markdown("")
    needhint = st.checkbox("Need a hint?", value=False)

    if needhint:

        st.markdown("""Here is the secret sound.  Can you find it hidden in the
        red noise above?
        """)

        st.audio(make_audio_file(secret), format='audio/wav')

        st.markdown("""You can also make the sound easier to hear by 
        clicking the 'Louder' option in the menu at left
        """)
        
if page == 5:
    # st.markdown("## 5: Whitening")

    st.markdown("""
    **Whitening** is a process that re-weights a signal, so that all
    frequency bins have a nearly equal amount of noise.  In our example,
    it is hard to hear the signal, because all of the low-frequency 
    noise covers it up.  By whitening the data,
    we can prevent the low-frequency noise from dominating what we hear. 
    
    :point_right: **Use the checkbox to whiten the data**
    """)

    
    whiten = st.checkbox("Whiten the data?", value=False)

    if whiten:
        whitemaze = maze.whiten()
    else:
        whitemaze = maze

    st.markdown("""
    After whitening, you can see the secret sound in the time domain.  You 
    may also notice that the whitened signal gently fades in at the beginning,
    and 
    out at the end - this gentle turn on / turn off is due to 
    **windowning**, and is important to apply in many signal processing
    applications.
    """)
    
    st.pyplot(whitemaze.plot())

    with lock:
        figwh = whitemaze.asd(fftlength=1).plot(ylim=[1e-12, 1], xlim=[30,fs/2], ylabel='Amplitude Spectral Density')
        st.pyplot(figwh)
    
    st.audio(make_audio_file(whitemaze), format='audio/wav')

    st.markdown("""Try using the checkbox to whiten the data.  Is it 
    easier to hear the secret sound with or without whitening?
    """)


if page == 6:

    # st.markdown("## 6: Gravitational Wave Data")

    st.markdown("""
    Finally, we'll try what we've learned on some real 
    gravitational-wave data from LIGO, around the binary black 
    hole signal GW150914.  We'll add one more element: 
    a **bandpass filter**.  A bandpass filter uses both a low frequency
    cutoff and a high frequency cutoff, and only passes signals in the 
    frequency band between these values. 

    :point_right: **Try using a whitening filter and a band-pass filter to reveal the
    gravitational wave signal in the data below.**  
    """)

    detector = 'H1'
    t0 = 1126259462.4   #-- GW150914

    st.text("Detector: {0}".format(detector))
    st.text("Time: {0} (GW150914)".format(t0))
    strain = load_gw(t0, detector)
    center = int(t0)
    strain = strain.crop(center-14, center+14)

    # -- Try whitened and band-passed plot
    # -- Whiten and bandpass data
    st.subheader('Whitened and Bandbassed Data')

    lowfreqreal, highfreqreal = st.slider("Band-pass filter cutoff (Hz)",
                                          1, 2000, value=(1,2000) )

    makewhite = st.checkbox("Apply whitening", value=False)

    if makewhite:
        white_data = strain.whiten()
    else:
        white_data = strain

    bp_data = white_data.bandpass(lowfreqreal, highfreqreal)

    st.markdown("""
    With the right filtering, you might be able to see the signal in the time domain plot.
    """)

    with lock:
        fig3 = bp_data.plot(xlim=[t0-0.1, t0+0.1])
        st.pyplot(fig3)

    # -- PSD of whitened data
    # -- Plot psd
    with lock:
        psdfig = bp_data.asd(fftlength=4).plot(xlim=[10, 1800], ylabel='Amplitude Spectral Density')    
        ax = psdfig.gca()
        ax.axvspan(1, lowfreqreal, color='red', alpha=0.3, label='Removed by filter')
        ax.axvspan(highfreqreal, 1800, color='red', alpha=0.3, label='Removed by filter')
        st.pyplot(psdfig)

    # -- Audio
    st.audio(make_audio_file(bp_data.crop(t0-1, t0+1)), format='audio/wav')

    # -- Close all open figures
    plt.close('all')

    st.markdown("""With the right filtering, you might be able to hear
    the black hole signal.  It doesn't sound like much - just a quick thump.  
 """)

    st.markdown("")
    hint = st.checkbox('Need a hint?')

    if hint:

        st.markdown("""
        Hint: Try using a band pass from 30 to 400 Hz, with whitening on.
        This is similar to what was used for Figure 1 of the 
        [GW150914 discovery paper](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.116.061102), also shown below:
        """)
        
        st.image('https://journals.aps.org/prl/article/10.1103/PhysRevLett.116.061102/figures/1/large')



st.markdown("""## About this app

This app displays data from LIGO, Virgo, and GEO downloaded from the Gravitational Wave Open Science Center at 
[https://gw-openscience.org](https://gw-osc.org).
""")
