# SigmaTransform-Cpp
Implements the Sigma Transform in C++

This repository shows an exemplary implementation of the "SigmaTransform", as defined in the thesis _"Quantum Frames and Uncertainty Principles arising from Symplectomorphisms"_, written in C/C++.

The code was compiled and tested with g++ (GCC) 4.8.1, on Windows 7 and g++ (Debian 4.7.2-5) 4.7.2 on Debian Linux and uses the C++11 standard, as well as the FFTW library (http://fftw.org/), which should be installed prior to compilation. For Windows, a current shared-library version of FFTW3 (libfftw3-3.dll) is provided in the subdirectory ./FFTW. On a Linux machine, the library should be installed 
via a package manager, e.g. using

    $ sudo apt-get install fftw3 fftw3-dev
    $ sudo pacman -S fftw3 fftw3-dev
    $ sudo yum install fftw3 fftw-devel
    $ sudo rpm -i fftw3 fftw3-devel
    
depending on your distribution, or compiled from scratch using the *--enable-threads* flag.

The examples

    Example1D_ConstantQ.cpp     # The 1D Constant-Q Transform
    Example1D_STFT.cpp          # The 1D Short-Time Fourier Transform
    Example1D_Wavelet.cpp       # The 1D Wavelet Transform
    Example1D_async.cpp         # Using the implementation asynchronously
    Example1D_inline.cpp        # Using the implementation inline
    Example1D_threads.cpp       # Using multiple threads/parallel processing
    Example2D_Curvelet.cpp      # The 2D Curvelet Transform
    Example2D_NPShearlet.cpp    # The Non-Parabolic Shearlet Transform
    Example2D_SIM2.cpp          # The SIM(2)-Transform
    Example2D_STFT.cpp          # The 2D Short-Time Fourier Transform
    Example2D_Wavelet.cpp,      # The 2D Wavelet Transform

located in the *./Examples* subdirectory show how to use the implementation, along with some special cases. The provided makefile should compile and link all examples - on Windows as well as Linux with the appropriate tools and libraries installed -, as well as the Code for the SigmaTransform itself. The binaries will be put into the subdirectory ./bin.

Note that this library is not intended to show maximal performance, but show the usability of the universal interface of the "Sigma Transform" to perform well-known signal processing transforms / algorithms, differing only by single paramater - the "spectral diffeomorphism".
