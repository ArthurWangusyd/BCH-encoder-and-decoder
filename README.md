# BCH encoder and decoder (QPSK)
This project implements BCH encoding and decoding without using MATLAB's built-in equations. 
The modulation method is QPSK.
The project uses two methods for the decoding part. The first is in C and the other method is implemented entirely in matlab. The decoder compiled in C has better performance. 
And the M1&2&3 files show three ways to use it, comparing pe VS BER, SNR VS BER, etc. The mathematical method used for the decoder is the Bredekamp's iterative method.
