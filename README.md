# T-wave_end_ECG
The algorithm can identify the end of T-wave in an ECG signal.

Files:
- 'sel33.dat' contains the ECG signal;
- 'sel33.hea' explain the content of the signal;
- 'sel33.q1c' contains manual annotations related to a segment of signal (about the start and end of P, T, QRS waves). T-peaks are identified as 't'.

Pre processing:
- high-pass filtering (2 Hz)
- low-pass filtering (30 Hz)

Elaboration:
- Find 'xm', the minimun of the derivative of the signal after the maximum of a positive T-wave, searched in a window of 200 ms after T-peak.
- Find 'xr', a point in the range between 200 and 400 ms after T peak, in the isoelectric segment, where the derivative is under a certein thresholde (near zero). If there is no a point that satisfy this requirement, choose the central point in this range.
- Compute the area of a trapezoid, where the major area is from the negative peak of the derivative of T-wave (xm) and a far point in the isoelectric (xr), while the minor area is from the same far point (xr) and a point in the descent front of P-wave (xend)
- The final point of T-wave (xend) is find as the point that maximize the area of the trapezoid.

       xm,y(xm)       xend,y(xm)                xr,y(xm)
         .________________.________________________. 
          \                                        |
           \                                       |
             \                                     |
               \                                   |
                  \                                |
                     \                             |
                        \.________________________.|
                     xend,y(xend)                    xr,y(xend)
