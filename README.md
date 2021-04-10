## Analysis of spin-wave dispersion in magnonic structures
A cluster of .ovf (OOMMF vector format) files containing 
spatio-temporal information of spin-waves are post-processed
to plot the spin-wave dispersion curve. The post-processing 
algorithm is adopted from [[1]](#1).

Additionally, space-time variation of spin-waves can be also
be plotted and theoretical calculation of dispersion function
and group velocity of the lowest order mode of the three classes of
spin-waves namely Backward volume spin-waves (BVSWs), Forward
volume spin-waves (FVSWs) and Surface spin-waves (SSWs) are done with and without the exchange effcet.

## Units of the parameters
| Parameters        | Units|
| ------------- |-------------|
| time      | ns|
| length      | nm|
| magnetic field | kA/m |
| frequency | GHz|
| angular frequency| Grad/s|
| wave-number| rad/nm|
| group velocity| km/s|
| frequency normalized by <img src="https://render.githubusercontent.com/render/math?math=f_{\text{M}} = \gamma_0 M_{\text{S}}, \overline{f}"> | 1 |
| group velocity normalized by <img src="https://render.githubusercontent.com/render/math?math=f_{\text{M}} = \gamma_0 M_{\text{S}}, \overline{v_{\text{g}}}"> | nm/rad |
| <img src="https://render.githubusercontent.com/render/math?math=\gamma_0 = \frac{\gamma \mu_0}{2 \pi}">| GHz.m/kA|

## References
<a id="1">[1]</a> 
Kumar, D., Dmytriiev, O., Ponraj, S. and Barman, A., 2011. Numerical calculation of spin wave dispersions in magnetic nanostructures. Journal of Physics D: Applied Physics, 45(1), p.015001.
