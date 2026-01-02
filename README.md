# moire-exciton-solver
The Continuum model is the standard theoretical framework for describing
Moiré systems with small ($`\theta < 5^\circ`$) twist angle. Instead of
simulating millions of atoms such as in a hubbard model simulation, we
approximate the atomic lattice potential into a continuous field, and
can select for desirable particles to perform additional operations on
in this field. Definitions 2 and 3 are the building blocks of this
model.

<div class="definition">

**Definition 1**. Excitons are quasiparticles that are electron-hole
pairs bonded by coulomb forces. These quasiparticles define the physics
of described device since coulomb forces are magnified by the quantum
confinement seen in 2D materials. Intuitively, this is due to
confinement forcing carriers closer together and thus increasing the
magnitude of charge between them. More of the hole and electron
wavefunctions overlap. These can be of multiple types but only type II
will be investigated.  
In XSe2-XS2 bilayers, a type II band gap is established, due to the
relationship
``` math
E_{C,XS_2}<E_{C,XSe_2},E_{V,XS_2}<E_{V,XSe_2}
```
This leads to the valence band being hosted by the XS2 layer, and the
conduction band by the XSe2 layer, forming the bandgap
``` math
{E_g=E}_{V,XS_2 }-E_{C,XSe_2}
```
The resulting interlayer charge separation leads to the formation of
interlayer excitons.  

</div>

<div class="definition">

**Definition 2**. The moiré superlattice is determined by the
interference between the two hexagonal lattices twisted by angle
$`\theta`$ and lattice mismatch $`\delta = \frac{|a_1-a_2|}{a_{avg}}`$.
The interference pattern is sinusoidal, and the superlattice period is
defined by the equation
``` math
\begin{equation}
    L_m = \frac{(1+\delta)a_{min}}{\sqrt{2(1+\delta)(1-cos\theta) + \delta^2}}
\end{equation}
```
This can be simplified down for our purposes, with small lattice
mismatch and twist angle to
``` math
\begin{equation}
    L_m \approx \frac{a_{min}}{\sqrt{\delta^2 + \theta^2}}
\end{equation}
```
The potential landscape over this can be derived from the local atomic
registry. It matches the $`L_m`$ value and is approximated by the first
harmonic expansion of the Moiré reciprocal lattice vectors $`G_j`$.
``` math
\begin{equation}
    V(r) = 2V_0 \sum^3_{j=1}cos(G_j.r + \psi)
\end{equation}
```
where $`V_0`$ is the potential amplitude and $`\psi`$ is the stacking
phase parameter.

</div>

<div class="corollary">

**Corollary 1**. *Holes have a higher effective mass than electrons,
which leads to the physics landscape being dominated by hole physics. As
a result, we can rely on the hole wavefunction to get satisfactorily
accurate measurements for our system. This can be more rigorously
derived from the bandwidths of the valence and conduction bands. The
bandwidth of a given band is given by
``` math
\begin{equation}
        W \propto \frac{\hbar^2}{m^* L_m^2}
\end{equation}
```
Since hole mass is greater than electron mass, the hole band forms the
localized trap needed for a qubit while the electron stays more
delocalized.  
At small angles ($`\theta \approx 1^\circ`$) in $`WS_2/WSe_2`$ bilayers,
the hole bandwidth W becomes extrememly narrow. This is the flat band
regime, where holes stop moving and get localized in the moiré potential
wells, forming an ordered Moiré-Wigner crystal. These are the quantum
dots alluded to previously.  
Additionally, type-II excitons such as in this system are more stable,
since the electron and hole are localized in different layers. This is
due to the fact that the energy barrier for recombination is greater
than in monolayer Moiré systems. This brings the decoherence time from a
few picoseconds to a relatively usable nanosecond timeframe.*

</div>

<div class="definition">

**Definition 3**. For a charge carrier in layer $`i`$, the effective
mass Hamiltonian is:
``` math
\begin{equation}
H = -\frac{\hbar^2}{2m^*} \nabla^2 + V(\mathbf{r})
\end{equation}
```

Kinetic Energy (Finite Difference Method): We discretize the
wavefunction $`\psi(x,y)`$ on a grid with spacing $`\Delta x`$. The
Laplacian $`\nabla^2`$ is approximated by a 5-point stencil:
``` math
\begin{equation}
\nabla^2 \psi_{i,j} \approx \frac{\psi_{i+1,j} + \psi_{i-1,j} + \psi_{i,j+1} + \psi_{i,j-1} - 4\psi_{i,j}}{\Delta x^2}
\end{equation}
```
Substituting this into the kinetic energy term
$`T = -\frac{\hbar^2}{2m^*} \nabla^2`$, we get hopping terms $`t`$ and
an on-site energy $`\epsilon`$:
``` math
\begin{equation}
t = -\frac{\hbar^2}{2m^* \Delta x^2}, \quad \epsilon = -4t
\end{equation}
```

</div>

<div class="definition">

**Definition 4**. Qubit Control (Rabi Oscillations): When driving the
qubit with a classical laser field $`\Omega \cos(\omega_L t)`$, the
hamiltonian in the Rotating Wave Approximation is:
``` math
\begin{equation}
H_{RWA} = \frac{\hbar \Delta}{2} \sigma_z + \frac{\hbar \Omega}{2} \sigma_x
\end{equation}
```
This allows for standard universal quantum gates (X, Z, H) also known as
the Pauli gate set.  
$`\Delta`$ is the detuning coefficient, and describes the frequency
mismatch of the laser. If it is off-resonance, the qubit vector rotates
around the Z axis and accumulates phase.  
$`\Omega \propto d.\epsilon_0`$ is the Rabi frequency proportional to
laser power. This rotates the qubit vector around the X axis and
transitions between the $`|0\rangle`$ and $`|1\rangle`$ states.  
Careful selection of these parameters can allow for the implementation
of the clifford group of quantum gates but that is left to the user.

</div>

<div class="definition">

**Definition 5**. Quantum Optics (Jaynes-Cummings): To simulate single
photon emissions, the system is coupled to a quantized cavity. The state
occupies the Hilbert space
$`\mathcal{H} = \mathcal{H}_{atom} \otimes \mathcal{H}_{photon}`$ and
the time dependent hamiltonian is given by.
``` math
\begin{equation}
        H_{JC} = \hbar \omega_c a^\dagger a + \frac{\hbar \omega_q}{2} \sigma_z + \hbar g (a^\dagger \sigma_- + a \sigma_+)
\end{equation}
```
where $`a^\dagger a`$ is the number of photons in the cavity,
$`a^\dagger(a) \sigma_{+(-)}`$ describes the emission (absorption)
process. (so for no. of photons n, $`n \rightarrow n +(-) 1`$). The
vacuum Rabi coupling factor g describes the speed of this energy
exchange.

</div>

<div class="definition">

**Definition 6**. Decoherence and Decay (Lindblad Master Equation): The
Gorini-Kossakowski-Sudarshan-Lindblad master equation is the most
general generator of Markovian dynamics in quantum systems, and thus is
the system of choice for modeling decay and decoherence in our system.
For a density matrix $`\rho`$ we have
``` math
\begin{equation}
        \frac{d{\rho}}{dt} = \frac{-i}{\hbar}[H, \rho] + \kappa \mathcal{D}[a]\rho + \gamma \mathcal{D}[\sigma_-]\rho
\end{equation}
```
where the dissipator is described by
``` math
\begin{equation}
         D[L]\rho = L\rho L^\dagger - \frac{1}{2}\{L^\dagger L, \rho\}
\end{equation}
```
$`\kappa`$ is cavity leakage, describing the amount of photon leakage
out of the cavity. This is desirable for single-photon sources since it
is the output signal.  
$`\gamma`$ is non-radioactive decay, describing the rate at which the
exciton decomposes into phonons.

</div>

<div class="definition">

**Definition 7**. The Purcell effect is a key concept in quantum
photonics, due to its utility in controlling the rate of spontaneous
emission. This helps enhance device efficiency and makes single photon
sources possible.  
Given that a photon decays at a rate of $`\Gamma_0`$ in free space, the
expression changes inside a cavity due to confinement modifying the
density of states available to the photon. If the cavity is resonant
($`\omega_q = \omega_e`$) the decay is enhanced by the Purcell Factor
$`F_p`$.
``` math
\begin{equation}
        \Gamma_{cav} = F_p \Gamma_0 = \frac{3}{4\pi^2} (\frac{\lambda}{n})^3 \frac{Q}{V_{mode}}\Gamma_0
\end{equation}
```
This enhancement improves the odds of the photon being emitted into the
cavity before it is lost to non-radioactive decay. Furthermore, it leads
to a brighter emission since less energy is lost before a photon is
emitted.
