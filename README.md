# Interactive Moiré Exciton Simulation Platform

### Overview
Atomically thin [transition metal dichalcogenides (TMDcs)](#background) have gained merit over the last two decades since the discovery of direct bandgaps in 2010. In 2D materials, electrons and holes are confined and form quasiparticles called [excitons](#def-exciton). These can be used as a quantum computation platform.

With lattice mismatch, a periodic potential landscape forms, known as a [moiré potential landscape](#def-moire). This potential landscape acts as an array of quantum dots, trapping interlayer excitons. This is tunable via twist angle, and in some materials (unimplemented) with applied potential.

This repo implements a toy moiré qubit by calculating material properties using a [Continuum Model](#def-continuum) to determine bandstructure, simulating Quantum Control via [Rabi oscillation](#def-rabi) implementations of Quantum gates and [Decoherence](#def-decoherence), and modelling Quantum Optics via the [Jaynes-Cummings Model](#def-jc) to simulate Cavity QED effects like the [Purcell Effect](#def-purcell) for single-photon emissions.

This Project features a Python Backend for the physics calculations and React Frontend for real-time interactive visualization.

# Background

Atomically thin transition metal dichalcogenides (TMDcs) have gained
merit over the last two decades since the discovery of direct bandgaps
in 2010. In 2D materials, electrons and holes are confined to just two
dimensions, increasing the coulomb attraction experienced. This causes
electrons and holes to be bound to each other, creating excitons. In
these materials, exciton binding energies $`E_b`$ are high, in the order
of 500 meV, as opposed to a typical bulk semiconductor material’s
$`E_b`$ of 10. This is high enough for them to be stable at room
temperature, since the thermal energy $`k_bT`$ at room temperature is
only 25.8meV. Yet, $`E_b`$ is low enough that excitonic effects dominate
all transport and optical properties in TMDcs.  
  
One promising family of materials is the group VI family of MX2, where M
is Mo or W, and X is S or Se. These materials are stable at room
temperature and have unique optical properties that can be engineered
with. They have an indirect band gap in bulk, but in atomically thin
layers, their band extrema are located at the K+ and K- locations on the
Brillouin zone, forming a direct bandgap.  
  
Furthermore, since the discovery of superconductivity in twisted bilayer
Graphene, twisted TMD’s have been in vogue, with some angles forming
periodic potential landscapes named Moiré superlattices. These
potentials act as quantum dots, trapping interlayer excitons. These dots
are tunable via their twist angle, making them promising for quantum
information theoretic and quantum optoelectronic applications.  
  
For the aforementioned reasons, the simulation of their physics is
attractive. A more accurate model that can accurately describe moiré
excitonic physics is of utility to researchers who wish to hone in on
optimal device parameters without wasting expensive $`WS_2/WSe_2`$
wafers. This area is active and is in a transitory period, which will be
touched on in a later section.  
  
Finally, there is a substantial lack of tools that help visualize moiré
excitons and the physics associated with them. No comparable tool to the
one outlined in this paper is known to the author at the time of
writing. A student looking to gain an understanding would have to
painstakingly work through hundreds of pages worth of dense research
papers to gain an intuition for which knobs to turn to induce desired
properties in the system. Any such reader is implored to set up the tool
and play around with it.


# Implemented Physics

This section lays out the base physics required to build our model. For
clarity, we focus on the $`WS_2/WSe_2`$ heterobilayer, but model
generality allows for alternate materials given they are also TMD
heterobilayers.  
The Continuum model is the standard theoretical framework for describing
Moiré systems with small ($`\theta < 5^\circ`$) twist angle. Instead of
simulating millions of atoms such as in a hubbard model simulation, we
approximate the atomic lattice potential into a continuous field, and
can select for desirable particles to perform additional operations on
in this field. Definitions 2 and 3 are the building blocks of this
model.

<div class="definition">
<a id="def-exciton"></a>
  
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
\begin{equation}
E_{C,XS_2} < E_{C,XSe_2},E_{V,XS_2} < E_{V,XSe_2}
\end{equation}
```
This leads to the valence band being hosted by the XS2 layer, and the
conduction band by the XSe2 layer, forming the bandgap
``` math
\begin{equation}
{E_g=E}_{V,XS_2 }-E_{C,XSe_2}
\end{equation}
```
The resulting interlayer charge separation leads to the formation of
interlayer excitons.  

</div>

<div class="definition">
<a id="def-moire"></a>
  
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
<a id="def-continuum"></a>
  
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
<a id="def-rabi"></a>
  
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
<a id="def-jc"></a>
  
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
<a id="def-decoherence"></a>
  
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
<a id="def-purcell"></a>
  
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

</div>

# Computational Methods

The physics simulation is written in python with a flask frontend. Each
of the above definitions are built into functions that a user may
interact with directly. Some considerations have been made to make
simulation computationally viable and the frontend real-time. The grid
and bandstructure generators are $`O(n^3)`$, with the two body
hamiltonian used for the qubit being $`O(n^8)`$ where n is grid size.
The former is computable while the latter is infeasible even for small
grids. A 20x20 grid would require terabytes of space to hold, and the
simulation platform targets 10 points per nm, which leads to a 78x78
grid at $`\theta = 0.5^\circ`$. This was alleviated by parallelizing the
hamiltonian calculation over multiple cores, and iteratively solving for
the eigenvalues of the hamiltonian.  
Qubit visualization was handled with custom bloch sphere object in react and
other data used basic data visualization routines.

# Results

Due to the simple model used and placeholder material parameters, a
quantitative analysis of experimental data would be unfruitful. DFT
parameter extraction is required for this to be possible, something that
may be done at a later date.  
Qualitatively viable flat-bands are derived with gap between the first
two states below 100meV, with visible dissipative effects as $`\theta`$
increases.  
Physically accurate potential landscape and real lattice are generated.
The behavior is as expected with changes in $`\theta`$.

# Conclusion and State of Current Research

This paper outlines a satisfactorily accurate simulation platform for
moiré excitons, with an interactive user interface that lends itself
well to education.  
Experimental research on this topic is thriving but slow despite being
well funded due to the high degree of difficulty associated with
maintaining consistent parameters across supercells. Recent research is
predominantly on graphene based Moiré structures, and most papers
concern themselves with mapping physics over building viable devices.  
The simulation landscape is more vibrant, with a recent shift away from
Hubbard models, which is the standard model used for many body physics.
The model outlined in this paper generalizes to a Hubbard model at small
twist angles. It was seen that standard Bosonic Bose-Hubbard models
predicted a greater number of interparticle interactions, while real
moiré excitons were never seen to have over four interactions in TMD’s.
To correct this, fermionic models were proposed. The first few models
are in peer review at the time of writing.


# References

1. G. Wang et al., “Colloquium : Excitons in atomically thin transition
metal dichalcogenides,” Reviews of Modern Physics, vol. 90, no. 2, Apr.
2018, doi: https://doi.org/10.1103/revmodphys.90.021001  
  
2. K. F. Mak and J. Shan, “Photonics and optoelectronics of 2D
semiconductor transition metal dichalcogenides,” Nature Photonics, vol.
10, no. 4, pp. 216–226, Mar. 2016, doi: https://doi.org/10.1038/nphoton.2015.282  
  
3. H. Baek et al., “Highly energy-tunable quantum light from
moiré-trapped excitons,” Science advances, vol. 6, no. 37, Sep. 2020, doi: https://doi.org/10.1126/sciadv.aba8526  
  
4. H. Guo, X. Zhang, and G. Lu, “Shedding light on moiré excitons: A
first-principles perspective,” Science Advances, vol. 6, no. 42, Oct.
2020, doi: https://doi.org/10.1126/sciadv.abc5638.  
  
5. T.-S. Huang, P. Lunts, and M. Hafezi, “Nonbosonic Moiré Excitons,”
Physical Review Letters, vol. 132, no. 18, pp. 186202–186202, Apr.
2024, doi: https://doi.org/10.1103/physrevlett.132.186202.  
  
6. L. Yuan et al., “Twist-angle-dependent interlayer exciton diffusion
in WS2–WSe2 heterobilayers,” Nature Materials, vol. 19, no. 6, pp.
617–623, May 2020, doi: https://doi.org/10.1038/s41563-020-0670-3.  
  
7. Y. Liu et al., ““Ideal” Topological Heavy Fermion Model in
Two-dimensional Moiré Heterostructures with Type-II Band Alignment”,
preprint, July 2025, arxiv: https://arxiv.org/pdf/2507.06168  
  
8. Y. T. Wang et al., “Moiré cavity quantum electrodynamics,” Science
Advances, vol. 11, no. 21, May 2025, doi: https://doi.org/10.1126/sciadv.adv8115.  
  
9. H. Wang et al., “Quantum coherence and interference of a single moiré
exciton in nano-fabricated twisted monolayer semiconductor
heterobilayers,” Nature Communications, vol. 15, no. 1, pp. 4905–4905,
Jun. 2024, doi: https://doi.org/10.1038/s41467-024-48623-4.  
  
10. M. Xie, M. Hafezi, and S. Das Sarma, “Long-Lived Topological
Flatband Excitons in Semiconductor Moiré Heterostructures: A Bosonic
Kane-Mele Model Platform,” Physical Review Letters, vol. 133, no. 13,
Sep. 2024, doi: https://doi.org/10.1103/physrevlett.133.136403.
