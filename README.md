# Kr pickup in helium nanodroplets

A small Python script for calculating the **mean number of captured Kr
atoms** using either a fixed-cross-section Poisson model or a Kuma-type
analytical droplet-shrinkage model.

## Files

- `Kuma_Poisson.py`: standalone calculation and plot.
- `requirements.txt`: Python dependencies.
- `.gitignore`: excludes local environments, caches and IDE settings.
- `references.bib`: reference for the analytical shrinkage model.

## Installation and execution

Use Python 3.10 or newer. In the project directory, run:

```bash
python -m pip install -r requirements.txt
python Kuma_Poisson.py
```

Alternatively, open and run `Kuma_Poisson.py` in PyCharm using an interpreter
with NumPy and Matplotlib installed.

## Settings

Edit the settings at the top of the script:

```python
MODEL = "kuma"       # or "poisson"
P_input = 1.3e-3     # pressure for the printed result, in mbar
P_max = 3.0e-3       # maximum plotting pressure, in mbar
```

Input pressures are **local Kr pressures in the pickup cell**, already
calibrated and background-subtracted. They are not raw chamber-gauge readings.

The script prints the mean pickup and intermediate physical quantities,
then opens a plot. It does not sample individual pickup events.

## Models and assumptions

For pressure `P` in Pa, the Poisson mean is `X * P`.
The Kuma-type mean is:

```text
mean_pickup = (N/c) * [1 - (1 - c*X*P/(3*N))^3]
```

Here `N` is the initial He-atom number and `c` is the assumed constant number
of He atoms evaporated per captured Kr atom. The droplet is spherical, with
a fixed He number density and geometrical cross section.

The relative-speed prescription is:

```text
v_dop = sqrt(3*kB*T/m)
v_rel = sqrt(v_drop^2 + v_dop^2)
```

This uses the RMS thermal speed of the gas and is a prescribed approximation,
not an exact average of the relative speed over the thermal distribution.

The energy deposited per captured atom is approximated by the sum of
`3*kB*T/2 + m*v_drop^2/2`, He-Kr binding energy, and bulk-like Kr aggregation
energy. Dividing by the energy per evaporated He atom gives `c`.

The assumed binding energies are 217 and 932 cm^-1; the evaporation energy
is 5 cm^-1 per He atom. These are model inputs, not universal constants.
In particular, the constant bulk-like aggregation energy does not resolve
individual small-cluster growth steps.

For Kuma, the complete-evaporation pressure is `3*N/(c*X)` in Pa. The script
rejects input pressures above this limit and truncates the plotted curve
there. It does not describe further growth of a bare Kr cluster.
The Poisson model neglects shrinkage, so it does not impose this bound.

Default example parameters are `N = 9678`, `T = 300 K`, `v_drop = 345 m/s`,
`L = 18 mm`, and `rho_He = 2.18e28 atoms/m^3`.

This standalone script is **not HeDopant** and is not a verified reproduction
of a particular manuscript figure. It does not include a droplet-size
distribution, velocity evolution, pressure gradients, external-chamber
pickup, skimmer transmission, or detector weighting.

## Reference

S. Kuma, H. Goto, M. N. Slipchenko, A. F. Vilesov, A. Khramov, and T. Momose,
"Laser induced fluorescence of Mg-phthalocyanine in He droplets: Evidence
for fluxionality of large H2 clusters at 0.38 K,"
*J. Chem. Phys.* **127**, 214301 (2007).
https://doi.org/10.1063/1.2802303

This implementation uses an analytical shrinkage expression attributed to
that work with the relative-speed prescription and energy assumptions
specified above. It is not the authors' original software.
