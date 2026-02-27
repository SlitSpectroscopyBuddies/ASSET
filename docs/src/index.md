# ASSET Documentation
```@contents
Pages = ["documentation.md","estimation.md"]
```

```@docs
AbstractBkg
BkgMdl
extract_spectrum!
ChromaticSeriesExpansionsInterpolator
CalibratedData
```
- link to [`ASSET.AbstractBkg`](@ref)
- link to [`ASSET.BkgMdl`](@ref)

```@index
```

# Theoretical Background

The method used in the [`ASSET`](https://github.com/SlitSpectroscopyBuddies/ASSET) package requires the following maps as inputs: 

- data maps $(d_\ell)_{\ell \in {1:L}}$, where $L$ is the amount of dithers/acquisitions/frames;
- weights maps $w_\ell$, where each element can be computed as the inverse variance of the pixel, forming the matrix $W_{\ell}= \Sigma_{\ell}^{-2}$. We assume that a defective pixel or artifacts have an infinite variance, i.e. a zero entry in $W_{\ell}$;
- spatial coordinate maps $X_\ell$ where $0$ should correspond to the center of the studied object;
- spectral coordinate maps $\Lambda_\ell$.

The method outputs are the extracted spectrum $z$, sampled over a given regular wavelength grid $(\lambda_n)_{n \in 1:N}$, and the parameters $\theta$ of the fitted PSF model. They are obtained by solving:  

$
z,\theta \in \mathrm{arg min} \Big\{\sum_\ell\Vert d_\ell - (m_\ell(z,\theta) + b)\Vert_{W_{\ell}}^2 + \mu_z\mathcal{R}_z(z) + \mu_\theta\mathcal{R}_\theta(\theta) + \mu_b\mathcal{R}_b(b), \Big\}
$

where $\mathcal{R}_z$, $\mathcal{R}_\theta$ and $\mathcal{R}_b$ are respectively the regularization of the extracted spectra, of the PSF parameters if required and of the background, with hyperparameters $\mu_z$, $\mu_\theta$ and $\mu_b$. The spatial distribution maps $X$ and the background map $b$ are auto-calibrated in the process.  The model of the data

$
m_\ell(z,\theta) = \alpha_\ell Z(\Lambda_\ell,z) \odot H(\theta, X_\ell, \Lambda_\ell)$

is the Hadamard (element-wise) product of the spectrum interpolated in the camera plane

$
Z(\Lambda,z)_{j,k}=\sum_n\phi\Big(\frac{\Lambda_{j,k}-\lambda_n}{\delta_\lambda}\Big)z_n
$

with $\phi$ an interpolation kernel, and of the chromatic PSF $H$ (parametric or non-parametric). The package provide several `ParametricPSF` and `NonParametricPSF` and the users can easily implement their own. A `ParametricPSF`  $H$ is a function parametrized by a few unknown variables $\theta_m$, thus with low degrees of freedom, *e.g.* a Gaussian chromatic with a minimum width model:

$
H(\theta, X_\ell, \Lambda_\ell)_{j,k} = \big(2 \pi(\theta_1 \Lambda_{j,k\ell}^2 + \theta_2)\big)^{-1}\exp\Big(-\frac{X_{j,k,\ell}^2}{2(\theta_1\Lambda_{j,k,\ell}^2+\theta_2)}\Big).
$

It does not require any regularization. See @The:2023 and @Denneulin:2023 for more details. 

A `NonParametricPSF` is parametrized directly by some profile $(\theta_m)_{m \in  1:M}$, thus with high degrees of freedom , *e.g.* taking $o$ order of the speckles expansion model [@Devaney:2017] which is the interpolation of the profiles $\theta_o$ in a reference plane of the spatial coordinates $(x_m)_{m \in 1:M}$,:

$
H(\theta, X_\ell, \Lambda_\ell)_{j,k} = \sum_o \gamma(\Lambda_\ell)^o \sum_m \phi\Big(\frac{\gamma(\Lambda_\ell)X_\ell - x_m}{\delta x}\Big)\theta_{m,o}(x_m).
$

For such a PSF, $\theta$ must be regularized. The hyperparameter is auto-calibrated in the method (see @The:2023 for and references therein for more details). 
