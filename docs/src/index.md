# ASSET Documentation

Slit spectroscopy acquisitions are two-dimensional data, consisting of one spatial and one spectral dimension. The most classical way to extract spectrum from spectroscopic data is the aperture extraction. In the context of slit data, it consists in integrating the photons, for each wavelength value in the spectral dimension, in a window of a given height. This method is suboptimal as it does not take into account the full width at half maximum (FWHM) of the instruments point spread function (PSF) or the noise statistic. Moreover, it can integrate artifacts (bad pixels, cosmic ray) that can pollute the spectral extraction.

Inverse problems framework is widely used in astrophysics to tackle these
issues, such as @Michalewicz:2023 and @Berdeu:2024.

The method used in the [`ASSET`](https://github.com/SlitSpectroscopyBuddies/ASSET) package was proposed by @The:2023 for the ESO/VLT SPHERE LSS, then adapted in @Denneulin:2023 for the JWST/NIRSpec FS. It requires the knowledge of the following maps as inputs: 

- the L data maps $(d_\ell)_{\ell \in {1:L}}$
- the L weights maps $w_\ell$, which each element can be computed as the inverse variance of the pixel, forming the matrix $W_{\ell}= \Sigma_{\ell}^{-2}$. We assume that a defective pixel or artifacts have an infinite variance, i.e. a zero entry in $W_{\ell}$.
- the L spatial distribution maps $X_\ell$ where $0$ should correspond to the center of the studied object. #TODO: implement it, it should not be too long even if it seems a big deal
- the L spectral distribution maps $\Lambda_\ell$.

Such maps must be stored in a `CalibratedData` structure using the constructor `CalibratedData(d, w, ρ_map, λ_map)`.

The outputs of the `extract_spectrum!` method are the extracted spectrum $z$, sampled over a given regular wavelength grid $(\lambda_n)_{n \in 1:N}$, and the parameters $\theta$ of the fitted PSF model, stored in a its corresponding structure. The spatial distribution maps $X$, stored in the  `CalibratedData` structure,  and the background map $b$, stored in a `BkgMdl` structure, are auto-calibrated in place.  

`ASSET` uses an inverse problem approach to extract the spectra by fitting a chromatic PSF model on the data.  
$$
z,\theta \in \mathrm{arg min} \Big\{\sum_\ell\Vert d_\ell - (m_\ell(z,\theta) + b)\Vert_{W_{\ell}}^2 + \mu_z\mathcal{R}_z(z) + \mu_\theta\mathcal{R}_\theta(\theta) + \mu_b\mathcal{R}_b(b), \Big\}
$$
where $\mathcal{R}_z$, $\mathcal{R}_\theta$ and $\mathcal{R}_b$ are respectively the regularization of the extracted spectra, of the PSF parameters if required and of the background, with hyperparameters $\mu_z$, $\mu_\theta$ and $\mu_b$. 
$$
m_\ell(z,\theta) = \alpha_\ell Z(\Lambda_\ell,z) \odot H(\theta, X_\ell, \Lambda_\ell)
$$
is the model of the data. It is the Hadamard (element-wise) product of the spectrum interpolated in the camera plane
$$
Z(\Lambda,z)_{j,k}=\sum_n\phi\Big(\frac{\Lambda_{j,k}-\lambda_n}{\delta_\lambda}\Big)z_n
$$

with $\phi$ an interpolation kernel, and of the chromatic PSF $H$. This PSF can have different degrees of freedom depending of the chosen model. In `ASSET` framework, we distinguish  `ParametricPSF` , with low degree of freedom, and `NonParametricPSF` , with hight degree of freedom.  A `ParametricPSF`  $H$ is a function parametrized by a few unknown variables $\theta_n$, *e.g.* a Gaussian chromatic with a minimum width model:
$$
H(\theta, X_\ell, \Lambda_\ell)_{j,k} = \big(2 \pi(\theta_1 \Lambda_{j,k\ell}^2 + \theta_2)\big)^{-1}\exp\Big(-\frac{X_{j,k,\ell}^2}{2(\theta_1\Lambda_{j,k,\ell}^2+\theta_2)}\Big).
$$

It does not require any regularization. See  @The:2023 and @Denneulin:2023 for more details. 

A `NonParametricPSF` is parametrized directly by some profile $(\theta_m)_{m \in  1:M}$,  *e.g.* taking $o$ order of the speckles expansion model [@Devaney:2017] which is the interpolation of the profiles $\theta_o$ in a reference plane of the spatial coordinates $(x_m)_{m \in 1:M}$,:
$$
H(\theta, X_\ell, \Lambda_\ell)_{j,k} = \sum_o \gamma(\Lambda_\ell)^o \sum_m \phi\Big(\frac{\gamma(\Lambda_\ell)X_\ell - x_m}{\delta x}\Big)\theta_{m,o}(x_m).
$$

For such a PSF, $\theta$ must be regularized. The regularization type and the initial hyperparameter must be given in the `NonParametricPSF` structure. The hyperparameter is  auto-calibrated in the method, see @The:2023 for and references therein for more details.

The package provide several `ParametricPSF` and `NonParametricPSF` and the users can easily implement their own. 


```@contents
```

```@docs
ChromaticSeriesExpansionsInterpolator
CalibratedData
```

```@index
```
