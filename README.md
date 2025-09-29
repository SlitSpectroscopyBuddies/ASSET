# Adaptable Slit Spectral Extraction Tool (ASSET)

This package aims at giving low-level general tools to the extraction of spectrum of an unresolved object with slit spectroscopy.

## Installation

Many of the dependencies are registered in the `EmmtRegistry`. In the Julia package manager, simply run:

``````julia
pkg> registry add https://github.com/emmt/EmmtRegistry
``````

To install the necessary `InverseProblem` package, run:

```julia
pkg> add https://github.com/SJJThe/InverseProblem
```

Then, to install ASSET run:

```julia
pkg> add https://github.com/SlitSpectroscopyBuddies/ASSET.git
```

## Usage

You can find in the `paper/` directory a submission to JOSS describing the package scientific background as well as results examples on JWST data. The notebook in `demo/` demonstrate how to use the package with a toy model.


## References

Additionally, the `ASSET` algorithm is described in the following papers:

*Characterization of stellar companions from high-contrast long-slit spectroscopy data -  The EXtraction Of SPEctrum of COmpanion (EXOSPECO) algorithm*, S.  Thé, É.  Thiébaut, L.  Denis, T.  Wanner, R. Thiébaut, M.  Langlois and F.  Soulez, A&A, 678 ,(2023) , A77 (**DOI:** https://doi.org/10.1051/0004-6361/202245565)

*An improved spectral extraction method for JWST/NIRSpec fixed slit observations*, L.  Denneulin, A.  Guilbert-Lepoutre, M.  Langlois, S.  Thé, E.  Thiébaut, B. J.  Holler and P.  Ferruit, A&A, 679  (2023) A63 (**DOI:** https://doi.org/10.1051/0004-6361/202346998)
