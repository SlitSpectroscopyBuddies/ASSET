#
# ASSET.jl --
#
# ASSET is an Adaptive Slit Spectal Extraction Tool. 
# This package aims at giving low-level general tools 
# to the extraction of spectrum of an unresolved object 
# with slit spectroscopy. 
#
#-------------------------------------------------------------------------------
# 
# Copyright 2023 SlitSpectroscopyBuddies organization (Samuel Thé & Laurence Denneulin) 
#
module ASSET

export 
    AbstractBkg,
    BkgMdl,
    bkg_step!,
    CalibratedData,
    chromGaussianPSF,
    chromwmwGaussianPSF,
    ChromaticSeriesExpansionsInterpolator,
    extract_spectrum!,
    fit_psf_center!,
    fit_psf_params!,
    fit_spectrum!,
    fit_psf_shift,
    loss,
    mask_object,
    fit_spectrum_and_psf!,
    oneDimensionalPSF,
    ParametrizedBkgMdl,
    psf_map!,
    psf_map,
    SeriesExpansionPSF,
    solve_analytic!,
    solve_vmlmb!

import Base: +, -, axes, size, eltype, show, copy
using AMORS
using InterpolationKernels
using InverseProblem
using LazyAlgebra
import LazyAlgebra: Mapping, vcreate, vcopy, apply!
using LinearInterpolators
using PyPlot
using OptimPackNextGen
import OptimPackNextGen.Brent
using PointSpreadFunctions
import PointSpreadFunctions:
    parameters, getfwhm
#using PowellMethods
using PRIMA
using SparseArrays



include("types.jl")
include("non_parametric_fitting.jl")
include("parametric_fitting.jl")
include("bkg_models.jl")
include("psf_models.jl")
include("utils.jl")
include("algorithms.jl")

end # module ASSET
