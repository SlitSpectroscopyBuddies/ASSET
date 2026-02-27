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
    CalibratedData,
    chromGaussianPSF,
    chromwmwGaussianPSF,
    chromMoffatPSF,
    chromwmwMoffatPSF,
    ChromaticSeriesExpansionsInterpolator,
    extract_spectrum!,
    fit_bkg!,
    fit_psf_center!,
    fit_psf_params,
    fit_psf_shift,
    fit_spectrum!,
    fit_spectrum_and_psf!,
    loss,
    mask_object,
    NonParametricPSF,
    psf_map!,
    psf_map,
    ParametricPSF,
    SeriesExpansionPSF,
    solve_analytic,
    solve_vmlmb

import Base: +, -, axes, size, eltype, show, copy
using AMORS
using InterpolationKernels
using InverseProblem
import InverseProblem: test_tol
using LazyAlgebra
import LazyAlgebra: Mapping, vcreate, vcopy, apply!
using LinearInterpolators
using OptimPackNextGen
import OptimPackNextGen.Brent
using PointSpreadFunctions
import PointSpreadFunctions:
    parameters, getfwhm
#using PowellMethods FIXME: Use PowellMethods instead of PRIMA, invest bound issues
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
