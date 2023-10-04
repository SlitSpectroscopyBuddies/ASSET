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
# FIXME : choisir une licence !!! Copyright 2023 Samuel Thé & Laurence Denneulin 
#
module ASSET

export 
    AbstractBkg,
    BkgMdl,
    bkg_step!,
    CalibratedData,
    extract_spectrum!,
    fit_psf_center!,
    fit_psf_params!,
    fit_spectrum!,
    loss,
    mask_object,
    object_step!,
    ParametrizedBkgMdl,
    psf_map!,
    solve_analytic!,
    solve_vmlmb!

import Base: +, -, axes, size, eltype, show, copy
using InterpolationKernels
using InverseProblem
using LazyAlgebra
using LinearInterpolators
using OptimPackNextGen
import OptimPackNextGen.Powell.Bobyqa
using PointSpreadFunctions
import PointSpreadFunctions:
    parameters, getfwhm
using SparseArrays



include("types.jl")
include("algorithms.jl")
include("bkg_models.jl")
include("psf_models.jl")


end # module ASSET
