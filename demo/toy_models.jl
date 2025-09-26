using DelimitedFiles
using InterpolationKernels
using LinearInterpolators

using ASSET

ker = CatmullRomSpline(Float64, Flat)


hg=ASSET.chromGaussianPSF(readdlm(target_LRS*"gaussian_parameters.txt")[1])
hm=ASSET.chromMoffatPSF(readdlm(target_LRS*"moffat_parameters.txt")[:])

h1=readdlm(target_LRS*"nonparametric_parameters_at_order_1.txt")
order1=size(h1)[2]
shift1=0.7745977082683768
psf_size=size(h1)[1]
hph=1e0
wh=cat(ones(psf_size,order1),zeros(psf_size,order1),dims=3)
Rh=WeightedTikhonov(hph, wh);
hnp1=ASSET.SeriesExpansionPSF(h1, shift1,ker ,Rh)

h2=readdlm(target_LRS*"nonparametric_parameters_at_order_2.txt")
order2=size(h2)[2]
shift2=0.5855828478647512
psf_size=size(h2)[1]
hph=1e0
wh=cat(ones(psf_size,order2),zeros(psf_size,order2),dims=3)
Rh=WeightedTikhonov(hph, wh);
hnp2=ASSET.SeriesExpansionPSF(h2, shift2,ker ,Rh)

(l,z) = readdlm("p330e_for_toymodel.txt")
