using Documenter, ASSET, DocumenterInterLinks

makedocs(sitename="ASSET documentation", 
         repo=Remotes.GitHub("SlitSpectroscopyBuddies", "ASSET"))
         
deploydocs(
    repo = "github.com/SlitSpectroscopyBuddies/ASSET.git",
)


links = InterLinks(
    "PointSpreadFunctions" => "https://emmt.github.io/PointSpreadFunctions.jl/"
)
