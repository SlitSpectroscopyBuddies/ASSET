using Documenter, ASSET, DocumenterInterLinks

makedocs(sitename="ASSET documentation", 
         repo=Remotes.GitHub("SlitSpectroscopyBuddies", "ASSET"))
         
deploydocs(
    repo = "github.com/SlitSpectroscopyBuddies/ASSET.git",
)


#links = InterLinks(
#    "PointSpreadFunctions" => "https://documenter.juliadocs.org/stable/objects.inv"
#)
