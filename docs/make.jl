using Documenter, ASSET

makedocs(sitename="ASSET documentation", 
         repo=Remotes.GitHub("SlitSpectroscopyBuddies", "ASSET"))
         
deploydocs(
    repo = "github.com/SlitSpectroscopyBuddies/ASSET.git",
)
