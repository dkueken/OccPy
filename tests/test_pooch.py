import pooch

# Download data from envidat and setup data folder
p = pooch.Pooch(base_url="doi:10.5281/zenodo.17750604/", path=r"../data/zenodo_test")
p.load_registry_from_doi()
p.fetch("MLS_demo.zip", processor=pooch.Unzip(members=["MLS_demo"]))
pc_dir = p.path / "MLS_demo.zip.unzip" / "LAZ"
print(pc_dir)
