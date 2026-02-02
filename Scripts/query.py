from GEO_fetch import GEOFetcher, GEOConfig

config = GEOConfig(
    query='"myopia"[All Fields] AND ("red light" OR photobiomodulation)',
    out_csv="geo_results.csv",
    download_dir="downloads",
    download_types=["matrix", "soft", "raw"],
    email="you@example.com",
    api_key="YOUR_KEY",
    retmax=200
)

fetcher = GEOFetcher(config)
fetcher.run()
