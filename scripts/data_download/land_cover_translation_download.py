import requests
from conf import glcc_dir, url_glcc_dict

glcc_dir.mkdir(parents=True, exist_ok=True)
for region, url in url_glcc_dict.items():
    print(region)
    print(url)
    response_query = requests.get(url)
    with open(glcc_dir / f"land_cover_translation_glcc_{region}.txt", "wb") as file:
        file.write(response_query.content)
