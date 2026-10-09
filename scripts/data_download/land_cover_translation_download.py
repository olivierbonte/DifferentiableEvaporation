import requests
from conf import glcc_dir, logger, url_glcc_dict


def main():
    glcc_dir.mkdir(parents=True, exist_ok=True)
    for region, url in url_glcc_dict.items():
        logger.info(region)
        logger.info(url)
        response_query = requests.get(url)
        with open(glcc_dir / f"land_cover_translation_glcc_{region}.txt", "wb") as file:
            file.write(response_query.content)


if __name__ == "__main__":
    main()
