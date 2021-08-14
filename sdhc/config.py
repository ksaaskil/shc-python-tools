import logging


def setup_logger():
    logging.basicConfig(
        level=logging.INFO, format="%(name)s [%(levelname)s]: %(message)s"
    )
    # logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


setup_logger()
logger = logging.getLogger("sdhc")
