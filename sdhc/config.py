# -*- coding: utf-8 -*-
# Kimmo Sääskilahti, 2021
import logging


def setup_logger():
    logging.basicConfig(
        level=logging.INFO, format="%(name)s [%(levelname)s]: %(message)s"
    )


setup_logger()
logger = logging.getLogger("sdhc")
