import os
import logging


def create_logger():
    """
        Allows for differnt log levels:
        logger.debug('This message should go to the log file')
        logger.info('So should this')
        logger.warning('And this, too')
        logger.error('And non-ASCII stuff, too, like Øresund and Malmö')

        log_level: logging.DEBUG or logging.INFO or logging.WARNING etc.
    """
    log_level = logging.INFO

    logger = logging.getLogger('Logger')
    logger.setLevel(log_level)

    logger.propagate = False

    ch = logging.StreamHandler()
    ch.setLevel(log_level)

    fh = logging.FileHandler('prototype.log', mode="w")
    fh.setLevel(log_level)

    formatter = logging.Formatter('%(levelname)s [%(asctime)s] [%(module)s:%(lineno)d]: %(message)s')

    ch.setFormatter(formatter)
    fh.setFormatter(formatter)

    logger.addHandler(ch)
    logger.addHandler(fh)

    return logger


logger = create_logger()
