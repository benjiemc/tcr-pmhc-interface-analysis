import logging

LOG_MAP = {
    'error': logging.ERROR,
    'warning': logging.WARNING,
    'info': logging.INFO,
    'debug': logging.DEBUG,
}


def setup_logger(logger, level='warning'):
    '''Setup logger for command line applications.
    
    Args:
		logger: logger object from logging.getLogger(__name__) call
		level: one of 'error', 'warning', 'info', or 'debug'

    '''
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s: %(message)s'))
    logger.addHandler(handler)

    logger.setLevel(LOG_MAP[level])
