import datetime


END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
TEAL = '\033[36m'

# ANSI color codes
#RED = '\033[31m'
#GREEN = '\033[32m'
#ORANGE = '\033[33m'
#BLUE = '\033[34m'
#MAGENTA = '\033[35m'


def bold(text):
    return BOLD + text + END_FORMATTING


def bold_teal(text):
    return TEAL + BOLD + text + END_FORMATTING


def bold_teal_underline(text):
    return TEAL + BOLD + UNDERLINE + text + END_FORMATTING


def get_timestamp():
    return '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())


def get_ascii_art():
    (
    bold_teal(r" _______     _       _______                  _                ") + '\n' +
    bold_teal(r"|__   __|   | |     |__   __|                | |               ") + '\n' +
    bold_teal(r"   | |  ___ | |  ___   | | _ __   __ _   ___ | | __  ___  _ __ ") + '\n' +
    bold_teal(r"   | | / _ \| | / _ \  | || '__| / _` | / __|| |/ / / _ \| '__|") + '\n' +
    bold_teal(r"   | ||  __/| || (_) | | || |   | (_| || (__ |   < |  __/| |   ") + '\n' +
    bold_teal(r"   |_| \___||_| \___/  |_||_|    \__,_| \___||_|\_\ \___||_|   ") + '\n'
    )
