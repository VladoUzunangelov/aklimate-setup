# FEB 2020	chrisw
# convert each row of tsv into stringified JSON lists

# imports
import logging
from optparse import OptionParser
import sys
import json
import os
import re

import datetime

# global vars
verbose = True

# methods and functions


def getOptions():
    "parse options"
    usage_text = []
    usage_text.append("%prog [options] < input")

    description_text = []
    description_text.append(
        "Output to stdout an AKLIMTE features file as stringified JSON object.")
    description_text.append(
        "Uniqueness of list members is not checked.")
    description_text.append(
        "spec: https://www.synapse.org/#!Synapse:syn8011998/wiki/600601")

    parser = OptionParser(usage="\n".join(usage_text),
                          description="\n".join(description_text))

    parser.add_option("-v", action="store_true", default=False,
                      dest="verbose", help="Switch for verbose mode.")

    (options, args) = parser.parse_args()

    return (options, args, parser)


def checkRequiredArguments(opts, parser):
    missing_options = []
    for option in parser.option_list:
        if re.match(r'^\[REQUIRED\]', option.help) and eval('opts.' + option.dest) is None:
            missing_options.extend(option._long_opts)
    if len(missing_options) > 0:
        parser.error('Missing REQUIRED parameters: ' + str(missing_options))


def getNow():
    """
    Get a datetime object for utc NOW.
    Convert to ISO 8601 format with datetime.isoformat()
    """
    now = datetime.datetime.utcnow()
    return now


def getTimeDelta(startDatetime):
    """
    get a timedelta object. Get seconds elapsed with timedelta.total_seconds().
    """
    endDatetime = datetime.datetime.utcnow()
    timedeltaObj = endDatetime - startDatetime
    return timedeltaObj


def compactJson(obj):
    s = json.dumps(obj, sort_keys=True, separators=(',', ':'))
    return s


def processLines(matrixFileLines):
    logging.debug("%d lines" % len(matrixFileLines))

    for line in matrixFileLines:
        fields = line.rstrip("\r\n").split("\t")

        s = compactJson(fields)
        sys.stdout.write("%s\n" % (s))
    return None


#:####################################


def main():
    startTime = getNow()
    (options, args, parser) = getOptions()

    checkRequiredArguments(options, parser)

    if options.verbose:
        logLevel = logging.DEBUG
    else:
        logLevel = logging.INFO
    # logfileName = os.path.basename(__file__).replace(".py", ".log")
    logFormat = "%(asctime)s %(levelname)s %(funcName)s:%(lineno)d %(message)s"
    logging.basicConfig(level=logLevel, format=logFormat)
    # setupLogging(logfileName, logFormat, logLevel)

    logging.debug('options:\t%s' % (str(options)))
    logging.debug('args:\t%s' % (str(args)))

    # DO SOMETHING HERE
    inputLines = sys.stdin.readlines()
    logging.debug("%d lines read from input" % (len(inputLines)))

    processLines(inputLines)

    runTime = getTimeDelta(startTime).total_seconds()
    logging.info("%s ran for %s s." %
                 (os.path.basename(__file__), str(runTime)))
    logging.shutdown()
    return None


# main program section
if __name__ == "__main__":
    main()
