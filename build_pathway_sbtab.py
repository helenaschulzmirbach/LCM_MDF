"""script to convert a simple network SBtab to a Pathway SBtab."""
# The MIT License (MIT)
#
# Copyright (c) 2013 Weizmann Institute of Science
# Copyright (c) 2018-2020 Institute for Molecular Systems Biology,
# ETH Zurich
# Copyright (c) 2018-2020 Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import argparse
from equilibrator_api import ComponentContribution, Q_
from equilibrator_pathway import Pathway


parser = argparse.ArgumentParser(
    description="Script for building a pathway configuration SBtab (for ECM or MDF) "
    "from a network SBtab."
)
parser.add_argument(
    "--ecm",
    action="store_true",
    help="make an ECM model (default: MDF)",
)
parser.add_argument("input_sbtab", type=str, help="Path to input network SBtab")
parser.add_argument("output_sbtab", type=str, help="Path to output pathway SBtab")
args = parser.parse_args()
comp_contrib = ComponentContribution()
pp = Pathway.from_network_sbtab(args.input_sbtab, comp_contrib=comp_contrib)
pp.to_sbtab().write(args.output_sbtab)
