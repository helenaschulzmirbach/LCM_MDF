"""Command-line script for pathway analysis."""
# The MIT License (MIT)
#
# Copyright (c) 2013 Weizmann Institute of Science
# Copyright (c) 2018 Institute for Molecular Systems Biology,
# ETH Zurich
# Copyright (c) 2018 Novo Nordisk Foundation Center for Biosustainability,
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
import logging

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sbtab.SBtab import SBtabDocument

from equilibrator_pathway import EnzymeCostModel, Pathway, ThermodynamicModel

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Calculate the Max-min Driving Force (MDF) of a pathway."
    )
    parser.add_argument(
        "infile",
        type=argparse.FileType("r"),
        help="path to input file containing reactions",
    )
    parser.add_argument(
        "outfile", type=str, help="prefix for the output TSV and PDF files"
    )
    logging.getLogger().setLevel(logging.WARNING)

    args = parser.parse_args()

    sbtabdoc = SBtabDocument(
        "pathway", args.infile.read(), "pathway.tsv", definitions_file=None
    )

    config_sbtab = sbtabdoc.get_sbtab_by_id("Configuration")

    pp = Pathway.from_sbtab(sbtabdoc)
    algorithm = pp.config_dict.get("algorithm", "MDF")
    if algorithm == "MDF":
        mdf_model = ThermodynamicModel.from_sbtab(sbtabdoc)
        mdf_sol = mdf_model.mdf_analysis()
        mdf_sol.to_sbtab().write(args.outfile + ".tsv")

        with PdfPages(open(args.outfile + ".pdf", "wb")) as pdf:
            fig1, ax = plt.subplots(1, 1, figsize=(10, 7))
            mdf_sol.plot_concentrations(ax=ax)
            ax.axes.yaxis.grid(True, which="major")
            fig1.tight_layout()
            pdf.savefig(fig1)

            fig2, ax = plt.subplots(1, 1, figsize=(10, 7))
            mdf_sol.plot_driving_forces(ax=ax)
            ax.axes.xaxis.grid(True, which="major")
            fig2.tight_layout()
            pdf.savefig(fig2)

            # find a way to "print" these tables to the PDF file
            net_rxn = pp.net_reaction
            rxn_df = mdf_sol.reaction_df
            cpd_df = mdf_sol.compound_df
            #added by Helena: add net reaction to PDF
            formula = pp.net_reaction_formula
            #print(formula)
            liste = list(formula)
            liste.insert(100,'\n')

            liste.insert(200,'\n')
            formula2 = ''.join(liste)

            fig3, ax = plt.subplots(1, 1, figsize=(10, 7))
            ax.text(0,1,formula2)
            ax.set_xlim([0,5])

            ax.set_ylim([1,1.2])
            ax.set_xticks([])
            ax.set_yticks([])
            ax.spines[['left','right','top','bottom']].set_visible(False)

            pdf.savefig(fig3)

    elif algorithm == "ECM":
        ecm_model = EnzymeCostModel.from_sbtab(sbtabdoc)
        ecm_sol = ecm_model.optimize_ecm()
        ecm_sol.to_sbtab().write(args.outfile + ".tsv")

        with PdfPages(open(args.outfile + ".pdf", "wb")) as pdf:
            fig1, ax = plt.subplots(1, 1, figsize=(7, 4))
            ax.set_title("ECM solution")
            ecm_sol.plot_enzyme_demand_breakdown(ax, plot_measured=True)
            ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
            ax.axes.yaxis.grid(True, which="major")
            fig1.tight_layout()
            pdf.savefig(fig1)

            fig2, ax = plt.subplots(1, 1, figsize=(5, 5))
            ecm_sol.plot_volumes_pie(ax=ax)
            fig2.tight_layout()
            pdf.savefig(fig2)
    else:
        raise ValueError(
            f"Unknown algorithm: {algorithm}. Please use 'MDF' " f"or 'ECM'"
        )
