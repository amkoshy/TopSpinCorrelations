#!/usr/bin/env python
import ROOT
import os
import math
import numpy as np
import argparse

def roundup(x, n=1):
    res = math.ceil(x/n)*n
    return res

#parser = argparse.ArgumentParser(description=__doc__)
#parser.add_argument('-y', '--era', dest='era',
#                        action='store', default='',
#                        help='Era (example: 2018, 2017, 2016')
#parser.add_argument('-c', '--channel', dest='channel',
#                        action='store', default='',
#                        help='Channel (example: ee, emu, mumu, combined')
#opts, opts_unknown = parser.parse_known_args()

eras = [
"2016preVFP",
"2016postVFP",
"2016",
"2017",
"2018",
"fullRun2UL",
]

channels = [
"ee",
"emu",
"mumu",
"combined",
]

vars2Dlist = ["TTBarMass", "ScatteringAngle_TTBarFrame", "ToppT", "ExtraJets"]



varsdict = {"ttbarmass":["0","450","600","800","inf"],"topscatteringangle":["m1","mhalf","0","phalf","p1"]}

controlplotlist = [
    
    # UnSymmetrized SpinCorr

    ['HypAntiLeptonBk', 'Standard', 'Events', '"cos#theta_{1}^{k}"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypLeptonBk', 'Standard', 'Events', '"cos#theta_{2}^{k}"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypAntiLeptonBj', 'Standard', 'Events', '"cos#theta_{1}^{k*}"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypLeptonBj', 'Standard', 'Events', '"cos#theta_{2}^{k*}"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypAntiLeptonBr', 'Standard', 'Events', '"cos#theta_{1}^{r}"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypLeptonBr', 'Standard', 'Events', '"cos#theta_{2}^{r}"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypAntiLeptonBq', 'Standard', 'Events', '"cos#theta_{1}^{r*}"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypLeptonBq', 'Standard', 'Events', '"cos#theta_{2}^{r*}"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypAntiLeptonBn', 'Standard', 'Events', '"cos#theta_{1}^{n}"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypLeptonBn', 'Standard', 'Events', '"cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],

    #['HypLLBarBPnn', 'Standard', 'Events', '"cos#theta_{1}^{n}+cos#theta_{2}^{n}"', '20', '1', '0', '0', '0', '40e3', '-2.0', '1.99999', '0', '0'],
    #['HypLLBarBMnn', 'Standard', 'Events', '"cos#theta_{1}^{n}-cos#theta_{2}^{n}"', '20', '1', '0', '0', '0', '40e3', '-2.0', '1.99999', '0', '0'],
    #['HypLLBarBPrr', 'Standard', 'Events', '"cos#theta_{1}^{r}+cos#theta_{2}^{r}"', '20', '1', '0', '0', '0', '40e3', '-2.0', '1.99999', '0', '0'],
    #['HypLLBarBMrr', 'Standard', 'Events', '"cos#theta_{1}^{r}-cos#theta_{2}^{r}"', '20', '1', '0', '0', '0', '40e3', '-2.0', '1.99999', '0', '0'],
    #['HypLLBarBPkk', 'Standard', 'Events', '"cos#theta_{1}^{k}+cos#theta_{2}^{k}"', '20', '1', '0', '0', '0', '40e3', '-2.0', '1.99999', '0', '0'],
    #['HypLLBarBMkk', 'Standard', 'Events', '"cos#theta_{1}^{k}-cos#theta_{2}^{k}"', '20', '1', '0', '0', '0', '40e3', '-2.0', '1.99999', '0', '0'],
    #['HypLLBarBPjj', 'Standard', 'Events', '"cos#theta_{1}^{k*}+cos#theta_{2}^{k*}"', '20', '1', '0', '0', '0', '40e3', '-2.0', '1.99999', '0', '0'],
    #['HypLLBarBMjj', 'Standard', 'Events', '"cos#theta_{1}^{k*}-cos#theta_{2}^{k*}"', '20', '1', '0', '0', '0', '40e3', '-2.0', '1.99999', '0', '0'],
    #['HypLLBarBPqq', 'Standard', 'Events', '"cos#theta_{1}^{r*}+cos#theta_{2}^{r*}"', '20', '1', '0', '0', '0', '40e3', '-2.0', '1.99999', '0', '0'],
    #['HypLLBarBMqq', 'Standard', 'Events', '"cos#theta_{1}^{r*}-cos#theta_{2}^{r*}"', '20', '1', '0', '0', '0', '40e3', '-2.0', '1.99999', '0', '0'],

    ['HypLLBarCkk', 'Standard', 'Events', '"cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypLLBarCrr', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypLLBarCnn', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypLLBarCkj', 'Standard', 'Events', '"cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k*}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypLLBarCrq', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r*}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],

    ['HypLLBarChan', 'Standard', 'Events', '"+cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypLLBarCsca', 'Standard', 'Events', '"-cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypLLBarCtra', 'Standard', 'Events', '"-cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypLLBarCkjL', 'Standard', 'Events', '"-cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypLLBarCrqL', 'Standard', 'Events', '"-cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r*}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],

    #['HypLLBarCrk', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    #['HypLLBarCkr', 'Standard', 'Events', '"cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{r}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    #['HypLLBarCnr', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    #['HypLLBarCrn', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    #['HypLLBarCnk', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    #['HypLLBarCkn', 'Standard', 'Events', '"cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    #['HypLLBarCrj', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k*}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    #['HypLLBarCjr', 'Standard', 'Events', '"cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{r}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],

    #['HypLLBarCqj', 'Standard', 'Events', '"cos#theta_{1}^{r*}#kern[-0.8]{ }cos#theta_{2}^{k*}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    #['HypLLBarCjq', 'Standard', 'Events', '"cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{r*}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    #['HypLLBarCnq', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r*}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    #['HypLLBarCqn', 'Standard', 'Events', '"cos#theta_{1}^{r*}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    #['HypLLBarCnj', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k*}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],

    ['HypLLBarCPrk', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{r}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    ['HypLLBarCMrk', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{r}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    ['HypLLBarCPnr', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    ['HypLLBarCMnr', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    ['HypLLBarCPnk', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    ['HypLLBarCMnk', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],

    ['HypLLBarCPrj', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{r}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    ['HypLLBarCMrj', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{r}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    #['HypLLBarCPqj', 'Standard', 'Events', '"cos#theta_{1}^{r*}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{r*}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    #['HypLLBarCMqj', 'Standard', 'Events', '"cos#theta_{1}^{r*}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{r*}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    #['HypLLBarCPnq', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r*}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{r*}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    #['HypLLBarCMnq', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r*}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{r*}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    #['HypLLBarCPnj', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    #['HypLLBarCMnj', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],

    ['HypLLBarCrkP', 'Standard', 'Events', '"-cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypLLBarCrkM', 'Standard', 'Events', '"-cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypLLBarCnrP', 'Standard', 'Events', '"-cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{n}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypLLBarCnrM', 'Standard', 'Events', '"-cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{n}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypLLBarCnkP', 'Standard', 'Events', '"-cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{n}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypLLBarCnkM', 'Standard', 'Events', '"-cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{n}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],

    #['HypLLBarCqjP', 'Standard', 'Events', '"-cos#theta_{1}^{r*}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{r*}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    #['HypLLBarCqjM', 'Standard', 'Events', '"-cos#theta_{1}^{r*}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{r*}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    #['HypLLBarCnqP', 'Standard', 'Events', '"-cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r*}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r*}#kern[-0.8]{ }cos#theta_{2}^{n}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{k*}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    #['HypLLBarCnqM', 'Standard', 'Events', '"-cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r*}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{r*}#kern[-0.8]{ }cos#theta_{2}^{n}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{k*}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    #['HypLLBarCnjP', 'Standard', 'Events', '"-cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{n}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r*}#kern[-0.8]{ }cos#theta_{2}^{r*}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    #['HypLLBarCnjP', 'Standard', 'Events', '"-cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{n}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r*}#kern[-0.8]{ }cos#theta_{2}^{r*}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],

    ['HypLLBarcHel', 'Standard', 'Events', '"cos#varphi"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99999', '0', '0'],
    ['HypLLBarcLab', 'Standard', 'Events', '"cos#varphi_{lab}"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99999', '0', '0'],

    ['HypLLBarkNorm', 'Standard', 'Events', '"sin#theta (k)"', '10', '1', '0', '0', '0', '0', '-2.0', '1.99999', '0', '0'],
    ['HypLLBarrNorm', 'Standard', 'Events', '"sin#theta (r)"', '10', '1', '0', '0', '0', '0', '-2.0', '1.99999', '0', '0'],
    ['HypLLBarnNorm', 'Standard', 'Events', '"sin#theta (n)"', '10', '1', '0', '0', '0', '0', '-2.0', '1.99999', '0', '0'],
    ['HypLLBarjNorm', 'Standard', 'Events', '"sin#theta (k*)"', '10', '1', '0', '0', '0', '0', '-2.0', '1.99999', '0', '0'],
    ['HypLLBarqNorm', 'Standard', 'Events', '"sin#theta (r*)"', '10', '1', '0', '0', '0', '0', '-2.0', '1.99999', '0', '0'],


    # Symmetrized SpinCorr

    ['HypSymAntiLeptonBk', 'Standard', 'Events', '"cos#theta_{1}^{k}"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypSymLeptonBk', 'Standard', 'Events', '"cos#theta_{2}^{k}"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypSymAntiLeptonBj', 'Standard', 'Events', '"cos#theta_{1}^{k*}"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypSymLeptonBj', 'Standard', 'Events', '"cos#theta_{2}^{k*}"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypSymAntiLeptonBr', 'Standard', 'Events', '"cos#theta_{1}^{r}"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypSymLeptonBr', 'Standard', 'Events', '"cos#theta_{2}^{r}"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypSymAntiLeptonBq', 'Standard', 'Events', '"cos#theta_{1}^{r*}"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypSymLeptonBq', 'Standard', 'Events', '"cos#theta_{2}^{r*}"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypSymAntiLeptonBn', 'Standard', 'Events', '"cos#theta_{1}^{n}"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypSymLeptonBn', 'Standard', 'Events', '"cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],

    #['HypSymLLBarBPnn', 'Standard', 'Events', '"cos#theta_{1}^{n}+cos#theta_{2}^{n}"', '20', '1', '0', '0', '0', '40e3', '-2.0', '1.99999', '0', '0'],
    #['HypSymLLBarBMnn', 'Standard', 'Events', '"cos#theta_{1}^{n}-cos#theta_{2}^{n}"', '20', '1', '0', '0', '0', '40e3', '-2.0', '1.99999', '0', '0'],
    #['HypSymLLBarBPrr', 'Standard', 'Events', '"cos#theta_{1}^{r}+cos#theta_{2}^{r}"', '20', '1', '0', '0', '0', '40e3', '-2.0', '1.99999', '0', '0'],
    #['HypSymLLBarBMrr', 'Standard', 'Events', '"cos#theta_{1}^{r}-cos#theta_{2}^{r}"', '20', '1', '0', '0', '0', '40e3', '-2.0', '1.99999', '0', '0'],
    #['HypSymLLBarBPkk', 'Standard', 'Events', '"cos#theta_{1}^{k}+cos#theta_{2}^{k}"', '20', '1', '0', '0', '0', '40e3', '-2.0', '1.99999', '0', '0'],
    #['HypSymLLBarBMkk', 'Standard', 'Events', '"cos#theta_{1}^{k}-cos#theta_{2}^{k}"', '20', '1', '0', '0', '0', '40e3', '-2.0', '1.99999', '0', '0'],
    #['HypSymLLBarBPjj', 'Standard', 'Events', '"cos#theta_{1}^{k*}+cos#theta_{2}^{k*}"', '20', '1', '0', '0', '0', '40e3', '-2.0', '1.99999', '0', '0'],
    #['HypSymLLBarBMjj', 'Standard', 'Events', '"cos#theta_{1}^{k*}-cos#theta_{2}^{k*}"', '20', '1', '0', '0', '0', '40e3', '-2.0', '1.99999', '0', '0'],
    #['HypSymLLBarBPqq', 'Standard', 'Events', '"cos#theta_{1}^{r*}+cos#theta_{2}^{r*}"', '20', '1', '0', '0', '0', '40e3', '-2.0', '1.99999', '0', '0'],
    #['HypSymLLBarBMqq', 'Standard', 'Events', '"cos#theta_{1}^{r*}-cos#theta_{2}^{r*}"', '20', '1', '0', '0', '0', '40e3', '-2.0', '1.99999', '0', '0'],

    ['HypSymLLBarCkk', 'Standard', 'Events', '"cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCrr', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCnn', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCkj', 'Standard', 'Events', '"cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k*}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCrq', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r*}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],

    ['HypSymLLBarChan', 'Standard', 'Events', '"+cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCsca', 'Standard', 'Events', '"-cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCtra', 'Standard', 'Events', '"-cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCkjL', 'Standard', 'Events', '"-cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCrqL', 'Standard', 'Events', '"-cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r*}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],

    ['HypSymLLBarCrk', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCkr', 'Standard', 'Events', '"cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{r}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCnr', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCrn', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCnk', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCkn', 'Standard', 'Events', '"cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCrj', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k*}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCjr', 'Standard', 'Events', '"cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{r}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],

    #['HypSymLLBarCqj', 'Standard', 'Events', '"cos#theta_{1}^{r*}#kern[-0.8]{ }cos#theta_{2}^{k*}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    #['HypSymLLBarCjq', 'Standard', 'Events', '"cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{r*}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    #['HypSymLLBarCnq', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r*}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    #['HypSymLLBarCqn', 'Standard', 'Events', '"cos#theta_{1}^{r*}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    #['HypSymLLBarCnj', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k*}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],

    ['HypSymLLBarCPrk', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{r}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCMrk', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{r}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCPnr', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCMnr', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCPnk', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCMnk', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],

    ['HypSymLLBarCPrj', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{r}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCMrj', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{r}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    #['HypSymLLBarCPqj', 'Standard', 'Events', '"cos#theta_{1}^{r*}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{r*}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    #['HypSymLLBarCMqj', 'Standard', 'Events', '"cos#theta_{1}^{r*}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{r*}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    #['HypSymLLBarCPnq', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r*}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{r*}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    #['HypSymLLBarCMnq', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r*}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{r*}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    #['HypSymLLBarCPnj', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    #['HypSymLLBarCMnj', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],

    ['HypSymLLBarCrkP', 'Standard', 'Events', '"-cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCrkM', 'Standard', 'Events', '"-cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCnrP', 'Standard', 'Events', '"-cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{n}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCnrM', 'Standard', 'Events', '"-cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{n}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCnkP', 'Standard', 'Events', '"-cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{n}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarCnkM', 'Standard', 'Events', '"-cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{n}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],

    #['HypSymLLBarCqjP', 'Standard', 'Events', '"-cos#theta_{1}^{r*}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{r*}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    #['HypSymLLBarCqjM', 'Standard', 'Events', '"-cos#theta_{1}^{r*}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{r*}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    #['HypSymLLBarCnqP', 'Standard', 'Events', '"-cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r*}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r*}#kern[-0.8]{ }cos#theta_{2}^{n}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{k*}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    #['HypSymLLBarCnqM', 'Standard', 'Events', '"-cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r*}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{r*}#kern[-0.8]{ }cos#theta_{2}^{n}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{k*}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    #['HypSymLLBarCnjP', 'Standard', 'Events', '"-cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{n}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r*}#kern[-0.8]{ }cos#theta_{2}^{r*}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    #['HypSymLLBarCnjP', 'Standard', 'Events', '"-cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{n}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r*}#kern[-0.8]{ }cos#theta_{2}^{r*}"', '10', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],

    ['HypSymLLBarcHel', 'Standard', 'Events', '"cos#varphi"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99999', '0', '0'],
    ['HypSymLLBarcLab', 'Standard', 'Events', '"cos#varphi_{lab}"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99999', '0', '0'],

    ['HypSymLLBarkNorm', 'Standard', 'Events', '"sin#theta (k)"', '10', '1', '0', '0', '0', '0', '-2.0', '1.99999', '0', '0'],
    ['HypSymLLBarrNorm', 'Standard', 'Events', '"sin#theta (r)"', '10', '1', '0', '0', '0', '0', '-2.0', '1.99999', '0', '0'],
    ['HypSymLLBarnNorm', 'Standard', 'Events', '"sin#theta (n)"', '10', '1', '0', '0', '0', '0', '-2.0', '1.99999', '0', '0'],
    ['HypSymLLBarjNorm', 'Standard', 'Events', '"sin#theta (k*)"', '10', '1', '0', '0', '0', '0', '-2.0', '1.99999', '0', '0'],
    ['HypSymLLBarqNorm', 'Standard', 'Events', '"sin#theta (r*)"', '10', '1', '0', '0', '0', '0', '-2.0', '1.99999', '0', '0'],

    #['#', 'Name,', 'Extra,', 'axis', 'labels', '(y,x),', 'rebin,', 'do_dyscale,', 'logx,', 'logy,', 'ymin,', 'ymax,', 'xmin,', 'xmax,', 'nbins,', 'xbins,', 'bcs'],

    # Number of vertices
    #['HypvertMulti', 'Standard', 'Events', 'N_{vtx}', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0'],

    # Other

    ['HypLLBarDPhi', 'Standard', 'Events', '"|#Delta#phi_{l#bar{l}}|"', '157', '1', '0', '0', '0', '30e3', '0', '3.13999', '0', '0'],
    ['HypLLBarDEta', 'Standard', 'Events', '"|#Delta#eta_{l#bar{l}}|"', '12', '1', '0', '1', '0', '70000.0', '0', '4.99999', '0', '0'],

    ['HypTTBarMass', 'Standard', 'Toppairs', 'm_{t#bar{t}}', '60', '1', '0', '0', '0', '150e3', '300', '1999', '0', '0'],

    ['HypScatteringAngle_TTBarFrame', 'Standard', 'Events', '"cos#Theta"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99999', '0', '0'],
    ['HypScatteringAngle_LabFrame', 'Standard', 'Events', '"cos#Theta_{lab}"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99999', '0', '0'],

]



controlplotlist_fromtrees = [

    # UnSymmetrized SpinCorr

    ['Hyp_AntiLeptonBk', 'Standard', 'Events', '"cos#theta_{1}^{k}"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['Hyp_LeptonBk', 'Standard', 'Events', '"cos#theta_{2}^{k}"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['Hyp_AntiLeptonBj', 'Standard', 'Events', '"cos#theta_{1}^{k*}"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['Hyp_LeptonBj', 'Standard', 'Events', '"cos#theta_{2}^{k*}"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['Hyp_AntiLeptonBr', 'Standard', 'Events', '"cos#theta_{1}^{r}"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['Hyp_LeptonBr', 'Standard', 'Events', '"cos#theta_{2}^{r}"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['Hyp_AntiLeptonBq', 'Standard', 'Events', '"cos#theta_{1}^{r*}"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['Hyp_LeptonBq', 'Standard', 'Events', '"cos#theta_{2}^{r*}"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['Hyp_AntiLeptonBn', 'Standard', 'Events', '"cos#theta_{1}^{n}"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['Hyp_LeptonBn', 'Standard', 'Events', '"cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],

    ['Hyp_LLBarCkk', 'Standard', 'Events', '"cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['Hyp_LLBarCrr', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['Hyp_LLBarCnn', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['Hyp_LLBarCkj', 'Standard', 'Events', '"cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k*}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['Hyp_LLBarCrq', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r*}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],

    ['Hyp_LLBarChan', 'Standard', 'Events', '"+cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['Hyp_LLBarCsca', 'Standard', 'Events', '"-cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['Hyp_LLBarCtra', 'Standard', 'Events', '"-cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['Hyp_LLBarCkjL', 'Standard', 'Events', '"-cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['Hyp_LLBarCrqL', 'Standard', 'Events', '"-cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r*}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],

    ['Hyp_LLBarCPrk', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{r}"', '1', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    ['Hyp_LLBarCMrk', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{r}"', '1', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    ['Hyp_LLBarCPnr', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    ['Hyp_LLBarCMnr', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    ['Hyp_LLBarCPnk', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    ['Hyp_LLBarCMnk', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],

    ['Hyp_LLBarCPrj', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{r}"', '1', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    ['Hyp_LLBarCMrj', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{r}"', '1', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],

    ['Hyp_LLBarCrkP', 'Standard', 'Events', '"-cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['Hyp_LLBarCrkM', 'Standard', 'Events', '"-cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['Hyp_LLBarCnrP', 'Standard', 'Events', '"-cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{n}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['Hyp_LLBarCnrM', 'Standard', 'Events', '"-cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{n}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['Hyp_LLBarCnkP', 'Standard', 'Events', '"-cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{n}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['Hyp_LLBarCnkM', 'Standard', 'Events', '"-cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{n}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],

    ['Hyp_LLBarcHel', 'Standard', 'Events', '"cos#varphi"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99999', '0', '0'],
    ['Hyp_LLBarcLab', 'Standard', 'Events', '"cos#varphi_{lab}"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99999', '0', '0'],

    ['Hyp_LLBarkNorm', 'Standard', 'Events', '"sin#theta (k)"', '1', '1', '0', '0', '0', '0', '-2.0', '1.99999', '0', '0'],
    ['Hyp_LLBarrNorm', 'Standard', 'Events', '"sin#theta (r)"', '1', '1', '0', '0', '0', '0', '-2.0', '1.99999', '0', '0'],
    ['Hyp_LLBarnNorm', 'Standard', 'Events', '"sin#theta (n)"', '1', '1', '0', '0', '0', '0', '-2.0', '1.99999', '0', '0'],
    ['Hyp_LLBarjNorm', 'Standard', 'Events', '"sin#theta (k*)"', '1', '1', '0', '0', '0', '0', '-2.0', '1.99999', '0', '0'],
    ['Hyp_LLBarqNorm', 'Standard', 'Events', '"sin#theta (r*)"', '1', '1', '0', '0', '0', '0', '-2.0', '1.99999', '0', '0'],


    # Symmetrized SpinCorr

    ['HypSym_AntiLeptonBk', 'Standard', 'Events', '"cos#theta_{1}^{k}"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypSym_LeptonBk', 'Standard', 'Events', '"cos#theta_{2}^{k}"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypSym_AntiLeptonBj', 'Standard', 'Events', '"cos#theta_{1}^{k*}"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypSym_LeptonBj', 'Standard', 'Events', '"cos#theta_{2}^{k*}"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypSym_AntiLeptonBr', 'Standard', 'Events', '"cos#theta_{1}^{r}"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypSym_LeptonBr', 'Standard', 'Events', '"cos#theta_{2}^{r}"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypSym_AntiLeptonBq', 'Standard', 'Events', '"cos#theta_{1}^{r*}"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypSym_LeptonBq', 'Standard', 'Events', '"cos#theta_{2}^{r*}"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypSym_AntiLeptonBn', 'Standard', 'Events', '"cos#theta_{1}^{n}"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],
    ['HypSym_LeptonBn', 'Standard', 'Events', '"cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99', '0', '0'],

    ['HypSym_LLBarCkk', 'Standard', 'Events', '"cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCrr', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCnn', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCkj', 'Standard', 'Events', '"cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k*}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCrq', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r*}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],

    ['HypSym_LLBarChan', 'Standard', 'Events', '"+cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCsca', 'Standard', 'Events', '"-cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCtra', 'Standard', 'Events', '"-cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCkjL', 'Standard', 'Events', '"-cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCrqL', 'Standard', 'Events', '"-cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r*}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],

    ['HypSym_LLBarCrk', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCkr', 'Standard', 'Events', '"cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{r}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCnr', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCrn', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCnk', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCkn', 'Standard', 'Events', '"cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCrj', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k*}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCjr', 'Standard', 'Events', '"cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{r}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],

    ['HypSym_LLBarCPrk', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{r}"', '1', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCMrk', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{r}"', '1', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCPnr', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCMnr', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCPnk', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCMnk', 'Standard', 'Events', '"cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],

    ['HypSym_LLBarCPrj', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{r}"', '1', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCMrj', 'Standard', 'Events', '"cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k*}#kern[-0.6]{ }-#kern[-0.4]{ }cos#theta_{1}^{k*}#kern[-0.8]{ }cos#theta_{2}^{r}"', '1', '1', '0', '0', '0', '40e3', '-1.0', '.99999', '0', '0'],

    ['HypSym_LLBarCrkP', 'Standard', 'Events', '"-cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCrkM', 'Standard', 'Events', '"-cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{n}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCnrP', 'Standard', 'Events', '"-cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{n}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCnrM', 'Standard', 'Events', '"-cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{r}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{n}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{k}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCnkP', 'Standard', 'Events', '"-cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{n}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarCnkM', 'Standard', 'Events', '"-cos#theta_{1}^{n}#kern[-0.8]{ }cos#theta_{2}^{k}#kern[-0.8]{ }+#kern[-0.6]{ }cos#theta_{1}^{k}#kern[-0.8]{ }cos#theta_{2}^{n}#kern[-0.8]{ }-#kern[-0.6]{ }cos#theta_{1}^{r}#kern[-0.8]{ }cos#theta_{2}^{r}"', '1', '1', '0', '0', '0', '80e3', '-1.0', '.99999', '0', '0'],

    ['HypSym_LLBarcHel', 'Standard', 'Events', '"cos#varphi"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99999', '0', '0'],
    ['HypSym_LLBarcLab', 'Standard', 'Events', '"cos#varphi_{lab}"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99999', '0', '0'],

    #['#', 'Name,', 'Extra,', 'axis', 'labels', '(y,x),', 'rebin,', 'do_dyscale,', 'logx,', 'logy,', 'ymin,', 'ymax,', 'xmin,', 'xmax,', 'nbins,', 'xbins,', 'bcs'],

    # Other

    ['Hyp_LLBarDPhi', 'Standard', 'Events', '"|#Delta#phi_{l#bar{l}}|"', '1', '1', '0', '0', '0', '30e3', '0', '3.13999', '0', '0'],
    ['Hyp_LLBarDEta', 'Standard', 'Events', '"|#Delta#eta_{l#bar{l}}|"', '1', '1', '0', '0', '0', '45000.0', '0', '4.99999', '0', '0'],

    ['Hyp_TTBarMass', 'Standard', 'Toppairs', 'm_{t#bar{t}}', '1', '1', '0', '0', '0', '150e3', '300', '1999', '0', '0'],
    ['Hyp_ScatteringAngle_TTBarFrame', 'Standard', 'Events', '"cos#Theta"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99999', '0', '0'],
    ['Hyp_ScatteringAngle_LabFrame', 'Standard', 'Events', '"cos#Theta_{lab}"', '1', '1', '0', '0', '0', '30e3', '-1.0', '.99999', '0', '0'],
    ['Hyp_ToppT', 'Standard', 'Topquarks', 'p_{T}^{t}', '1', '1', '0', '0', '0', '140e3', '0', '549', '0', '0'],

]



# from AnalyzerSpinCorr

for era in eras:
    print(era)
    
    for channel in channels:
        print(channel)

        runNominalfile = open("scripts/LoadPlotter_Nominal_" + era + ".sh","w") 
        runSystsfile = open("scripts/LoadPlotter_Systs_" + era + ".sh","w") 
        runEnvelopesfile = open("scripts/LoadPlotter_Envelopes_" + era + ".sh","w") 
        runBandsfile = open("scripts/LoadPlotter_Bands_" + era + ".sh","w") 

        file = open("data/binnings/HistoList_CP_" + era + "_" + channel,"w") 

        file.write('# Name, Extra, axis labels (y,x), rebin, do_dyscale, logx, logy, ymin, ymax, xmin, xmax, nbins, xbins, bcs\n')
        file.write('\n')

        for cp in controlplotlist:
            #print(cp)

#            for var in varsdict.keys():
#                for i in range(len(varsdict[var]) - 1):
#                    rootfilepath = 'Plots_'+era+'/Nominal/' + channel + '/'+cp[0]+'_step8_'+varsdict[var][i]+'to'+varsdict[var][i+1]+'_'+var+'_source.root'
#                    if os.path.exists(rootfilepath):
#                        rootfile = ROOT.TFile.Open(rootfilepath, "read")
#                       hist = rootfile.Get(cp[0]+'_step8_'+varsdict[var][i]+'to'+varsdict[var][i+1]+'_'+var+'_data')
#                       #print(varsdict[var][i]+'to'+varsdict[var][i+1])
#                       ymax = 0
#                       for j in range(1,hist.GetNbinsX()-1):
#                           if hist.GetBinContent(j) > ymax:
#                               ymax = hist.GetBinContent(j)
#                               #print(ymax)
#
#                   else:
#                       ymax = float(cp[9])
#                       #print(ymax)
#
#                   if float(cp[7]) == 0:
#                       ymax = roundup(ymax * 1.5,5000)
#                       #print(ymax)
#                   elif float(cp[7]) == 1:
#                       ymax = math.pow(10,3+roundup(math.log10(ymax * 1.5)))
#                       #print(ymax)
#
#                   file.write(cp[0]+'_step8_'+varsdict[var][i]+'to'+varsdict[var][i+1]+'_'+var+' '+cp[1]+' '+cp[2]+' '+cp[3]+' '+cp[4]+' '+cp[5]+' '+cp[6]+' '+cp[7]+' '+cp[8]+' '+str(ymax)+' '+cp[10]+' '+cp[11]+' '+cp[12]+' '+cp[13]+'\n')
#
#                   if(channel == "combined"):
#                       runNominalfile.write('nohup ./install/bin/load_Plotter -y '+era+' -t cp -s Nominal -p +'+cp[0]+'_step8_'+varsdict[var][i]+'to'+varsdict[var][i+1]+'_'+var+' &> plotslogs/'+cp[0]+'_step8_'+varsdict[var][i]+'to'+varsdict[var][i+1]+'_'+var+'_'+era+'.out \n')
#                       runSystsfile.write('nohup ./install/bin/load_Plotter -y '+era+' -t cp -s all -p +'+cp[0]+'_step8_'+varsdict[var][i]+'to'+varsdict[var][i+1]+'_'+var+' &> plotslogs/'+cp[0]+'_step8_'+varsdict[var][i]+'to'+varsdict[var][i+1]+'_'+var+'_'+era+'.out \n')
#                       runEnvelopesfile.write('nohup ./install/bin/load_Plotter -y '+era+' -t cp -m -s all -p +'+cp[0]+'_step8_'+varsdict[var][i]+'to'+varsdict[var][i+1]+'_'+var+' &> plotslogs/'+cp[0]+'_step8_'+varsdict[var][i]+'to'+varsdict[var][i+1]+'_'+var+'_'+era+'.out \n')
#                       runBandsfile.write('nohup ./install/bin/load_Plotter -y '+era+' -t cp -b -s Nominal -p +'+cp[0]+'_step8_'+varsdict[var][i]+'to'+varsdict[var][i+1]+'_'+var+' &> plotslogs/'+cp[0]+'_step8_'+varsdict[var][i]+'to'+varsdict[var][i+1]+'_'+var+'_'+era+'.out \n')
#                   else: continue

            rootfilepath = 'Plots_'+era+'/Nominal/' + channel + '/'+cp[0]+'_step8_source.root'
            if os.path.exists(rootfilepath):
                rootfile = ROOT.TFile.Open(rootfilepath, "read")
                hist = rootfile.Get(cp[0]+'_step8_data')
                ymax = 0
                for j in range(1,hist.GetNbinsX()-1):
                    if hist.GetBinContent(j) > ymax:
                        ymax = hist.GetBinContent(j)

            else:
                ymax = float(cp[9])
                #print(ymax)

            if float(cp[7]) == 0:
                ymax = roundup(ymax * 1.5,5000)
            elif float(cp[7]) == 1 and ymax != 0:
                ymax = math.pow(10,3+roundup(math.log10(ymax * 1.5)))
            file.write(cp[0]+'_step8 '+cp[1]+' '+cp[2]+' '+cp[3]+' '+cp[4]+' '+cp[5]+' '+cp[6]+' '+cp[7]+' '+cp[8]+' '+str(ymax)+' '+cp[10]+' '+cp[11]+' '+cp[12]+' '+cp[13]+'\n')

            if(channel == "combined"):
                runNominalfile.write('nohup ./install/bin/load_Plotter -y '+era+' -t cp -s Nominal -p +'+cp[0]+'_step8'+' &> plotslogs/'+cp[0]+'_step8'+'_'+era+'.out \n')
                runSystsfile.write('nohup ./install/bin/load_Plotter -y '+era+' -t cp -s all -p +'+cp[0]+'_step8'+' &> plotslogs/'+cp[0]+'_step8'+'_'+era+'.out \n')
                runEnvelopesfile.write('nohup ./install/bin/load_Plotter -y '+era+' -t cp -m -s all -p +'+cp[0]+'_step8'+' &> plotslogs/'+cp[0]+'_step8'+'_'+era+'.out \n')
                runBandsfile.write('nohup ./install/bin/load_Plotter -y '+era+' -t cp -b -s Nominal -p +'+cp[0]+'_step8'+' &> plotslogs/'+cp[0]+'_step8'+'_'+era+'.out \n')
            else: continue
                
        file.write("\n")

        file.close() 
        runNominalfile.close()
        runSystsfile.close()
        runEnvelopesfile.close()
        runBandsfile.close()




# from minitrees

for era in eras:
    print(era)
    
    for channel in channels:
        print(channel)

        runNominalfile = open("scripts/LoadPlotter_Nominal_" + era + ".sh","a") 
        runSystsfile = open("scripts/LoadPlotter_Systs_" + era + ".sh","a") 
        runEnvelopesfile = open("scripts/LoadPlotter_Envelopes_" + era + ".sh","a") 
        runBandsfile = open("scripts/LoadPlotter_Bands_" + era + ".sh","a") 

        file = open("data/binnings/HistoList_CP_" + era + "_" + channel,"a") 

        file.write('# Name, Extra, axis labels (y,x), rebin, do_dyscale, logx, logy, ymin, ymax, xmin, xmax, nbins, xbins, bcs\n')
        file.write('\n')

        for cp in controlplotlist_fromtrees:
            #print(cp)

            for var in vars2Dlist:

                rootfilepath = 'Plots_'+era+'/Nominal/' + channel + '/'+cp[0]+'_vs_'+var+'_source.root'
                if os.path.exists(rootfilepath):
                    rootfile = ROOT.TFile.Open(rootfilepath, "read")
                    hist = rootfile.Get(cp[0]+'_vs_'+var+'_data')
                    #print(varsdict[var][i]+'to'+varsdict[var][i+1])
                    ymax = 0
                    for j in range(1,hist.GetNbinsX()-1):
                        if hist.GetBinContent(j) > ymax:
                            ymax = hist.GetBinContent(j)
                            #print(ymax)

                else:
                    ymax = float(cp[9])
                    #print(ymax)

                if float(cp[7]) == 0:
                    ymax = roundup(ymax * 1.5,5000)
                    #print(ymax)
                elif float(cp[7]) == 1:
                    ymax = math.pow(10,3+roundup(math.log10(ymax * 1.5)))
                    #print(ymax)

                file.write(cp[0]+'_vs_'+var+' '+cp[1]+' '+cp[2]+' '+cp[3]+' '+cp[4]+' '+cp[5]+' '+cp[6]+' '+cp[7]+' '+cp[8]+' '+str(ymax)+' '+cp[10]+' '+cp[11]+' '+cp[12]+' '+cp[13]+'\n')

                if(channel == "combined"):
                    runNominalfile.write('nohup ./install/bin/load_Plotter -y '+era+' -t cp -s Nominal -p +'+cp[0]+'_vs_'+var+' &> plotslogs/'+cp[0]+'_vs_'+var+'_'+era+'.out \n')
                    runSystsfile.write('nohup ./install/bin/load_Plotter -y '+era+' -t cp -s all -p +'+cp[0]+'_vs_'+var+' &> plotslogs/'+cp[0]+'_vs_'+var+'_'+era+'.out \n')
                    runEnvelopesfile.write('nohup ./install/bin/load_Plotter -y '+era+' -t cp -m -s all -p +'+cp[0]+'_vs_'+var+' &> plotslogs/'+cp[0]+'_vs_'+var+'_'+era+'.out \n')
                    runBandsfile.write('nohup ./install/bin/load_Plotter -y '+era+' -t cp -b -s Nominal -p +'+cp[0]+'_vs_'+var+' &> plotslogs/'+cp[0]+'_vs_'+var+'_'+era+'.out \n')
                else: continue

            rootfilepath = 'Plots_'+era+'/Nominal/' + channel + '/'+cp[0]+'_source.root'
            if os.path.exists(rootfilepath):
                rootfile = ROOT.TFile.Open(rootfilepath, "read")
                hist = rootfile.Get(cp[0]+'_data')
                ymax = 0
                for j in range(1,hist.GetNbinsX()-1):
                    if hist.GetBinContent(j) > ymax:
                        ymax = hist.GetBinContent(j)

            else:
                ymax = float(cp[9])
                #print(ymax)

            if float(cp[7]) == 0:
                ymax = roundup(ymax * 1.5,5000)
            elif float(cp[7]) == 1 and ymax != 0:
                ymax = math.pow(10,3+roundup(math.log10(ymax * 1.5)))
            file.write(cp[0]+' '+cp[1]+' '+cp[2]+' '+cp[3]+' '+cp[4]+' '+cp[5]+' '+cp[6]+' '+cp[7]+' '+cp[8]+' '+str(ymax)+' '+cp[10]+' '+cp[11]+' '+cp[12]+' '+cp[13]+'\n')

            if(channel == "combined"):
                runNominalfile.write('nohup ./install/bin/load_Plotter -y '+era+' -t cp -s Nominal -p +'+cp[0]+''+' &> plotslogs/'+cp[0]+''+'_'+era+'.out \n')
                runSystsfile.write('nohup ./install/bin/load_Plotter -y '+era+' -t cp -s all -p +'+cp[0]+''+' &> plotslogs/'+cp[0]+''+'_'+era+'.out \n')
                runEnvelopesfile.write('nohup ./install/bin/load_Plotter -y '+era+' -t cp -m -s all -p +'+cp[0]+''+' &> plotslogs/'+cp[0]+''+'_'+era+'.out \n')
                runBandsfile.write('nohup ./install/bin/load_Plotter -y '+era+' -t cp -b -s Nominal -p +'+cp[0]+''+' &> plotslogs/'+cp[0]+''+'_'+era+'.out \n')
            else: continue

        file.write("\n")

        file.close() 
        runNominalfile.close()
        runSystsfile.close()
        runEnvelopesfile.close()
        runBandsfile.close()


controlplototherlist = [

#['#', 'Name,', 'Extra,', 'axis', 'labels', '(y,x),', 'rebin,', 'do_dyscale,', 'logx,', 'logy,', 'ymin,', 'ymax,', 'xmin,', 'xmax,', 'nbins,', 'xbins,', 'bcs'],


# Pile-Up Reweighting (Step 3 and Step 8)
['vertMulti_step3', 'Standard', 'Events', 'N_{vtx}', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0'],
['vertMulti_noPU_step3', 'Standard', 'Events', 'N_{vtx}', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0'],
['vertMulti_step8', 'Standard', 'Events', 'N_{vtx}', '5', '1', '0', '0', '0', '0', '0', '0', '0', '0'],
['vertMulti_noPU_step8', 'Standard', 'Events', 'N_{vtx}', '5', '1', '0', '0', '0', '0', '0', '0', '0', '0'],

# Full Dilepton System / pre-Zpeak cut (Step 3)
['DIMFull', 'Standard', 'Events', 'm_{l,#bar{l}}', '2', '0', '0', '1', '800', '0', '0', '0', '0', '0'],

# Dilepton System / pre-jetMult cut (Step 4)
['HypjetMulti_diLep', 'Standard', 'Events', 'N_{jets}', '1', '0', '0', '1', '11', '1e10', '0', '6', '0', '0'],
['MuonpT_diLep', 'Standard', 'Muons', 'p_{T,#mu}', '2', '0', '0', '1', '16', '0', '20', '399', '0', '0'],
['ElectronpT_diLep', 'Standard', 'Electrons', 'p_{T,e}', '2', '0', '0', '1', '22', '0', '20', '399', '0', '0'],
['MET_diLep', 'Standard', 'Events', 'MET', '10', '0', '0', '0', '0', '0', '0', '299', '0', '0'],

# pre-MET cut (Step 5)
['MET_preMETcut', 'Standard', 'Events', 'MET', '10', '0', '0', '0', '0', '0', '0', '299', '0', '0'],

# pre-btag cut (Step 6)
['MuonpT_postMETcut', 'Standard', 'Muons', 'p_{T,#mu}', '2', '1', '0', '1', '11', '0', '20', '399', '0', '0'],
['ElectronpT_postMETcut', 'Standard', 'Electrons', 'p_{T,e}', '2', '1', '0', '1', '1', '0', '20', '399', '0', '0'],
['MuonEta_postMETcut', 'Standard', 'Muons', '#eta_{#mu}', '2', '1', '0', '0', '0', '30e3', '-2.4', '2.39', '0', '0'],
['ElectronEta_postMETcut', 'Standard', 'Electrons', '#eta_{e}', '2', '1', '0', '0', '0', '30e3', '-2.4', '2.39', '0', '0'],
['LeptonpT_postMETcut', 'Standard', 'Leptons', 'p_{T}^{l}', '2', '1', '0', '1', '11', '0', '20', '399', '0', '0'],
['AntiLeptonpT_postMETcut', 'Standard', 'Leptons', 'p_{T}^{l}', '2', '1', '0', '1', '11', '0', '20', '399', '0', '0'],
['LeptonEta_postMETcut', 'Standard', 'Leptons', '#eta_{l}', '2', '1', '0', '0', '0', '25e3', '-2.4', '2.39', '0', '0'],
['AntiLeptonEta_postMETcut', 'Standard', 'Leptons', '#eta_{l}', '2', '1', '0', '0', '0', '25e3', '-2.4', '2.39', '0', '0'],
['HypjetMulti_noBTag', 'Standard', 'Events', 'N_{jets}', '1', '1', '0', '1', '4', '1e8', '2', '6', '0', '0'],
['HypBjetMulti_noBTag', 'Standard', 'Events', 'N_{bjets}', '1', '1', '0', '1', '0.6', '1e9', '0', '5', '0', '0'],


# pre-KinReco (Step 7)
['MuonpT_BeforeKinReco', 'Standard', 'Muons', 'p_{T,#mu}', '2', '1', '0', '1', '11', '0', '20', '399', '0', '0'],
['ElectronpT_BeforeKinReco', 'Standard', 'Electrons', 'p_{T,e}', '2', '1', '0', '1', '1', '0', '20', '399', '0', '0'],
['MuonEta_BeforeKinReco', 'Standard', 'Muons', '#eta_{#mu}', '2', '1', '0', '0', '0', '30e3', '-2.4', '2.39', '0', '0'],
['ElectronEta_BeforeKinReco', 'Standard', 'Electrons', '#eta_{e}', '2', '1', '0', '0', '0', '30e3', '-2.4', '2.39', '0', '0'],
['LeptonpT_AfterBTag', 'Standard', 'Leptons', 'p_{T}^{l}', '1', '1', '0', '0', '0', '0', '20', '159', '0', '0'],
['LeptonpT_BeforeKinReco', 'Standard', 'Leptons', 'p_{T}^{l}', '1', '1', '0', '0', '0', '0', '20', '159', '0', '0'],
['LeptonEta_BeforeKinReco', 'Standard', 'Leptons', '#eta_{l}', '2', '1', '0', '0', '0', '3.5e5', '-2.4', '2.39', '0', '0'],
['HypjetMulti_step7', 'Standard', 'Events', 'N_{jets}', '1', '1', '0', '1', '3', '1e8', '2', '6', '0', '0'],
['jetHT', 'Standard', 'Jets', 'H_{T}', '2', '1', '0', '0', '0', '0', '0', '799', '0', '0'],
['jetpT', 'Standard', 'Jets', 'p_{T}', '2', '1', '0', '1', '90', '0', '30', '399', '0', '0'],
['MET_step7', 'Standard', 'Events', 'MET', '10', '1', '0', '0', '0', '0', '0', '299', '0', '0'],
['HypBjetMulti_step7', 'Standard', 'Events', 'N_{bjets}', '1', '1', '0', '1', '0.6', '1e9', '1', '5', '0', '0'],
['BjetpT_BeforeKinReco', 'Standard', 'b-jets', 'p_{T}^{b}', '2', '1', '0', '1', '90', '0', '30', '399', '0', '0'],

# post-KinReco (Step 8)

['DIMFull_fullSel', 'Standard', 'Events', 'm_{l,#bar{l}}', '2', '1', '0', '0', '0', '0', '0', '0', '0', '0'],
['HypJetMultpt30', 'Standard', 'Events', 'Jets', '1', '1', '0', '1', '3', '1e8', '2', '9', '0', '0'],
['HypToppTLead', 'Standard', 'Topquarks', 'p_{T}^{leadingt}', '20', '1', '0', '0', '0', '35e3', '0', '559', '0', '0'],
['HypToppTNLead', 'Standard', 'Topquarks', 'p_{T}^{trailingt}', '20', '1', '0', '0', '0', '42e3', '0', '559', '0', '0'],
    
# Leptons
['HypLeptonpT', 'Standard', 'Events', 'p_{T}^{l}', '10', '1', '0', '0', '0', '0', '20', '159', '0', '0'],
['HypAntiLeptonpT', 'Standard', 'Events', 'p_{T}^{#bar{l}}', '10', '1', '0', '0', '0', '0', '20', '159', '0', '0'],
['HypLeptonEta', 'Standard', 'Events', '#eta^{l}', '10', '1', '0', '0', '0', '0', '-3', '3', '0', '0'],
['HypAntiLeptonEta', 'Standard', 'Events', '#eta^{#bar{l}}', '10', '1', '0', '0', '0', '0', '-3', '3', '0', '0'],

#['HypLeptonpT_AntiTopFrame', 'Standard', 'Events', 'p_{T}^{l}', '10', '1', '0', '0', '0', '0', '20', '159', '0', '0'],
#['HypAntiLeptonpT_TopFrame', 'Standard', 'Events', 'p_{T}^{#bar{l}}', '10', '1', '0', '0', '0', '0', '20', '159', '0', '0'],
#['HypLeptonEta_AntiTopFrame', 'Standard', 'Events', '#eta^{l}', '10', '1', '0', '0', '0', '0', '-3', '3', '0', '0'],
#['HypAntiLeptonEta_TopFrame', 'Standard', 'Events', '#eta^{#bar{l}}', '10', '1', '0', '0', '0', '0', '-3', '3', '0', '0'],

['HypLLBarMass', 'Standard', 'Events', 'm^{l^{+}l^{-}}', '10', '1', '0', '0', '0', '0', '0', '500', '0', '0'],
['HypLLBarpT', 'Standard', 'Events', 'p_{T}^{l^{+}l^{-}}', '4', '1', '0', '0', '0', '0', '0', '500', '0', '0'],
['HypLLBarDPhi', 'Standard', 'Events', '"|#Delta#phi_{l#bar{l}}|"', '157', '1', '0', '0', '0', '30e3', '0', '3.13999', '0', '0'],
#['HypLLBarDEta', 'Standard', 'Events', '"#Delta#eta{l,#bar{l}}"', '6', '1', '0', '0', '0', '0', '-3.0', '3.0', '0', '0'],

# Jets 

['HypJetMulti', 'Standard', 'Events', 'N_{jets}', '1', '1', '0', '1', '3', '1e8', '2', '9', '0', '0'],
#['HypJetpT', 'Standard', 'Events', 'p_{T}', '2', '1', '0', '0', '0', '0', '30', '399', '0', '0'],
#['HypJetHT', 'Standard', 'Events', 'H_{T}', '2', '1', '0', '0', '0', '0', '0', '799', '0', '0'],
['HypexjetMulti', 'Standard', 'Events', 'N_{extra#kern[-0.8]jets}', '1', '1', '0', '1', '3', '1e8', '0', '5', '0', '0'],

# Extra Jets 

#['HypExtraJetMulti', 'Standard', 'Events', 'N_{jets}', '1', '1', '0', '1', '3', '1e8', '2', '9', '0', '0'],
#['HypExtraJetpT', 'Standard', 'Events', 'p_{T}', '2', '1', '0', '0', '0', '0', '30', '399', '0', '0'],
#['HypExtraJetHT', 'Standard', 'Events', 'H_{T}', '2', '1', '0', '0', '0', '0', '0', '799', '0', '0'],

# B Jets

#['HypBjetMulti', 'Standard', 'Events', 'N_{bjets}', '1', '1', '0', '1', '0.6', '1e9', '1', '6', '0', '0'],
['HypBJetpT', 'Standard', 'Events', 'p_{T}^{b}', '10', '1', '0', '0', '0', '0', '39', '399', '0', '0'],
['HypBJetEta', 'Standard', 'Events', '#eta_{b}', '4', '1', '0', '0', '0', '63e3', '-2.4', '2.39', '0', '0'],
['HypAntiBJetpT', 'Standard', 'Events', 'p_{T}^{#bar{b}}', '10', '1', '0', '0', '0', '0', '39', '399', '0', '0'],
['HypAntiBJetEta', 'Standard', 'Events', '#eta_{#bar{b}}', '4', '1', '0', '0', '0', '63e3', '-2.4', '2.39', '0', '0'],

#['HypBBBarpT', 'Standard', 'Events', 'p_{T}^{b#bar{b}}', '10', '1', '0', '0', '0', '0', '0', '350', '0', '0'],
#['HypBBBarMass', 'Standard', 'Events', 'm^{b#bar{b}}', '10', '1', '0', '0', '0', '0', '0', '600', '0', '0'],
#['HypBBBarDPhi', 'Standard', 'Events', '|#Delta#phi_{b,#bar{b}}|', '10', '1', '0', '0', '0', '27999', '0.0', '3.19', '0', '0'],

['HypLeptonBjetMass', 'Standard', 'Events', 'm_{l#bar{b}}', '5', '1', '0', '0', '0', '130e3', '0', '179', '0', '0'],
['HypAntiLeptonBjetMass', 'Standard', 'Events', 'm_{#bar{l}b}', '5', '1', '0', '0', '0', '130e3', '0', '179', '0', '0'],

# MET and Neutrinos

#['HypMET', 'Standard', 'Events', 'MET', '10', '0', '0', '0', '0', '0', '0', '300', '0', '0'],

# Top and Anti-Top

['HypToppT', 'Standard', 'Topquarks', 'p_{T}^{t}', '20', '1', '0', '0', '0', '140e3', '0', '559', '0', '0'],
['HypTopEta', 'Standard', 'Topquarks', '#eta_{t}', '8', '1', '0', '0', '0', '150e3', '-2.4', '2.39', '0', '0'],
['HypTopMass', 'Standard', 'Topquarks', 'm_{t}', '1', '1', '0', '0', '0', '0', '170', '175', '0', '0'],
['HypTopRapidity', 'Standard', 'Topquarks', 'y_{t}', '8', '1', '0', '0', '0', '150e3', '-2.6', '2.59', '0', '0'],
['HypTopRapidityAbs', 'Standard', 'Topquarks', '|y_{t}|', '4', '1', '0', '0', '0', '0', '0', '2.59', '0', '0'],
#['HypToppT_TTBarFrame', 'Standard', 'Topquarks', 'p_{T}^{t}', '20', '1', '0', '0', '0', '140e3', '0', '559', '0', '0'],
#['HypTopEta_TTBarFrame', 'Standard', 'Topquarks', '#eta_{t}', '8', '1', '0', '0', '0', '150e3', '-2.4', '2.39', '0', '0'],
#['HypTopRapidity_TTBarFrame', 'Standard', 'Topquarks', 'y_{t}', '8', '1', '0', '0', '0', '150e3', '-2.6', '2.59', '0', '0'],

['HypAntiToppT', 'Standard', 'Topquarks', 'p_{T}^{t}', '20', '1', '0', '0', '0', '140e3', '0', '559', '0', '0'],
['HypAntiTopEta', 'Standard', 'Topquarks', '#eta_{t}', '8', '1', '0', '0', '0', '150e3', '-2.4', '2.39', '0', '0'],
['HypAntiTopMass', 'Standard', 'Topquarks', 'm_{t}', '1', '1', '0', '0', '0', '0', '170', '175', '0', '0'],
['HypAntiTopRapidity', 'Standard', 'Topquarks', 'y_{t}', '8', '1', '0', '0', '0', '150e3', '-2.6', '2.59', '0', '0'],
['HypAntiTopRapidityAbs', 'Standard', 'Topquarks', '|y_{t}|', '4', '1', '0', '0', '0', '0', '0', '2.59', '0', '0'],
#['HypAntiToppT_TTBarFrame', 'Standard', 'Topquarks', 'p_{T}^{t}', '20', '1', '0', '0', '0', '140e3', '0', '559', '0', '0'],
#['HypAntiTopEta_TTBarFrame', 'Standard', 'Topquarks', '#eta_{t}', '8', '1', '0', '0', '0', '150e3', '-2.4', '2.39', '0', '0'],
#['HypAntiTopRapidity_TTBarFrame', 'Standard', 'Topquarks', 'y_{t}', '8', '1', '0', '0', '0', '150e3', '-2.6', '2.59', '0', '0'],

['HypScatteringAngle_TTBarFrame', 'Standard', 'Events', '"cos#Theta"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99999', '0', '0'],
['HypScatteringAngle_LabFrame', 'Standard', 'Events', '"cos#Theta_{lab}"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99999', '0', '0'],

['HypTTBarMass', 'Standard', 'Toppairs', 'm_{t#bar{t}}', '60', '1', '0', '0', '0', '150e3', '300', '1999', '0', '0'],
['HypTTBarRapidity', 'Standard', 'Toppairs', 'y_{t#bar{t}}', '8', '1', '0', '0', '0', '100e3', '-2.6', '2.59', '0', '0'],
['HypTTBarRapidityAbs', 'Standard', 'Toppairs', '|y_{t#bar{t}}|', '4', '1', '0', '0', '0', '50e3', '0', '2.59', '0', '0'],
['HypTTBarpT', 'Standard', 'Toppairs', 'p_{T}^{t#bar{t}}', '20', '1', '0', '0', '0', '1.5e5', '0', '499', '0', '0'],
#['HypTTBarDeltaPhi', 'Standard', 'N_{Events}', '#Delta#phi_{t,#bar{t}}', '10', '1', '0', '0', '0', '0', '0', '3.15', '0', '0'],
#['HypTTBarDeltaRapidity', 'Standard', 'N_{Events}', '|y_{t}|-|y_{#bar{t}}|', '8', '1', '0', '0', '0', '70e3', '-2.6', '2.59', '0', '0'],

#['HypScatteringAngle_TTBarFrame', 'Standard', 'Events', '"cos#Theta"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99999', '0', '0'],
#['HypScatteringAngle_LabFrame', 'Standard', 'Events', '"cos#Theta_{lab}"', '10', '1', '0', '0', '0', '30e3', '-1.0', '.99999', '0', '0'],

#['HypAbsDeltaPhiExtraJet12', 'Standard', 'Events', '|#Delta#phi_{jet1,jet2}|', '1', '1', '0', '0', '0', '10000', '0.0', '3.15', '0', '0'],
#['HypAbsDeltaPhiExtraJet12_eta1', 'Standard', 'Events', '|#Delta#phi_{jet1,jet2}|', '1', '1', '0', '0', '0', '1500', '0.0', '3.15', '0', '0'],
#['HypAbsDeltaPhiExtraJet12_eta2', 'Standard', 'Events', '|#Delta#phi_{jet1,jet2}|', '1', '1', '0', '0', '0', '1000', '0.0', '3.15', '0', '0'],
#['HypAbsDeltaPhiExtraJet12_eta3', 'Standard', 'Events', '|#Delta#phi_{jet1,jet2}|', '1', '1', '0', '0', '0', '500', '0.0', '3.15', '0', '0'],
#['HypAbsDeltaPhiExtraJet12_eta4', 'Standard', 'Events', '|#Delta#phi_{jet1,jet2}|', '1', '1', '0', '0', '0', '100', '0.0', '3.15', '0', '0'],
#['HypjetMultiXSec', 'Standard', 'N_{Events}', 'N_{jets}', '1', '1', '0', '0', '0', '0', '2', '9', '0', '0'],

#['HypBJetpTLead', 'Standard', 'b-jets', 'p_{T}(lead.b-jet)', '10', '1', '0', '0', '0', '0', '20', '299', '0', '0'],
#['HypBJetpTNLead', 'Standard', 'b-jets', 'p_{T}(2^{nd}lead.b-jet)', '10', '1', '0', '0', '0', '0', '20', '299', '0', '0'],
#['HypExtraJetpT', 'Standard', '1^{st}-add-jet', 'p_{T}', '10', '1', '0', '0', '0', '0', '20', '499', '0', '0'],
#['HypExtraJetpT2', 'Standard', '2^{nd}-add-jet', 'p_{T}', '10', '1', '0', '0', '0', '0', '20', '499', '0', '0'],
#['HypExtraJetpT3', 'Standard', '3^{th}-add-jet', 'p_{T}', '10', '1', '0', '0', '0', '0', '20', '499', '0', '0'],
#['HypExtraJetpT4', 'Standard', '4^{th}-add-jet', 'p_{T}', '10', '1', '0', '0', '0', '0', '20', '499', '0', '0'],

#['#', 'Name,', 'Extra,', 'axis', 'labels', '(y,x),', 'rebin,', 'do_dyscale,', 'logx,', 'logy,', 'ymin,', 'ymax,', 'xmin,', 'xmax,', 'nbins,', 'xbins,', 'bcs'],

['events_weighted_step3', 'Standard', 'Events', '-', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0'],
['events_weighted_step4', 'Standard', 'Events', '-', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0'],
['events_weighted_step5', 'Standard', 'Events', '-', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0'],
['events_weighted_step6', 'Standard', 'Events', '-', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0'],
['events_weighted_step7', 'Standard', 'Events', '-', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0'],
['events_weighted_step7L', 'Standard', 'Events', '-', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0'],
['events_weighted_step8', 'Standard', 'Events', '-', '1', '1', '0', '0', '0', '0', '0', '0', '0', '0'],

#['#', 'Name,', 'Extra,', 'axis', 'labels', '(y,x),', 'rebin,', 'do_dyscale,', 'logx,', 'logy,', 'ymin,', 'ymax,', 'xmin,', 'xmax,', 'nbins,', 'xbins,', 'bcs'],

['dyScaling_Allh1_step4', 'Standard', 'Entries', 'm^{ll}', '1', '0', '0', '1', '10', '0', '0', '0', '0', '0'],
['dyScaling_Allh1_step5', 'Standard', 'Entries', 'm^{ll}', '1', '0', '0', '1', '10', '0', '0', '0', '0', '0'],
['dyScaling_Allh1_step6', 'Standard', 'Entries', 'm^{ll}', '1', '0', '0', '1', '10', '0', '0', '0', '0', '0'],
['dyScaling_Allh1_step7', 'Standard', 'Entries', 'm^{ll}', '1', '0', '0', '1', '10', '0', '0', '0', '0', '0'],
['dyScaling_Allh1_step8', 'Standard', 'Entries', 'm^{ll}', '1', '0', '0', '1', '10', '0', '0', '0', '0', '0'],

['dyScaling_TTh1_step4', 'Standard', 'Entries', 'm^{ll}', '1', '0', '0', '1', '10', '0', '0', '0', '0', '0'],
['dyScaling_TTh1_step5', 'Standard', 'Entries', 'm^{ll}', '1', '0', '0', '1', '10', '0', '0', '0', '0', '0'],
['dyScaling_TTh1_step6', 'Standard', 'Entries', 'm^{ll}', '1', '0', '0', '1', '10', '0', '0', '0', '0', '0'],
['dyScaling_TTh1_step7', 'Standard', 'Entries', 'm^{ll}', '1', '0', '0', '1', '10', '0', '0', '0', '0', '0'],
['dyScaling_TTh1_step8', 'Standard', 'Entries', 'm^{ll}', '1', '0', '0', '1', '10', '0', '0', '0', '0', '0'],

['dyScaling_Zh1_step4zWindow', 'Standard', 'Entries', 'm^{ll}', '1', '0', '0', '1', '10', '0', '60', '120', '0', '0'],
['dyScaling_Zh1_step5zWindow', 'Standard', 'Entries', 'm^{ll}', '1', '0', '0', '1', '10', '0', '60', '120', '0', '0'],
['dyScaling_Zh1_step6zWindow', 'Standard', 'Entries', 'm^{ll}', '1', '0', '0', '1', '10', '0', '60', '120', '0', '0'],
['dyScaling_Zh1_step7zWindow', 'Standard', 'Entries', 'm^{ll}', '1', '0', '0', '1', '10', '0', '60', '120', '0', '0'],
['dyScaling_Zh1_step8zWindow', 'Standard', 'Entries', 'm^{ll}', '1', '0', '0', '1', '10', '0', '60', '120', '0', '0'],


#['#', 'Name,', 'Extra,', 'axis', 'labels', '(y,x),', 'rebin,', 'do_dyscale,', 'logx,', 'logy,', 'ymin,', 'ymax,', 'xmin,', 'xmax,', 'nbins,', 'xbins,', 'bcs'],

#['basic_vertMulti_step3', 'Standard', 'Events', 'N_{vtx}', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0'],
#['basic_vertMulti_noPU_step3', 'Standard', 'Events', 'N_{vtx}', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0'],
#['basic_vertMulti_step4', 'Standard', 'Events', 'N_{vtx}', '2', '0', '0', '0', '0', '0', '0', '0', '0', '0'],
#['basic_vertMulti_noPU_step4', 'Standard', 'Events', 'N_{vtx}', '2', '0', '0', '0', '0', '0', '0', '0', '0', '0'],
#['basic_vertMulti_step5', 'Standard', 'Events', 'N_{vtx}', '2', '1', '0', '0', '0', '0', '0', '0', '0', '0'],
#['basic_vertMulti_noPU_step5', 'Standard', 'Events', 'N_{vtx}', '2', '1', '0', '0', '0', '0', '0', '0', '0', '0'],
#['basic_vertMulti_step6', 'Standard', 'Events', 'N_{vtx}', '2', '1', '0', '0', '0', '0', '0', '0', '0', '0'],
#['basic_vertMulti_noPU_step6', 'Standard', 'Events', 'N_{vtx}', '2', '1', '0', '0', '0', '0', '0', '0', '0', '0'],
#['basic_vertMulti_step7', 'Standard', 'Events', 'N_{vtx}', '5', '1', '0', '0', '0', '0', '0', '0', '0', '0'],
#['basic_vertMulti_noPU_step7', 'Standard', 'Events', 'N_{vtx}', '5', '1', '0', '0', '0', '0', '0', '0', '0', '0'],
#['basic_vertMulti_step8', 'Standard', 'Events', 'N_{vtx}', '5', '1', '0', '0', '0', '0', '0', '0', '0', '0'],
#['basic_vertMulti_noPU_step8', 'Standard', 'Events', 'N_{vtx}', '5', '1', '0', '0', '0', '0', '0', '0', '0', '0'],

#['#', 'Name,', 'Extra,', 'axis', 'labels', '(y,x),', 'rebin,', 'do_dyscale,', 'logx,', 'logy,', 'ymin,', 'ymax,', 'xmin,', 'xmax,', 'nbins,', 'xbins,', 'bcs'],
#['###', 'All', 'selected', 'lepton', 'properties']
#['basic_dilepton_mass_step3', 'Standard', 'Entries', 'm_{ll}', '1', '0', '0', '1', '0.1', '0', '20', '299', '0', '0'],
#['basic_dilepton_mass_step4', 'Standard', 'Entries', 'm_{ll}', '1', '0', '0', '1', '0.1', '0', '20', '299', '0', '0'],
#['basic_dilepton_mass_step4zWindow', 'Standard', 'Entries', 'm_{ll}', '1', '0', '0', '1', '0.1', '0', '40', '159', '0', '0'],
#['basic_dilepton_mass_step5', 'Standard', 'Entries', 'm_{ll}', '1', '1', '0', '1', '0.1', '0', '20', '299', '0', '0'],
#['basic_dilepton_mass_step5zWindow', 'Standard', 'Entries', 'm_{ll}', '1', '0', '0', '1', '0.1', '0', '40', '159', '0', '0'],
#['basic_dilepton_mass_step6', 'Standard', 'Entries', 'm_{ll}', '1', '1', '0', '1', '0.1', '0', '20', '299', '0', '0'],
#['basic_dilepton_mass_step6zWindow', 'Standard', 'Entries', 'm_{ll}', '1', '0', '0', '1', '0.1', '0', '40', '159', '0', '0'],
#['basic_dilepton_mass_step7', 'Standard', 'Entries', 'm_{ll}', '1', '1', '0', '1', '0.1', '0', '20', '299', '0', '0'],
#['basic_dilepton_mass_step7zWindow', 'Standard', 'Entries', 'm_{ll}', '1', '0', '0', '1', '0.1', '0', '40', '159', '0', '0'],
#['basic_dilepton_mass_step8', 'Standard', 'Entries', 'm_{ll}', '1', '1', '0', '1', '0.1', '0', '20', '299', '0', '0'],
#['basic_dilepton_mass_step8zWindow', 'Standard', 'Entries', 'm_{ll}', '1', '0', '0', '1', '0.1', '0', '40', '159', '0', '0'],

#['basic_dilepton_pt_step2', 'Standard', 'Entries', 'p_{T}^{ll}', '2', '0', '0', '0', '0', '0', '0', '0', '0', '0'],
#['basic_dilepton_pt_step3', 'Standard', 'Entries', 'p_{T}^{ll}', '2', '0', '0', '0', '0', '0', '0', '0', '0', '0'],
#['basic_dilepton_pt_step4', 'Standard', 'Entries', 'p_{T}^{ll}', '2', '0', '0', '0', '0', '0', '0', '0', '0', '0'],
#['basic_dilepton_pt_step5', 'Standard', 'Entries', 'p_{T}^{ll}', '2', '1', '0', '0', '0', '0', '0', '0', '0', '0'],
#['basic_dilepton_pt_step6', 'Standard', 'Entries', 'p_{T}^{ll}', '2', '1', '0', '0', '0', '0', '0', '250', '0', '0'],
#['basic_dilepton_pt_step7', 'Standard', 'Entries', 'p_{T}^{ll}', '2', '1', '0', '0', '0', '0', '0', '250', '0', '0'],
#['basic_dilepton_pt_step7L', 'Standard', 'Entries', 'p_{T}^{ll}', '2', '1', '0', '0', '0', '0', '0', '250', '0', '0'],
#['basic_dilepton_pt_step8', 'Standard', 'Entries', 'p_{T}^{ll}', '2', '1', '0', '0', '0', '0', '0', '250', '0', '0'],

#['basic_lepton_multiplicity_step2', 'Standard', 'Events', 'N_{leptons}', '1', '0', '0', '1', '0.1', '0', '2', '5', '0', '0'],
#['basic_lepton_multiplicity_step3', 'Standard', 'Events', 'N_{leptons}', '1', '0', '0', '1', '0.1', '0', '2', '5', '0', '0'],
#['basic_lepton_multiplicity_step4', 'Standard', 'Events', 'N_{leptons}', '1', '0', '0', '1', '0.1', '0', '2', '5', '0', '0'],
#['basic_lepton_multiplicity_step5', 'Standard', 'Events', 'N_{leptons}', '1', '1', '0', '1', '0.1', '0', '2', '5', '0', '0'],
#['basic_lepton_multiplicity_step6', 'Standard', 'Events', 'N_{leptons}', '1', '1', '0', '1', '0.1', '0', '2', '5', '0', '0'],
#['basic_lepton_multiplicity_step7', 'Standard', 'Events', 'N_{leptons}', '1', '1', '0', '1', '0.1', '0', '2', '5', '0', '0'],
#['basic_lepton_multiplicity_step7L', 'Standard', 'Events', 'N_{leptons}', '1', '1', '0', '1', '0.1', '0', '2', '5', '0', '0'],
#['basic_lepton_multiplicity_step8', 'Standard', 'Events', 'N_{leptons}', '1', '1', '0', '1', '0.1', '0', '2', '5', '0', '0'],

#['basic_lepton_pt_step3', 'Standard', 'Leptons', 'p_{T}^{l}', '4', '0', '0', '1', '0.1', '0', '20', '299', '0', '0'],
#['basic_lepton_pt_step4', 'Standard', 'Leptons', 'p_{T}^{l}', '4', '0', '0', '1', '0.1', '0', '20', '299', '0', '0'],
#['basic_lepton_pt_step4zWindow', 'Standard', 'Leptons', 'p_{T}^{l}', '4', '0', '0', '1', '0.1', '0', '20', '299', '0', '0'],
#['basic_lepton_pt_step5', 'Standard', 'Leptons', 'p_{T}^{l}', '4', '1', '0', '1', '0.1', '0', '20', '299', '0', '0'],
#['basic_lepton_pt_step5zWindow', 'Standard', 'Leptons', 'p_{T}^{l}', '4', '1', '0', '1', '0.1', '0', '20', '299', '0', '0'],
#['basic_lepton_pt_step6', 'Standard', 'Leptons', 'p_{T}^{l}', '4', '1', '0', '1', '0.1', '0', '20', '299', '0', '0'],
#['basic_lepton_pt_step6zWindow', 'Standard', 'Leptons', 'p_{T}^{l}', '4', '1', '0', '1', '0.1', '0', '20', '299', '0', '0'],
#['basic_lepton_pt_step7', 'Standard', 'Leptons', 'p_{T}^{l}', '3', '1', '0', '0', '0', '0', '20', '159', '0', '0'],
#['basic_lepton_pt_step7zWindow', 'Standard', 'Leptons', 'p_{T}^{l}', '3', '1', '0', '0', '0', '0', '20', '159', '0', '0'],
#['basic_lepton_pt_step8', 'Standard', 'Leptons', 'p_{T}^{l}', '4', '1', '0', '1', '0.1', '0', '20', '299', '0', '0'],
#['basic_lepton_pt_step8zWindow', 'Standard', 'Leptons', 'p_{T}^{l}', '4', '1', '0', '1', '0.1', '0', '20', '299', '0', '0'],

#['basic_lepton_eta_step3', 'Standard', 'Leptons', '#eta_{l}', '2', '0', '0', '0', '0', '0', '-2.4', '2.39', '0', '0'],
#['basic_lepton_eta_step4', 'Standard', 'Leptons', '#eta_{l}', '2', '0', '0', '0', '0', '0', '-2.4', '2.39', '0', '0'],
#['basic_lepton_eta_step4zWindow', 'Standard', 'Leptons', '#eta_{l}', '2', '0', '0', '0', '0', '0', '-2.4', '2.39', '0', '0'],
#['basic_lepton_eta_step5', 'Standard', 'Leptons', '#eta_{l}', '2', '1', '0', '0', '0', '0', '-2.4', '2.39', '0', '0'],
#['basic_lepton_eta_step5zWindow', 'Standard', 'Leptons', '#eta_{l}', '2', '1', '0', '0', '0', '0', '-2.4', '2.39', '0', '0'],
#['basic_lepton_eta_step6', 'Standard', 'Leptons', '#eta_{l}', '2', '1', '0', '0', '0', '0', '-2.4', '2.39', '0', '0'],
#['basic_lepton_eta_step6zWindow', 'Standard', 'Leptons', '#eta_{l}', '2', '1', '0', '0', '0', '0', '-2.4', '2.39', '0', '0'],
#['basic_lepton_eta_step7', 'Standard', 'Leptons', '#eta_{l}', '2', '1', '0', '0', '0', '0', '-2.4', '2.39', '0', '0'],
#['basic_lepton_eta_step7zWindow', 'Standard', 'Leptons', '#eta_{l}', '2', '1', '0', '0', '0', '0', '-2.4', '2.39', '0', '0'],
#['basic_lepton_eta_step8', 'Standard', 'Leptons', '#eta_{l}', '2', '1', '0', '0', '0', '0', '-2.4', '2.39', '0', '0'],
#['basic_lepton_eta_step8zWindow', 'Standard', 'Leptons', '#eta_{l}', '2', '1', '0', '0', '0', '0', '-2.4', '2.39', '0', '0'],

#['#', 'Name,', 'Extra,', 'axis', 'labels', '(y,x),', 'rebin,', 'do_dyscale,', 'logx,', 'logy,', 'ymin,', 'ymax,', 'xmin,', 'xmax,', 'nbins,', 'xbins,', 'bcs'],
#['###', 'All', 'selected', 'jet', 'properties']

#['basic_jet_multiplicity_step3', 'Standard', 'Events', 'N_{jets}', '1', '0', '0', '1', '1e2', '0', '0', '9', '0', '0'],
#['basic_jet_multiplicity_step4', 'Standard', 'Events', 'N_{jets}', '1', '0', '0', '1', '1e2', '0', '0', '9', '0', '0'],
#['basic_jet_multiplicity_step4zWindow', 'Standard', 'Events', 'N_{jets}', '1', '0', '0', '1', '1e2', '0', '0', '9', '0', '0'],
#['basic_jet_multiplicity_step5', 'Standard', 'Events', 'N_{jets}', '1', '1', '0', '0', '0', '0', '0', '9', '0', '0'],
#['basic_jet_multiplicity_step5zWindow', 'Standard', 'Events', 'N_{jets}', '1', '1', '0', '0', '0', '0', '0', '9', '0', '0'],
#['basic_jet_multiplicity_step6', 'Standard', 'Events', 'N_{jets}', '1', '1', '0', '1', '1e2', '0', '2', '9', '0', '0'],
#['basic_jet_multiplicity_step6zWindow', 'Standard', 'Events', 'N_{jets}', '1', '1', '0', '1', '1e2', '0', '2', '9', '0', '0'],
#['basic_jet_multiplicity_step7', 'Standard', 'Events', 'N_{jets}', '1', '1', '0', '0', '0', '0', '2', '9', '0', '0'],
#['basic_jet_multiplicity_step7zWindow', 'Standard', 'Events', 'N_{jets}', '1', '1', '0', '0', '0', '0', '2', '9', '0', '0'],
#['basic_jet_multiplicity_step8', 'Standard', 'Events', 'N_{jets}', '1', '1', '0', '1', '1e2', '0', '2', '9', '0', '0'],
#['basic_jet_multiplicity_step8zWindow', 'Standard', 'Events', 'N_{jets}', '1', '1', '0', '1', '1e2', '0', '2', '9', '0', '0'],

#['basic_jet_pt_step3', 'Standard', 'Jets', 'p_{T}^{jet}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet_pt_step3zWindow', 'Standard', 'Jets', 'p_{T}^{jet}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet_pt_step4', 'Standard', 'Jets', 'p_{T}^{jet}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet_pt_step4zWindow', 'Standard', 'Jets', 'p_{T}^{jet}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet_pt_step5', 'Standard', 'Jets', 'p_{T}^{jet}', '5', '1', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet_pt_step5zWindow', 'Standard', 'Jets', 'p_{T}^{jet}', '5', '1', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet_pt_step6', 'Standard', 'Jets', 'p_{T}^{jet}', '5', '1', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet_pt_step6zWindow', 'Standard', 'Jets', 'p_{T}^{jet}', '5', '1', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet_pt_step7', 'Standard', 'Jets', 'p_{T}^{jet}', '5', '1', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet_pt_step7zWindow', 'Standard', 'Jets', 'p_{T}^{jet}', '5', '1', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet_pt_step8', 'Standard', 'Jets', 'p_{T}^{jet}', '5', '1', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet_pt_step8zWindow', 'Standard', 'Jets', 'p_{T}^{jet}', '5', '1', '0', '1', '0.1', '0', '30', '299', '0', '0'],

#['basic_jet_eta_step4', 'Standard', 'Jets', '#eta_{jet}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet_eta_step4zWindow', 'Standard', 'Jets', '#eta_{jet}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet_eta_step5', 'Standard', 'Jets', '#eta_{jet}', '2', '1', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet_eta_step5zWindow', 'Standard', 'Jets', '#eta_{jet}', '2', '1', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet_eta_step6', 'Standard', 'Jets', '#eta_{jet}', '2', '1', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet_eta_step6zWindow', 'Standard', 'Jets', '#eta_{jet}', '2', '1', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet_eta_step7', 'Standard', 'Jets', '#eta_{jet}', '2', '1', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet_eta_step7zWindow', 'Standard', 'Jets', '#eta_{jet}', '2', '1', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet_eta_step8', 'Standard', 'Jets', '#eta_{jet}', '2', '1', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet_eta_step8zWindow', 'Standard', 'Jets', '#eta_{jet}', '2', '1', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],

#['basic_jet2lead_pt_step2', 'Standard', 'Jets', 'p_{T}^{jet}', '2', '0', '0', '1', '1', '0', '30', '299', '0', '0'],
#['basic_jet2lead_pt_step3', 'Standard', 'Jets', 'p_{T}^{jet}', '2', '0', '0', '1', '1', '0', '30', '299', '0', '0'],
#['basic_jet2lead_pt_step4', 'Standard', 'Jets', 'p_{T}^{jet}', '2', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet2lead_pt_step5', 'Standard', 'Jets', 'p_{T}^{jet}', '2', '1', '0', '0', '0', '0', '0', '349', '0', '0'],
#['basic_jet2lead_pt_step6', 'Standard', 'Jets', 'p_{T}^{jet}', '5', '1', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet2lead_pt_step7', 'Standard', 'Jets', 'p_{T}^{jet}', '5', '1', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet2lead_pt_step8', 'Standard', 'Jets', 'p_{T}^{jet}', '5', '1', '0', '1', '0.1', '0', '30', '299', '0', '0'],

#['basic_jet1_pt_step4', 'Standard', 'Jets', 'p_{T}^{jet1}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet2_pt_step4', 'Standard', 'Jets', 'p_{T}^{jet2}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet3_pt_step4', 'Standard', 'Jets', 'p_{T}^{jet3}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet4_pt_step4', 'Standard', 'Jets', 'p_{T}^{jet4}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet5_pt_step4', 'Standard', 'Jets', 'p_{T}^{jet5}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet6_pt_step4', 'Standard', 'Jets', 'p_{T}^{jet6}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet1_pt_step4zWindow', 'Standard', 'Jets', 'p_{T}^{jet1}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet2_pt_step4zWindow', 'Standard', 'Jets', 'p_{T}^{jet2}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet3_pt_step4zWindow', 'Standard', 'Jets', 'p_{T}^{jet3}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet4_pt_step4zWindow', 'Standard', 'Jets', 'p_{T}^{jet4}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet5_pt_step4zWindow', 'Standard', 'Jets', 'p_{T}^{jet5}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet6_pt_step4zWindow', 'Standard', 'Jets', 'p_{T}^{jet6}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet1_pt_step5', 'Standard', 'Jets', 'p_{T}^{jet1}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet2_pt_step5', 'Standard', 'Jets', 'p_{T}^{jet2}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet3_pt_step5', 'Standard', 'Jets', 'p_{T}^{jet3}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet4_pt_step5', 'Standard', 'Jets', 'p_{T}^{jet4}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet5_pt_step5', 'Standard', 'Jets', 'p_{T}^{jet5}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet6_pt_step5', 'Standard', 'Jets', 'p_{T}^{jet6}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet1_pt_step5zWindow', 'Standard', 'Jets', 'p_{T}^{jet1}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet2_pt_step5zWindow', 'Standard', 'Jets', 'p_{T}^{jet2}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet3_pt_step5zWindow', 'Standard', 'Jets', 'p_{T}^{jet3}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet4_pt_step5zWindow', 'Standard', 'Jets', 'p_{T}^{jet4}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet5_pt_step5zWindow', 'Standard', 'Jets', 'p_{T}^{jet5}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet6_pt_step5zWindow', 'Standard', 'Jets', 'p_{T}^{jet6}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet1_pt_step6', 'Standard', 'Jets', 'p_{T}^{jet1}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet2_pt_step6', 'Standard', 'Jets', 'p_{T}^{jet2}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet3_pt_step6', 'Standard', 'Jets', 'p_{T}^{jet3}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet4_pt_step6', 'Standard', 'Jets', 'p_{T}^{jet4}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet5_pt_step6', 'Standard', 'Jets', 'p_{T}^{jet5}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet6_pt_step6', 'Standard', 'Jets', 'p_{T}^{jet6}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet1_pt_step6zWindow', 'Standard', 'Jets', 'p_{T}^{jet1}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet2_pt_step6zWindow', 'Standard', 'Jets', 'p_{T}^{jet2}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet3_pt_step6zWindow', 'Standard', 'Jets', 'p_{T}^{jet3}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet4_pt_step6zWindow', 'Standard', 'Jets', 'p_{T}^{jet4}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet5_pt_step6zWindow', 'Standard', 'Jets', 'p_{T}^{jet5}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet6_pt_step6zWindow', 'Standard', 'Jets', 'p_{T}^{jet6}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet1_pt_step7', 'Standard', 'Jets', 'p_{T}^{jet1}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet2_pt_step7', 'Standard', 'Jets', 'p_{T}^{jet2}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet3_pt_step7', 'Standard', 'Jets', 'p_{T}^{jet3}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet4_pt_step7', 'Standard', 'Jets', 'p_{T}^{jet4}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet5_pt_step7', 'Standard', 'Jets', 'p_{T}^{jet5}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet6_pt_step7', 'Standard', 'Jets', 'p_{T}^{jet6}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet1_pt_step7zWindow', 'Standard', 'Jets', 'p_{T}^{jet1}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet2_pt_step7zWindow', 'Standard', 'Jets', 'p_{T}^{jet2}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet3_pt_step7zWindow', 'Standard', 'Jets', 'p_{T}^{jet3}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet4_pt_step7zWindow', 'Standard', 'Jets', 'p_{T}^{jet4}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet5_pt_step7zWindow', 'Standard', 'Jets', 'p_{T}^{jet5}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet6_pt_step7zWindow', 'Standard', 'Jets', 'p_{T}^{jet6}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet1_pt_step8', 'Standard', 'Jets', 'p_{T}^{jet1}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet2_pt_step8', 'Standard', 'Jets', 'p_{T}^{jet2}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet3_pt_step8', 'Standard', 'Jets', 'p_{T}^{jet3}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet4_pt_step8', 'Standard', 'Jets', 'p_{T}^{jet4}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet5_pt_step8', 'Standard', 'Jets', 'p_{T}^{jet5}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet6_pt_step8', 'Standard', 'Jets', 'p_{T}^{jet6}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet1_pt_step8zWindow', 'Standard', 'Jets', 'p_{T}^{jet1}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet2_pt_step8zWindow', 'Standard', 'Jets', 'p_{T}^{jet2}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet3_pt_step8zWindow', 'Standard', 'Jets', 'p_{T}^{jet3}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet4_pt_step8zWindow', 'Standard', 'Jets', 'p_{T}^{jet4}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet5_pt_step8zWindow', 'Standard', 'Jets', 'p_{T}^{jet5}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],
#['basic_jet6_pt_step8zWindow', 'Standard', 'Jets', 'p_{T}^{jet6}', '5', '0', '0', '1', '0.1', '0', '30', '299', '0', '0'],

#['basic_jet1_eta_step3', 'Standard', 'Jets', '#eta_{jet1}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet2_eta_step3', 'Standard', 'Jets', '#eta_{jet2}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet3_eta_step3', 'Standard', 'Jets', '#eta_{jet3}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet4_eta_step3', 'Standard', 'Jets', '#eta_{jet4}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet5_eta_step3', 'Standard', 'Jets', '#eta_{jet5}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet6_eta_step3', 'Standard', 'Jets', '#eta_{jet6}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet1_eta_step4', 'Standard', 'Jets', '#eta_{jet1}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet2_eta_step4', 'Standard', 'Jets', '#eta_{jet2}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet3_eta_step4', 'Standard', 'Jets', '#eta_{jet3}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet4_eta_step4', 'Standard', 'Jets', '#eta_{jet4}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet5_eta_step4', 'Standard', 'Jets', '#eta_{jet5}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet6_eta_step4', 'Standard', 'Jets', '#eta_{jet6}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet1_eta_step5', 'Standard', 'Jets', '#eta_{jet1}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet2_eta_step5', 'Standard', 'Jets', '#eta_{jet2}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet3_eta_step5', 'Standard', 'Jets', '#eta_{jet3}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet4_eta_step5', 'Standard', 'Jets', '#eta_{jet4}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet5_eta_step5', 'Standard', 'Jets', '#eta_{jet5}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet6_eta_step5', 'Standard', 'Jets', '#eta_{jet6}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet1_eta_step6', 'Standard', 'Jets', '#eta_{jet1}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet2_eta_step6', 'Standard', 'Jets', '#eta_{jet2}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet3_eta_step6', 'Standard', 'Jets', '#eta_{jet3}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet4_eta_step6', 'Standard', 'Jets', '#eta_{jet4}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet5_eta_step6', 'Standard', 'Jets', '#eta_{jet5}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet6_eta_step6', 'Standard', 'Jets', '#eta_{jet6}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet1_eta_step7', 'Standard', 'Jets', '#eta_{jet1}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet2_eta_step7', 'Standard', 'Jets', '#eta_{jet2}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet3_eta_step7', 'Standard', 'Jets', '#eta_{jet3}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet4_eta_step7', 'Standard', 'Jets', '#eta_{jet4}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet5_eta_step7', 'Standard', 'Jets', '#eta_{jet5}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet6_eta_step7', 'Standard', 'Jets', '#eta_{jet6}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet1_eta_step8', 'Standard', 'Jets', '#eta_{jet1}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet2_eta_step8', 'Standard', 'Jets', '#eta_{jet2}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet3_eta_step8', 'Standard', 'Jets', '#eta_{jet3}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet4_eta_step8', 'Standard', 'Jets', '#eta_{jet4}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet5_eta_step8', 'Standard', 'Jets', '#eta_{jet5}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],
#['basic_jet6_eta_step8', 'Standard', 'Jets', '#eta_{jet6}', '2', '0', '0', '1', '0.1', '0', '-2.4', '2.39', '0', '0'],

#['#', 'Name,', 'Extra,', 'axis', 'labels', '(y,x),', 'rebin,', 'do_dyscale,', 'logx,', 'logy,', 'ymin,', 'ymax,', 'xmin,', 'xmax,', 'nbins,', 'xbins,', 'bcs'],
#['###', 'b-tagged', 'jet', 'properties']

#['basic_bjet_multiplicity_step3', 'Standard', 'Events', 'N_{bjets}', '1', '0', '0', '1', '1e2', '0', '0', '9', '0', '0'],
#['basic_bjet_multiplicity_step4', 'Standard', 'Events', 'N_{bjets}', '1', '0', '0', '1', '1e2', '0', '0', '9', '0', '0'],
#['basic_bjet_multiplicity_step4zWindow', 'Standard', 'Events', 'N_{bjets}', '1', '0', '0', '1', '1e2', '0', '0', '9', '0', '0'],
#['basic_bjet_multiplicity_step5', 'Standard', 'Events', 'N_{bjets}', '1', '1', '0', '0', '0', '0', '0', '9', '0', '0'],
#['basic_bjet_multiplicity_step5zWindow', 'Standard', 'Events', 'N_{bjets}', '1', '1', '0', '0', '0', '0', '0', '9', '0', '0'],
#['basic_bjet_multiplicity_step6', 'Standard', 'Events', 'N_{bjets}', '1', '1', '0', '1', '1e2', '0', '0', '9', '0', '0'],
#['basic_bjet_multiplicity_step6zWindow', 'Standard', 'Events', 'N_{bjets}', '1', '1', '0', '1', '1e2', '0', '0', '9', '0', '0'],
#['basic_bjet_multiplicity_step7', 'Standard', 'Events', 'N_{bjets}', '1', '1', '0', '0', '0', '0', '0', '9', '0', '0'],
#['basic_bjet_multiplicity_step7zWindow', 'Standard', 'Events', 'N_{bjets}', '1', '1', '0', '0', '0', '0', '0', '9', '0', '0'],
#['basic_bjet_multiplicity_step8', 'Standard', 'Events', 'N_{bjets}', '1', '1', '0', '1', '1e2', '0', '0', '9', '0', '0'],
#['basic_bjet_multiplicity_step8zWindow', 'Standard', 'Events', 'N_{bjets}', '1', '1', '0', '1', '1e2', '0', '0', '9', '0', '0'],

#['basic_jet_btagDiscriminator_step6', 'Standard', 'Jets', '"b-tag.', 'discriminator"', '2', '1', '0', '0', '0', '0', '0', '0', '0', '0'],


#['#', 'Name,', 'Extra,', 'axis', 'labels', '(y,x),', 'rebin,', 'do_dyscale,', 'logx,', 'logy,', 'ymin,', 'ymax,', 'xmin,', 'xmax,', 'nbins,', 'xbins,', 'bcs'],
#['###', 'MET', 'properties']

#['basic_met_et_step3', 'Standard', 'Entries', 'E_{met}', '10', '0', '0', '0', '0', '0', '0', '299', '0', '0'],
#['basic_met_et_step4', 'Standard', 'Entries', 'E_{met}', '10', '0', '0', '0', '0', '0', '0', '299', '0', '0'],
#['basic_met_et_step4zWindow', 'Standard', 'Entries', 'E_{met}', '10', '0', '0', '0', '0', '0', '0', '299', '0', '0'],
#['basic_met_et_step5', 'Standard', 'Entries', 'E_{met}', '10', '1', '0', '0', '0', '0', '0', '299', '0', '0'],
#['basic_met_et_step5zWindow', 'Standard', 'Entries', 'E_{met}', '10', '1', '0', '0', '0', '0', '0', '299', '0', '0'],
#['basic_met_et_step6', 'Standard', 'Entries', 'E_{met}', '10', '1', '0', '0', '0', '0', '0', '299', '0', '0'],
#['basic_met_et_step6zWindow', 'Standard', 'Entries', 'E_{met}', '10', '1', '0', '0', '0', '0', '0', '299', '0', '0'],
#['basic_met_et_step7', 'Standard', 'Entries', 'E_{met}', '10', '1', '0', '0', '0', '0', '0', '299', '0', '0'],
#['basic_met_et_step7zWindow', 'Standard', 'Entries', 'E_{met}', '10', '1', '0', '0', '0', '0', '0', '299', '0', '0'],
#['basic_met_et_step8', 'Standard', 'Entries', 'E_{met}', '10', '1', '0', '0', '0', '0', '0', '299', '0', '0'],
#['basic_met_et_step8zWindow', 'Standard', 'Entries', 'E_{met}', '10', '1', '0', '0', '0', '0', '0', '299', '0', '0'],

#['basic_met_phi_step2', 'Standard', 'Entries', '\\phi_{met}', '2', '0', '0', '0', '0', '0', '0', '0', '0', '0'],
#['basic_met_phi_step3', 'Standard', 'Entries', '\\phi_{met}', '2', '0', '0', '0', '0', '4e6', '0', '0', '0', '0'],
#['basic_met_phi_step4', 'Standard', 'Entries', '\\phi_{met}', '2', '0', '0', '0', '0', '0', '0', '0', '0', '0'],
#['basic_met_phi_step5', 'Standard', 'Entries', '\\phi_{met}', '2', '1', '0', '0', '0', '0', '0', '0', '0', '0'],
#['basic_met_phi_step6', 'Standard', 'Entries', '\\phi_{met}', '5', '1', '0', '0', '0', '80e3', '0', '0', '0', '0'],
#['basic_met_phi_step7', 'Standard', 'Entries', '\\phi_{met}', '5', '1', '0', '0', '0', '65e3', '0', '0', '0', '0'],
#['basic_met_phi_step7L', 'Standard', 'Entries', '\\phi_{met}', '5', '1', '0', '0', '0', '65e3', '0', '0', '0', '0'],
#['basic_met_phi_step8', 'Standard', 'Entries', '\\phi_{met}', '5', '1', '0', '0', '0', '60e3', '0', '0', '0', '0'],

#['basic_met_E_step2', 'Standard', 'Jets', 'p_{t}^{met}', '5', '0', '0', '1', '1', '0', '0', '299', '0', '0'],
#['basic_met_E_step3', 'Standard', 'Jets', 'p_{t}^{met}', '5', '0', '0', '1', '1', '0', '0', '299', '0', '0'],
#['basic_met_E_step4', 'Standard', 'Jets', 'p_{t}^{met}', '5', '0', '0', '1', '0.1', '0', '0', '299', '0', '0'],
#['basic_met_E_step5', 'Standard', 'Jets', 'p_{t}^{met}', '5', '0', '0', '1', '0.1', '0', '0', '299', '0', '0'],
#['basic_met_E_step6', 'Standard', 'Jets', 'p_{t}^{met}', '5', '0', '0', '1', '0.1', '0', '0', '299', '0', '0'],
#['basic_met_E_step7', 'Standard', 'Jets', 'p_{t}^{met}', '5', '0', '0', '1', '0.1', '0', '0', '299', '0', '0'],

#['###', 'Additional', 'Control', 'Plots']

#['#', 'Name,', 'Extra,', 'axis', 'labels', '(y,x),', 'rebin,', 'do_dyscale,', 'logx,', 'logy,', 'ymin,', 'ymax,', 'xmin,', 'xmax,', 'nbins,', 'xbins,', 'bcs'],

]

for era in eras:
    print(era)
    
    for channel in channels:
        print(channel)

        runNominalfile = open("scripts/LoadPlotter_Nominal_" + era + ".sh","a") 
        runSystsfile = open("scripts/LoadPlotter_Systs_" + era + ".sh","a") 
        runEnvelopesfile = open("scripts/LoadPlotter_Envelopes_" + era + ".sh","a") 
        runBandsfile = open("scripts/LoadPlotter_Bands_" + era + ".sh","a") 

        file = open("data/binnings/HistoList_CP_" + era + "_" + channel,"a") 

        file.write('# Name, Extra, axis labels (y,x), rebin, do_dyscale, logx, logy, ymin, ymax, xmin, xmax, nbins, xbins, bcs\n')
        file.write('\n')

        for cp in controlplototherlist:
            #print(cp)
            rootfilepath = 'Plots_' + era + '/Nominal/' + channel + '/'+cp[0]+'_source.root'
            if os.path.exists(rootfilepath):
                rootfile = ROOT.TFile.Open(rootfilepath, "read")
                hist = rootfile.Get(cp[0]+'_data')
                ymax = 0
                for j in range(1,hist.GetNbinsX()-1):
                    if hist.GetBinContent(j) > ymax:
                        ymax = hist.GetBinContent(j)
                #print(ymax)
            else:
                ymax = float(cp[9])

            if float(cp[7]) == 1 and ymax != 0:
                ymax = math.pow(10,3+roundup(math.log10(ymax * 1.5)))
                #print(ymax)
            else:
                ymax = roundup(ymax * 1.5,5000)
                #print (ymax)
            file.write(cp[0]+' '+cp[1]+' '+cp[2]+' '+cp[3]+' '+cp[4]+' '+cp[5]+' '+cp[6]+' '+cp[7]+' '+cp[8]+' '+str(ymax)+' '+cp[10]+' '+cp[11]+' '+cp[12]+' '+cp[13]+'\n')

            if(channel == "combined"):
                runNominalfile.write('nohup ./install/bin/load_Plotter -y '+era+' -t cp -s Nominal -p +'+cp[0]+' &> plotslogs/'+cp[0]+'_'+era+'.out \n')
                runSystsfile.write('nohup ./install/bin/load_Plotter -y '+era+' -t cp -s all -p +'+cp[0]+' &> plotslogs/'+cp[0]+'_'+era+'.out \n')
                runEnvelopesfile.write('nohup ./install/bin/load_Plotter -y '+era+' -t cp -m -s all -p +'+cp[0]+' &> plotslogs/'+cp[0]+'_'+era+'.out \n')
                runBandsfile.write('nohup ./install/bin/load_Plotter -y '+era+' -t cp -b -s Nominal -p +'+cp[0]+' &> plotslogs/'+cp[0]+'_'+era+'.out \n')
            else: continue

        file.write('\n')

        file.close() 
        runNominalfile.close()
        runSystsfile.close()
        runEnvelopesfile.close()
        runBandsfile.close()
