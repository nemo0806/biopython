# -*- coding: utf-8 -*-
"""
Created on Sun Sep  9 21:46:14 2018

@author: tjk
"""

from Bio.Align.Applications import ClustalwCommandline

#clustalw_cline = ClustalwCommandline(r"E:\PC_Tools\Anaconda3\pkgs\biopython-1.72-py36h830ac7b_0\Lib\site-packages\Bio\Align\Applications\_Clustalw.py", infile="DNAseq\opuntia.fasta")
#stdout, stderr = clustalw_cline()
#print (clustalw_cline.value())

import os
clustalw_exe = r"C:\Program Files (x86)\ClustalW2\clustalw2.exe"
clustalw_cline = ClustalwCommandline(clustalw_exe, infile = r"C:\Users\tjk\Downloads\influenza.fasta")
assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
stdout, stderr = clustalw_cline()

from Bio import Phylo
tree = Phylo.read(r"C:\Users\tjk\Downloads\influenza.dnd", "newick")
tree.rooted = True
Phylo.draw(tree)
