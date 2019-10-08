# Grapevine-inflorescence-architecture
This code contains all Matlab functions needed to calculate persistence barcodes, 
bottleneck distances for persistent homology traits, simulation for berry potential, 
and other geometric features used in the study described in:

Mao Li, Laura L. Klein, Keith Duncan, Ni Jiang, Daniel H. Chitwood, Jason Londo, Allison J. Miller, Christopher N. Topp (2019)
Characterizing 3D inflorescence architecture in grapevine using X-ray imaging and advanced morphometrics: implications for understanding cluster density. bioRxiv 557819 and Journal of Experimental Botany, in press

Code is implemented by Mao Li (maoli0923@gmail.com, github:maoli0923).

A few remarks:

(a) The study was based on mesh.ply files from X-ray imaging.

(b) Computing persistence barcode needs to install Javaplex http://appliedtopology.github.io/javaplex/

Note that some functions such as 'edgeL2adj','networkComponents',and 'ply_read' are origially implemented by others. You can find these authors' names in each function file if available.
