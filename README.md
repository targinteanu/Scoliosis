# Scoliosis
BE 490 Independent Study - Quantifying Scoliosis Severity using Differential Geometry

The purpose of this project is to investigate quantities such as Writhe, Twist, and Torsion as potential metrics of scoliosis severity. 

files include a graphical user interface for segmenting vertebrae from axial CT scan and a script for analyzing segmented vertebrae to calculate twist and writhe. 

getTwist, getWrithe, getTosrion use discrete approximations of integrals/derivatives. levittWrithe uses a line-segment-chain approach to calculate Writhe that has previously been used on protein and DNA folding. lewinerTorsion uses a least-squares approximation to estimate Torsion that has previously been used on scoliotic spines. deturckTwist and deturckWrithe use vertebral endplates to estimate tangent vectors. 
