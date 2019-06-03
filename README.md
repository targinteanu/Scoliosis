# Scoliosis
BE 490 Independent Study - Quantifying Scoliosis Severity using Differential Geometry

The purpose of this project is to investigate quantities such as Writhe, Twist, and Torsion as potential metrics of scoliosis severity. 

Scripts metricsFromXL and TorsionFromXL were used to calculate metrics from a spreadsheet containing XYZ/rotation data and previous clustering information. GUI1 was used to segment axial CT images into XYZ/rotation data. Functions levittWrithe, getTwist, lewinerTorsion, and levittWritheAbs were used to get metrics from XYZ/rotation data. 

getTwist, getWrithe, getTosrion use discrete approximations of integrals/derivatives. levittWrithe uses a line-segment-chain approach to calculate Writhe that has previously been used on protein and DNA folding. lewinerTorsion uses a least-squares approximation to estimate Torsion that has previously been used on scoliotic spines. deturckTwist and deturckWrithe use vertebral endplates to estimate tangent vectors. 
