# Table of Contents

1.  [Presentation](#org81fa71d)
2.  [Installation](#org70fca4a)
3.  [Reduce Runtime](#org63fbd3e)
3.  [Remarks](#org23er53e)

<a id="org81fa71d"></a>

# Presentation

We simulate a signal travelling through a synaptic cleft. Neurotransmitters are released from a axon terminal and are bounded to 
receptors at a dendrite spine, located at the other end of the synaptic cleft. We look at 3D and a geometrical reduction to 2D. 
We introduce glia cells that contain transporters that react with neurotransmitters, to clear the synaptic cleft. We introduce the 
possible of cross-talk, where neurotransmitters can diffuse from one synaptic cleft to a neighbouring one. Finally we include an 
underlying flow in the fluid in the synaptic cleft.

The neurotransmitters are spatially resolved and can diffuse. The receptors are fixed at the dendrite spine.
The transporters are fixed at the glia cells.

The simulator setup is based on BattMo developed at SINTEF which provides generic tools for simulation electro
chemical system.

Results can be found in ScaledGifs folder.

<a id="org70fca4a"></a>

# Installation

-   Install [BattMo](https://github.com/BattMoTeam/BattMo). Follow the [installation instruction](https://github.com/BattMoTeam/BattMo#installation) from the webpage
-   Add the directories `Scaled`, 'Unscaled' and `Models` in your MATLAB path (or run `startup` script)
    
<a id="org63fbd3e"></a>
# Reduce runtime
All files has a nx, ny (nz in 3D) specifying number of gridpoints in each direction. These can be reduced to reduce the runtime.
All files also has a n, specifying the number of timesteps to be taken. This can also be reduced.

<a id="org23er53e"></a>
# Remarks
It should be noted that although GitHub says that 99.9% of the code was written in Jupyter Notebook, the entirety of the new attempt is written in Matlab. This percentage is because some of the jupyter notebook files contain large animations.
