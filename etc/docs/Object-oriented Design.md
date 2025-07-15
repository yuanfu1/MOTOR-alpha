# Object-oriented Design

This document focuses on the module designs.

## From Control variables to Observations
MOTOR-DA separates the control variables and observations by two abstract modules, C2M and M2O.

	C2M: this module transforms the control variables on to model states;

	M2O: this module transforms the model states into the requested observation, observation operators;

These setting seems providing flexibility for various applications, 3DVAR or even 4DVAR.

For examples:

1. For cumulus parameterization scheme, C2M implements an instantiation maps Q and T etc to model states. M2O does not even know it runs a cumulus parameterization or not.
2. For a 4DVAR, forecast models map initial and/or boundary conditions on to model states.

C2O implementation:
It constructs an abstract module, and implements different controls. This has not been implemented. We added all observation operators to this module and make it tedious and hard to understand.

For example of the wind speed and direction as control case, it implements a subclasse mapping wind speeds and directions on to U and V.

M2O implementation:
It constructs an abstract module and implements various observations in its subclasses. This has been implemented in the current MOTOR-DA.

For example of satellite data, it implements a mapping using RTTOV from model states to brightness temporature.

The minimization algorithm simply uses controls along with their gradients to solve.

	M2O: we need to define a unified model states so that any changes in the C2M would not affect the M2O. 
	For now, we set the following as the model states: $u$, $v$, $p$, $T$, $q$ (water vapor mixing ratio).

For a vorticity and divergence control variable case, C2O_Z module converts the controls into $u$ and $v$ by solving Poisson equations for streamfunction and velocity potential and then into $u$ and $v$. Note that for vorticity and divergence control variable case, the $u$ and $v¥ around the boundaries are control variables too. Thus, the solver of streamfunction and velocity potential is a mixed $\psi$ and $\chi$ boundary problem. We have already developed such a solver implemented by a multigrid.

```
MOTOR-DA
│
├──	 BkgdMD				% Multigrid background processing module
│
├──	 BMartix			% Background error covariance matrix
│
├──	 C2M				% define the controls to predefined model states
│
├──  M2O				% From predefined model states transfer to observation space
│
├──	 RMartix			% Observation error covariance
│
└──  States  			% Model state data structure

```

The C2M module closely resembles Ctl2State in its functionality. However, while Ctl2State consolidates all options into a single module, resulting in a cumbersome and complex structure that is difficult to maintain, C2M takes a more abstract approach. By allowing users to define their own controls and craft the forward, tangent linear, and adjoint methods according to their preferences, C2M offers a level of flexibility that prevents options from conflicting with each other. This design not only streamlines the module but also enhances its adaptability and ease of maintenance.

```
C2M offers a default control variable setting and users can inherit it to construct their own control case. 
	The default setting is:
	u,v,w,lnP,"rh" a exponential scaled water vapor density. Users can choose any profile to scale model 
	mixing raito into "rh".
```

### The current implementation of control variables spreads all of the places, including MGOpt.F90 e.g. for satellite. This makes our system hard to maintain and update, even debugging. I tried to use the C2M module separate all of the issues related to control variable transformation. 
#### The current setting: water vapor scalings are set in singleGrid; they are included in the observation error covariance. It mixes the control transformation with the covariance matrix. The new setting will remove the scaling into control variables explicitly and keep the error covariance alone.

#### Another improvement of this new structure helps to remove the Ctl2State, which is applied at different places. It is hard to understand and to maintain the system. Particularly, there are some hard coded variable transformations.