# SPECFEM-C 

This form of SPECFEM-X is for the development of SPECFEM-C prior to any integration of SPECFEM-X with SPECFEM Kokkos. It incorporates the ability to model sea-level change based on the rate-dependent formulation of [Crawford, et al., 2018](https://scholar.google.com/scholar?hl=en&as_sdt=0%2C31&q=al+attar+crawford+2018&btnG=&oq=craw). 


## To do list: 

- [ ] Add the final weak form term (ice contibution)
- [x] Write ocean function updater using ice and theta comparison 
- [ ] When calling ```set_petsc_stiffness_SL``` do we need to parse the SL matrices as args? arent they global?

### Time looping
- [ ] Stiffness matrix needs to be updated at each timestep...not just the first 1 (elastic) or 2 (viscoelastic)
- [ ] Need to store ```nodalu```, ```nodalphi```, ```nodalsl``` etc at each timestep (not to be overwritten) 
- [ ] Incorporate ice load into the RHS 
- [ ] Add in time-marching scheme to estimate  ${\phi_{t+1}}$, $m_{t+1}$,  $u_{t+1}$, and $\theta_{t+1}$,


## Clarify with HNG 
- [ ] **Do we need to incorporate SL terms into the true Kmat so that it influences the ```bcnodalv``` calculations?**
- [ ] ```apply_nonzero_bc``` - double check exactly what is happening here 
- [ ] Why do we keep setting lots of the loads' first elements to zero. e.g. ```bodyload(0) = ZERO ```
- [ ] Do we need to split up the memory integral $$\int_{M_s} 2\mu_0\Big[\dot{\mathbf{m}} : \tilde{\mathbf{m}} + \frac{1}{\tau} (\mathbf{d} - \mathbf{m}):(\tilde{\mathbf{d}} - \tilde{\mathbf{m}}) \Big]  dV$$ so that $$\int_{M_s} 2\mu_0\Big[ \frac{1}{\tau} (\mathbf{d} - \mathbf{m}):(\tilde{\mathbf{d}} - \tilde{\mathbf{m}}) \Big]  dV$$ is part of the RHS?
- [ ] What are the different loads being used? 
