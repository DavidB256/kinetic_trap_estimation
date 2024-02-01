## Kinetic trap estimation
Rotation project of David Bass in Margaret Johnson lab, Nov. 2023--Jan. 2024. Based on work initiated in Jhaveri and Loggia et al.'s in-review preprint.

### Contributions
I wrote a Julia framework for constructing ODE models of fully connected rate-growth chemical reaction networks of arbitrary size.
- Generalizing this framework to have other topologies or to make different assumptions would be non-trivial, but not difficult.

I wrote a Julia framework for estimating the rate constants that parameterize a 
- This framework needs a moderate amount of refactoring and documentation to be reproducible.
- Substituting the loss (aka objective) function, e.g. to change the type of data to which  is easy

I wrote a benchmarking wrapper for a routine that randomly samples true rate constants from whithin a narrow range, simulates experimental data, and then times repeated attempts at  simulates rate constant estimation.

### Potential next steps
- Re-implement more of the Kinetic_AssemblyAD repository in Julia.
- 

### Links
- (Original version of `rotation_talk.pdf` on Google Slides)[https://docs.google.com/presentation/d/1-Mb23PiFgSqMkGJpZ9rYqLjJCuRr13WviM6cUggBn4I/edit?usp=sharing]

### References (currently just a dump of my Zotero collection)
Baydin, A. G., Pearlmutter, B. A., Radul, A. A., & Siskind, J. M. (2018). Automatic differentiation in machine learning: A survey (arXiv:1502.05767). arXiv. https://doi.org/10.48550/arXiv.1502.05767

Faeder, J. R., Blinov, M. L., & Hlavacek, W. S. (2009). Rule-Based Modeling of Biochemical Systems with BioNetGen. In I. V. Maly (Ed.), Systems Biology (pp. 113–167). Humana Press. https://doi.org/10.1007/978-1-59745-525-1_5

Gartner, F. M., Graf, I. R., & Frey, E. (2022). The time complexity of self-assembly. Proceedings of the National Academy of Sciences, 119(4), e2116373119. https://doi.org/10.1073/pnas.2116373119

Hagan, M. F., Elrad, O. M., & Jack, R. L. (2011). Mechanisms of kinetic trapping in self-assembly and phase transformation. The Journal of Chemical Physics, 135(10), 104115. https://doi.org/10.1063/1.3635775

Jhaveri, A., Loggia, S., Qian, Y., & Johnson, M. E. (2023). Discovering optimal kinetic pathways for self-assembly using automatic differentiation. bioRxiv, 2023.08.30.555551. https://doi.org/10.1101/2023.08.30.555551

Kingma, D. P., & Ba, J. (2017). Adam: A Method for Stochastic Optimization (arXiv:1412.6980). arXiv. https://doi.org/10.48550/arXiv.1412.6980

Qian, Y., Evans, D., Mishra, B., Fu, Y., Liu, Z. H., Guo, S., & Johnson, M. E. (2023a). Temporal control by co-factors prevents kinetic trapping in retroviral Gag lattice assembly (p. 2023.02.08.527704). bioRxiv. https://doi.org/10.1101/2023.02.08.527704

Sánchez, A. D., Nicola, E. M., & Wio, H. S. (1997). Kinetics of Trapping Reactions with a Time Dependent Density of Traps. Physical Review Letters, 78(11), 2244–2247. https://doi.org/10.1103/PhysRevLett.78.2244
Varga, M. J., Fu, Y., Loggia, S., Yogurtcu, O. N., & Johnson, M. E. (2020). NERDSS: A Nonequilibrium Simulator for Multibody Self-Assembly at the Cellular Scale. Biophysical Journal, 118(12), 3026–3040. https://doi.org/10.1016/j.bpj.2020.05.002