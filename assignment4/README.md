# Assignment 4

Name: Haitong Shi

Legi-Nr: 20-960-340

## Required results
Edit this 'README.md' file to report all your results. You only need to update the tables in the reports section by adding screenshots and reporting results.

### Mandatory Tasks

1) Screenshots of the parameterizations and textured (checkerboard) models for all the implemented methods and boundary conditions (models: cathead.obj, hemisphere.off, hemisphere_non_convex_boundary.off, Octo_cut2.obj)

2) Several examples of the distortion visualizations.


## Reports
### (mandatory-1) parameterization and checkerboard texture models
#### cathead
| Method            | checkerboard textured models          |         Parameterization             |
| :--------------:  | ------------------------------------- |------------------------------------- |
| Uniform (fixed)   |<img align="center" src="./res/catuni.png" width="300">| <img align="center"  src="./res/catunip.png" width="300"> |
| Cotangent (fixed) |<img align="center" src="./res/catcot.png" width="300">| <img align="center"  src="./res/catcotp.png" width="300"> |
| LSCM (fixed)      |<img align="center" src="./res/catlscm.png" width="300">| <img align="center"  src="./res/catlscmp.png" width="300"> |
| ARAP (fixed)      |<img align="center" src="./res/catasap.png" width="300">| <img align="center"  src="./res/catasapp.png" width="300"> |
| LSCM (free)       |<img align="center" src="./res/catlscmf.png" width="300">| <img align="center"  src="./res/catlscmfp.png" width="300"> |
| ARAP (free)       |<img align="center" src="./res/catarapf.png" width="300">| <img align="center"  src="./res/catarapfp.png" width="300"> |

#### hemisphere
| Method            | checkerboard textured models          |         Parameterization             |
| :--------------:  | ------------------------------------- |------------------------------------- |
| Uniform (fixed)   |<img align="center" src="./res/hemiuni.png" width="300">| <img align="center"  src="./res/hemiunip.png" width="300"> |
| Cotangent (fixed) |<img align="center" src="./res/hemicot.png" width="300">| <img align="center"  src="./res/hemicotp.png" width="300"> |
| LSCM (fixed)      |<img align="center" src="./res/hemilscm.png" width="300">| <img align="center"  src="./res/hemilscmp.png" width="300"> |
| ARAP (fixed)      |<img align="center" src="./res/hemiarap.png" width="300">| <img align="center"  src="./res/hemiarapp.png" width="300"> |
| LSCM (free)       |<img align="center" src="./res/hemilscmf.png" width="300">| <img align="center"  src="./res/hemilscmfp.png" width="300"> |
| ARAP (free)       |<img align="center" src="./res/hemiarapf.png" width="300">| <img align="center"  src="./res/hemiarapfp.png" width="300"> |

#### hemisphere_non_convex_boundary
| Method            | checkerboard textured models          |         Parameterization             |
| :--------------:  | ------------------------------------- |------------------------------------- |
| Uniform (fixed)   |<img align="center" src="./res/heminuni.png" width="300">| <img align="center"  src="./res/heminunip.png" width="300"> |
| Cotangent (fixed) |<img align="center" src="./res/hemincot.png" width="300">| <img align="center"  src="./res/hemincotp.png" width="300"> |
| LSCM (fixed)      |<img align="center" src="./res/heminlscm.png" width="300">| <img align="center"  src="./res/heminlscmp.png" width="300"> |
| ARAP (fixed)      |<img align="center" src="./res/heminarap.png" width="300">| <img align="center"  src="./res/heminarapp.png" width="300"> |
| LSCM (free)       |<img align="center" src="./res/heminlscmf.png" width="300">| <img align="center"  src="./res/heminlscmfp.png" width="300"> |
| ARAP (free)       |<img align="center" src="./res/heminarapf.png" width="300">| <img align="center"  src="./res/heminarapfp.png" width="300"> |

#### Octo_cut2
| Method            | checkerboard textured models          |         Parameterization             |
| :--------------:  | ------------------------------------- |------------------------------------- |
| Uniform (fixed)   |<img align="center" src="./res/octouni.png" width="300">| <img align="center"  src="./res/octounip.png" width="300"> |
| Cotangent (fixed) |<img align="center" src="./res/octocot.png" width="300">| <img align="center"  src="./res/octocotp.png" width="300"> |
| LSCM (fixed)      |<img align="center" src="./res/octolscm.png" width="300">| <img align="center"  src="./res/octolscmp.png" width="300"> |
| ARAP (fixed)      |<img align="center" src="./res/octoarap.png" width="300">| <img align="center"  src="./res/octoarapp.png" width="300"> |
| LSCM (free)       |<img align="center" src="./res/octolscmf.png" width="300">| <img align="center"  src="./res/octolscmfp.png" width="300"> |
| ARAP (free)       |<img align="center" src="./res/octoarapfp.png" width="300">| <img align="center"  src="./res/octoarapf.png" width="300"> |


### (mandatory-2) distortion visualization
#### cathead
| mtd \ metric      | Conformal (angle) |    Authalic (area)  |  Isometric  (length)    |
| :--------------:  | ----------------- | ------------------- | ----------------------- |
| LSCM (free)       |<img align="center" src="./res/cat_lscm_1.png" width="300">| <img align="center"  src="./res/cat_lscm_2.png" width="300"> | <img align="center"  src="./res/cat_lscm_3.png" width="300"> |
| ARAP (free) |<img align="center" src="./res/cat_arap_1.png" width="300">| <img align="center"  src="./res/cat_arap_2.png" width="300"> |<img align="center"  src="./res/cat_arap_3.png" width="300"> |


#### hemisphere
| mtd \ metric      | Conformal (angle) |    Authalic (area)  |  Isometric  (length)    |
| :--------------:  | ----------------- | ------------------- | ----------------------- |
| LSCM (free)       |<img align="center" src="./res/hemi_lscm_1.png" width="300">| <img align="center"  src="./res/hemi_lscm_2.png" width="300"> | <img align="center"  src="./res/hemi_lscm_3.png" width="300"> |
| ARAP (free) |<img align="center" src="./res/hemi_arap_1.png" width="300">| <img align="center"  src="./res/hemi_arap_2.png" width="300"> |<img align="center"  src="./res/hemi_arap_3.png" width="300"> |


